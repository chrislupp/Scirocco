#!/usr/bin/env python3
"""Quick: solve mfoil at the comprehensive_reference base state and append
the assembled R_V matrix as /comprehensive/R_V. Avoids re-running the full
202-perturbation CS shape extraction.
"""
import os
import sys
import h5py
import numpy as np
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'mfoil'))
import mfoil  # noqa: E402

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'comprehensive_reference.h5')

# Same setup as extract_comprehensive_reference.py
M = mfoil.mfoil(naca='0012', npanel=100)
M.param.verb = 0
M.param.doplot = False
M.setoper(alpha=3.0, Re=1e6, Ma=0.0)


def recompute_geom(M):
    x = M.foil.x
    N = M.foil.N
    ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
    s = np.zeros(N, dtype=x.dtype)
    s[1:] = np.cumsum(ds)
    M.foil.s = s


# Match comprehensive_reference solve_base() exactly: ue-based xi
# convention (after commit 4051c2d): converge at spline foil.s, switch
# to cumulative, refresh sstag/xi via stagpoint_move at the new foil.s.
M.solve()
recompute_geom(M)
mfoil.build_gamma(M, 3.0)
mfoil.build_wake(M)
mfoil.identify_surfaces(M)
mfoil.set_wake_gap(M)
mfoil.calc_ue_m(M)
mfoil.stagpoint_move(M)
mfoil.build_glob_sys(M)

Nsys = M.glob.Nsys
NN = 4 * Nsys
ue = M.glob.U[3, :]
ds = M.glob.U[1, :]
Ids = slice(1, 4 * Nsys, 4)
Iue = slice(3, 4 * Nsys, 4)

R_V = sparse.lil_matrix((NN, NN))
R_V[0:3 * Nsys, 0:4 * Nsys] = M.glob.R_U
I = slice(3 * Nsys, 4 * Nsys, 1)
R_V[I, Iue] = sparse.identity(Nsys) - M.vsol.ue_m @ np.diag(ds)
R_V[I, Ids] = -M.vsol.ue_m @ np.diag(ue)
R_V = R_V.tocsr().toarray()

with h5py.File(OUTPUT_FILE, 'r+') as f:
    g = f['comprehensive']
    if 'R_V' in g:
        del g['R_V']
    g.create_dataset('R_V', data=R_V)
    # IMPORTANT: also save sstag_ue. mfoil's sstag_ue at this Frankenstein
    # state (Newton-converged U, recompute_geom to cumulative foil.s,
    # stagpoint_find to gamma-based xi) is the SPLINE-foil.s value set during
    # the initial M.solve(). It's NOT recomputed by recompute_geom or
    # stagpoint_find. cppfoil tests need this exact value to match mfoil's
    # R_V at the xi-stag-chain augmentation.
    if 'sstag_ue' in g:
        del g['sstag_ue']
    g.create_dataset('sstag_ue', data=np.array(M.isol.sstag_ue))
    if 'sstag_g' in g:
        del g['sstag_g']
    g.create_dataset('sstag_g', data=np.array(M.isol.sstag_g))
    print(f"cl_base saved: {float(g['cl_base'][()])}")
    print(f"cl computed:   {float(M.post.cl)}")
    print(f"sstag_ue saved: {np.array(M.isol.sstag_ue)}")
print(f"Wrote R_V {R_V.shape} to {OUTPUT_FILE}")
