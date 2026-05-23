#!/usr/bin/env python3
"""
Generate machine-precision reference data for the Phase 1 viscous geometric
chain functions.

Outputs: cppfoil/tests/geom_chain_reference.h5

HDF5 structure (group "geom_chain"):

  Base state (cumulative-panel-norm foil.s, matching cppfoil RecomputeArcLength):
    foil_x, foil_s, foil_gam, wake_x, wake_s, sstag, Istag, xi, wgap,
    chord, wakelen, alpha

  CS-derived chain sensitivities (foil.x perturbed complex, foil.s recomputed
  cumulative, full inviscid solve + wake build + set_wake_gap rerun):
    dfoil_s_dx [N, 2N]
    dsstag_dx [2N]
    dwake_s_dx [Nw, 2N]
    dwake_x_dx [2*Nw, 2N]   (row 2*j+d = wake.x(d, j))
    dwgap_dx [Nw, 2N]
    dxi_dx [N+Nw, 2N]

Test case: NACA 0012, 100 panels, alpha=3 (inviscid solve only — no Re).
"""

import sys
import os
import copy
import time
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', 'temp'))
import mfoil_cs  # noqa: F401, E402 — installs monkey patches into mfoil
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil  # noqa: E402

CS_STEP = 1e-30
NACA = "0012"
NPANEL = 100
ALPHA = 3.0

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'geom_chain_reference.h5')


def recompute_geom(M):
    """Cumulative-panel-norm foil.s (matches cppfoil RecomputeArcLength).

    foil.x is read as-is (may be complex); foil.s is written with matching
    dtype.
    """
    x = M.foil.x
    N = M.foil.N
    ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
    s = np.zeros(N, dtype=x.dtype)
    s[1:] = np.cumsum(ds)
    M.foil.s = s


def rebuild_geom_chain(M, alpha):
    """Re-run the geometric chain at the current foil.x (and updated foil.s):
       solve_inviscid -> build_wake -> stagpoint_find -> identify_surfaces ->
       set_wake_gap.

    Used by the BASE solve (real-valued) to establish a consistent starting
    state. Inside the CS loop we use rebuild_geom_chain_fixedU instead, which
    holds M.glob.U[3, :] fixed and runs stagpoint_move (ue-based) — matching
    cppfoil's production convention after commit 4051c2d."""
    mfoil.solve_inviscid(M)   # rebuilds AIC, gamref, gam, stagpoint_find (xi foil-only)
    mfoil.build_wake(M)
    mfoil.stagpoint_find(M)
    mfoil.identify_surfaces(M)
    mfoil.set_wake_gap(M)


def rebuild_geom_chain_fixedU(M, alpha):
    """Hybrid convention matching cppfoil's production:
      - sstag/xi: ue-based, with M.glob.U[3, :] HELD FIXED (matches
        cppfoil's DSstagDFoilX reading glob.U(3, Istag) as fixed inputs).
      - wake.x: includes the chain through gam → uewi (matches cppfoil's
        DWakeXDFoilX which chains through DGammaDFoilX).

    solve_inviscid updates gam/uewi but does NOT touch M.glob.U, so the
    subsequent stagpoint_move still sees the fixed U[3]."""
    mfoil.solve_inviscid(M)  # updates gam, AIC, gamref, uewi (does not touch glob.U)
    mfoil.build_wake(M)      # uses new uewi → new wake.x (gam chain captured)
    mfoil.stagpoint_move(M)  # ue-based, uses M.glob.U[3] which we hold fixed
    mfoil.identify_surfaces(M)
    mfoil.set_wake_gap(M)


def populate_glob_U_from_inviscid(M):
    """After solve_inviscid + build_wake + stagpoint_find, populate
    M.glob.U[3, :] with the inviscid ue values. Other rows zeroed —
    fixed-U CS only needs ue."""
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw
    M.glob.U = np.zeros((4, Nsys), dtype=M.foil.x.dtype)
    cs = np.array([np.cos(np.deg2rad(M.oper.alpha)),
                   np.sin(np.deg2rad(M.oper.alpha))],
                  dtype=M.foil.x.dtype)
    uea = M.isol.sgnue * (M.isol.gamref @ cs)
    M.glob.U[3, :N] = uea
    if Nw > 0:
        uew = M.isol.uewiref @ cs
        # First wake node = TE continuity
        uew[0] = uea[-1]
        M.glob.U[3, N:] = uew


def base_solve():
    """Build the base mfoil object, then run the full geometric chain at real
    foil.x with cumulative-panel-norm foil.s (consistent base for CS).

    Also populates M.glob.U[3, :] with inviscid ue (held fixed during CS).
    """
    M = mfoil.mfoil(naca=NACA, npanel=NPANEL)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=ALPHA)
    # Initial solve to populate things and chord/wakelen.
    M.solve()
    # Switch foil.s to cumulative-panel-norm so the base matches the
    # perturbation convention (recompute_geom inside the CS loop).
    recompute_geom(M)
    rebuild_geom_chain(M, ALPHA)
    populate_glob_U_from_inviscid(M)
    return M


def cs_perturb(M_base, k, d, step):
    """Deep-copy M_base, perturb foil.x[d, k] complex, recompute foil.s
    (cumulative), re-run the geom chain at FIXED U. Returns the perturbed
    mfoil object."""
    M = copy.deepcopy(M_base)
    # Promote required arrays to complex (M.glob.U held fixed but needs
    # complex dtype to participate in complex_step_mode arithmetic).
    M.foil.x = M.foil.x.astype(np.complex128)
    M.foil.s = M.foil.s.astype(np.complex128)
    M.glob.U = M.glob.U.astype(np.complex128)
    if isinstance(M.isol.gam, np.ndarray):
        M.isol.gam = M.isol.gam.astype(np.complex128)
    if isinstance(M.isol.gamref, np.ndarray):
        M.isol.gamref = M.isol.gamref.astype(np.complex128)
    if isinstance(M.isol.uewiref, np.ndarray):
        M.isol.uewiref = M.isol.uewiref.astype(np.complex128)
    M.foil.x[d, k] += 1j * step
    recompute_geom(M)
    with complex_step_mode():
        rebuild_geom_chain_fixedU(M, ALPHA)
    return M


def main():
    print(f"Generating Phase 1 geom chain reference: NACA {NACA}, "
          f"npanel={NPANEL}, alpha={ALPHA}")
    t0 = time.time()

    M_base = base_solve()
    N = M_base.foil.N
    Nw = M_base.wake.N
    Nsys = N + Nw
    twoN = 2 * N

    print(f"  N={N}, Nw={Nw}, Nsys={Nsys}")
    assert Nw > 0, "Wake build failed in base solve"

    # Capture base state
    foil_x = np.array(M_base.foil.x, dtype=float)
    foil_s = np.array(M_base.foil.s, dtype=float)
    foil_gam = np.array(M_base.isol.gam, dtype=float)
    wake_x = np.array(M_base.wake.x, dtype=float)
    wake_s = np.array(M_base.wake.s, dtype=float)
    sstag = float(M_base.isol.sstag)
    Istag = np.array(M_base.isol.Istag, dtype=np.int64)
    xi = np.array(M_base.isol.xi, dtype=float)
    wgap = np.array(M_base.vsol.wgap, dtype=float)
    chord = float(M_base.geom.chord)
    wakelen = float(M_base.geom.wakelen)

    # Allocate CS arrays
    dfoil_s = np.zeros((N, twoN))
    dsstag = np.zeros(twoN)
    dwake_s = np.zeros((Nw, twoN))
    dwake_x = np.zeros((2, Nw, twoN))
    dwgap = np.zeros((Nw, twoN))
    dxi = np.zeros((Nsys, twoN))

    for k in range(N):
        for d in range(2):
            M_p = cs_perturb(M_base, k, d, CS_STEP)
            col = 2 * k + d
            dfoil_s[:, col] = np.imag(M_p.foil.s) / CS_STEP
            dsstag[col] = float(np.imag(M_p.isol.sstag) / CS_STEP)
            dwake_s[:, col] = np.imag(M_p.wake.s) / CS_STEP
            dwake_x[:, :, col] = np.imag(M_p.wake.x) / CS_STEP
            dwgap[:, col] = np.imag(M_p.vsol.wgap) / CS_STEP
            dxi[:, col] = np.imag(M_p.isol.xi) / CS_STEP
        if (k + 1) % 10 == 0 or k == N - 1:
            print(f"  CS sweep: {k+1}/{N} nodes done "
                  f"(elapsed {time.time()-t0:.1f}s)")

    with h5py.File(OUTPUT_FILE, 'w') as f:
        g = f.create_group('geom_chain')
        g.create_dataset('foil_x', data=foil_x)
        g.create_dataset('foil_s', data=foil_s)
        g.create_dataset('foil_gam', data=foil_gam)
        g.create_dataset('wake_x', data=wake_x)
        g.create_dataset('wake_s', data=wake_s)
        g.create_dataset('sstag', data=sstag)
        g.create_dataset('Istag', data=Istag)
        g.create_dataset('xi', data=xi)
        g.create_dataset('wgap', data=wgap)
        g.create_dataset('chord', data=chord)
        g.create_dataset('wakelen', data=wakelen)
        g.create_dataset('alpha', data=float(ALPHA))
        g.create_dataset('dfoil_s_dx', data=dfoil_s)
        g.create_dataset('dsstag_dx', data=dsstag)
        g.create_dataset('dwake_s_dx', data=dwake_s)
        g.create_dataset('dwake_x_dx', data=dwake_x.reshape(2 * Nw, twoN))
        g.create_dataset('dwgap_dx', data=dwgap)
        g.create_dataset('dxi_dx', data=dxi)

    print(f"Wrote {OUTPUT_FILE} ({time.time()-t0:.1f}s total)")


if __name__ == '__main__':
    main()
