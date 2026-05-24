#!/usr/bin/env python3
"""
Extract mfoil forward-pass reference data: R, R_x, R_U from build_glob_sys
at the converged state. Used to find forward-code differences vs C++ BuildGlobSys.

Usage:
  cd mfoil
  python ../scirocco/tests/extract_forward_reference.py
"""

import sys
import os
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'mfoil'))
import mfoil

NACA = '0012'
NPANEL = 100
ALPHA = 3.0
RE = 1e6
MA = 0.0
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'forward_reference.h5')


def main():
    print("Solving base case...")
    M = mfoil.mfoil(naca=NACA, npanel=NPANEL)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=ALPHA, Re=RE, Ma=MA)
    M.solve()

    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw
    print(f"  N={N}, Nw={Nw}, Nsys={Nsys}")
    print(f"  cl={M.post.cl:.10f}, cd={M.post.cd:.10f}")

    # Rebuild BL system at converged state (force realloc to get fresh matrices)
    M.glob.realloc = True
    mfoil.build_glob_sys(M)

    # Extract assembled matrices
    R = M.glob.R.copy()
    R_x = M.glob.R_x.toarray()   # [3*Nsys, Nsys] dense
    R_U = M.glob.R_U.toarray()   # [3*Nsys, 4*Nsys] dense

    print(f"  |R| = {np.linalg.norm(R):.6e}")
    print(f"  R_x shape = {R_x.shape}, nnz = {np.count_nonzero(np.abs(R_x) > 1e-12)}")
    print(f"  R_U shape = {R_U.shape}, nnz = {np.count_nonzero(np.abs(R_U) > 1e-12)}")

    # Per-station R and R_x for fine-grained debugging
    station_R = []    # list of 3-vectors
    station_Rx = []   # list of 3x2 matrices
    station_info = [] # (surface, local_idx, global_idx, type_str)
    for si in range(3):
        Is = M.vsol.Is[si]
        xi_s = M.isol.xi[Is]
        Ns = len(Is)
        U = M.glob.U[:, Is]
        Aux = np.zeros(Ns)
        param = mfoil.build_param(M, si)
        if si == 2:
            Aux[:] = M.vsol.wgap

        i0 = 1 if (si < 2) and (xi_s[0] < 1e-8 * xi_s[-1]) else 0

        # First point
        if si < 2:
            Ip = [i0, i0 + 1]
            Ust, Ust_U, Ust_x, xst = mfoil.stagnation_state(U[:, Ip], xi_s[Ip])
            param.turb, param.simi = False, True
            R1, R1_Ut, R1_x_raw = mfoil.residual_station(
                param, np.r_[xst, xst],
                np.stack((Ust, Ust), axis=-1), Aux[[i0, i0]])
            param.simi = False
            R1_Ust = R1_Ut[:, 0:4] + R1_Ut[:, 4:8]
            R1_x = np.dot(R1_Ust, Ust_x)  # chained through stagnation state
            station_R.append(R1.copy())
            station_Rx.append(R1_x.copy())
            station_info.append((si, i0, Is[i0], 'stag'))
        else:
            R1, R1_U, J = mfoil.wake_sys(M, param)
            station_R.append(R1.copy())
            station_Rx.append(np.zeros((3, 2)))  # no xi dep for wake init
            station_info.append((si, 0, Is[0], 'wake_init'))
            param.turb, param.wake = True, True

        # March
        for i in range(i0 + 1, Ns):
            Ip = [i - 1, i]
            tran = M.vsol.turb[Is[i - 1]] ^ M.vsol.turb[Is[i]]
            if tran:
                Ri, Ri_U, Ri_x = mfoil.residual_transition(
                    M, param, xi_s[Ip], U[:, Ip], Aux[Ip])
                param.turb = True
                stype = 'tran'
            else:
                Ri, Ri_U, Ri_x = mfoil.residual_station(
                    param, xi_s[Ip], U[:, Ip], Aux[Ip])
                stype = 'turb' if param.turb else 'lam'

            station_R.append(Ri.copy())
            station_Rx.append(Ri_x.copy())
            station_info.append((si, i, Is[i], stype))

    n_stations = len(station_R)
    print(f"  {n_stations} stations extracted")

    # Pack per-station data into arrays
    all_R = np.array(station_R)          # [n_stations, 3]
    all_Rx = np.array(station_Rx)        # [n_stations, 3, 2]
    all_si = np.array([s[0] for s in station_info], dtype=np.int64)
    all_local = np.array([s[1] for s in station_info], dtype=np.int64)
    all_global = np.array([s[2] for s in station_info], dtype=np.int64)

    # Write HDF5
    print(f"Writing {OUTPUT_FILE}...")
    with h5py.File(OUTPUT_FILE, 'w') as f:
        g = f.create_group('forward')

        # Assembled BL system
        g.create_dataset('R', data=R)
        g.create_dataset('R_x', data=R_x)
        g.create_dataset('R_U', data=R_U)

        # State
        g.create_dataset('U', data=M.glob.U)
        g.create_dataset('xi', data=M.isol.xi)
        g.create_dataset('turb', data=np.array(M.vsol.turb, dtype=np.int64))
        g.create_dataset('wgap', data=M.vsol.wgap)

        # Geometry
        g.create_dataset('foil_x', data=M.foil.x[:, :N])
        g.create_dataset('foil_s', data=M.foil.s)

        # Surface indices
        g.create_dataset('Is0', data=np.array(M.vsol.Is[0], dtype=np.int64))
        g.create_dataset('Is1', data=np.array(M.vsol.Is[1], dtype=np.int64))
        g.create_dataset('Is2', data=np.array(M.vsol.Is[2], dtype=np.int64))

        # Stagnation info
        g.create_dataset('Istag', data=np.array(M.isol.Istag, dtype=np.int64))
        g.create_dataset('sstag_ue', data=np.array(M.isol.sstag_ue))
        g.create_dataset('sgnue', data=M.isol.sgnue)

        # Scalars
        g.create_dataset('cl_base', data=float(M.post.cl))
        g.create_dataset('cd_base', data=float(M.post.cd))

        # Per-station data
        sg = f.create_group('stations')
        sg.create_dataset('R', data=all_R)
        sg.create_dataset('R_x', data=all_Rx)
        sg.create_dataset('surface', data=all_si)
        sg.create_dataset('local_idx', data=all_local)
        sg.create_dataset('global_idx', data=all_global)

        # Primal (pre-VJP) quantities used by ForwardReferenceTest.PrimalComparison
        pg = f.create_group('primal')
        pg.create_dataset('wake_x', data=M.wake.x)
        pg.create_dataset('wake_t', data=M.wake.t)
        pg.create_dataset('ue_sigma', data=np.array(M.vsol.ue_sigma))
        # sigma_m may be sparse — densify
        sigma_m_dense = M.vsol.sigma_m.toarray() if hasattr(M.vsol.sigma_m, 'toarray') else np.asarray(M.vsol.sigma_m)
        pg.create_dataset('sigma_m', data=sigma_m_dense)
        pg.create_dataset('gamref', data=M.isol.gamref)
        pg.create_dataset('uewiref', data=M.isol.uewiref)

    print("Done.")


if __name__ == '__main__':
    main()
