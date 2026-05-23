#!/usr/bin/env python3
"""mfoil-CS complex-step ∂R_BL/∂Ma at fixed U for the ViscousMach test config
(NACA 0012, alpha=3, Re=1e6, Ma=0.3). Writes one HDF5 with dR_dMa[3*Nsys].

Compared against cppfoil's analytic glob.dR_dMa_times_Ma / Ma to localize the
remaining ~1e-5 discrepancy in ViscousMach dcl/dcm.
"""
import os
import sys
import numpy as np
import h5py
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'temp'))
import mfoil_cs  # noqa: E402
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil  # noqa: E402

CS_STEP = 1e-30
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'drv_dma_cs_reference.h5')

NACA = '0012'
NPANEL = 100
ALPHA = 3.0
RE = 1e6
MA = 0.3


def compute_R_BL_complex(M):
    """Build R_BL [3*Nsys] with complex dtype (CS-safe). Mirrors the helper
    in extract_comprehensive_reference.py."""
    Nsys = M.glob.Nsys
    M.glob.R = np.zeros(3 * Nsys, dtype=complex)
    M.glob.R_U = sparse.lil_matrix((3 * Nsys, 4 * Nsys), dtype=complex)
    M.glob.R_x = sparse.lil_matrix((3 * Nsys, Nsys), dtype=complex)
    M.glob.realloc = False
    mfoil.build_glob_sys(M)
    return M.glob.R.copy()


def main():
    M = mfoil.mfoil(naca=NACA, npanel=NPANEL)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=ALPHA, Re=RE, Ma=MA)
    M.solve()
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw
    print(f'Base: cl={M.post.cl:.10f} cd={M.post.cd:.10f} cm={M.post.cm:.10f}')
    print(f'N={N} Nw={Nw} Nsys={Nsys}')

    # Promote BL state to complex (held fixed across the Ma perturbation).
    M.glob.U = M.glob.U.astype(complex)
    M.isol.gam = M.isol.gam.astype(complex)
    M.isol.uewi = M.isol.uewi.astype(complex)
    M.isol.uewiref = M.isol.uewiref.astype(complex)
    M.isol.sgnue = M.isol.sgnue.astype(complex)

    Ma0 = M.oper.Ma
    with complex_step_mode():
        M.oper.Ma = Ma0 + 1j * CS_STEP
        mfoil.init_thermo(M)
        R_pert = compute_R_BL_complex(M)

    dR_dMa = np.imag(R_pert) / CS_STEP

    print(f'\ndR_dMa stats: |.|={np.linalg.norm(dR_dMa):.6e} '
          f'max|.|={np.max(np.abs(dR_dMa)):.6e}')
    # Top 10 entries
    idx = np.argsort(-np.abs(dR_dMa))[:10]
    for i in idx:
        print(f'  row {i:4d} station {i//3:3d} eq {i%3}: '
              f'dR/dMa = {dR_dMa[i]:.6e}')

    with h5py.File(OUTPUT_FILE, 'w') as f:
        f.create_dataset('N', data=N)
        f.create_dataset('Nw', data=Nw)
        f.create_dataset('Nsys', data=Nsys)
        f.create_dataset('Ma', data=float(Ma0))
        f.create_dataset('U', data=np.real(M.glob.U))
        f.create_dataset('dR_dMa', data=dR_dMa)
    print(f'Wrote {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
