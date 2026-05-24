#!/usr/bin/env python3
"""CS reference data for viscous_output_ping_test.cpp.

Replaces the in-test FD validations of:
  - cd_U (CdUMatchesFD): ∂cd/∂Uw at the wake-end station.
  - dR_V/∂alpha (DRVDAlphaMatchesFD): partial at fixed glob.U.

Both computed via mfoil_cs complex-step perturbation. Test config matches
the C++ fixture exactly: NACA 0012, alpha=3, Re=1e6, Ma=0.0, npanel=100.
"""
import os
import sys
import numpy as np
import h5py
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'mfoil'))
import mfoil_cs  # noqa: E402
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil  # noqa: E402

CS_STEP = 1e-30
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'viscous_output_ping_reference.h5')

NACA = '0012'
NPANEL = 100
ALPHA = 3.0
RE = 1e6
MA = 0.0


def solve_base():
    M = mfoil.mfoil(naca=NACA, npanel=NPANEL)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=ALPHA, Re=RE, Ma=MA)
    M.solve()
    return M


def compute_RV_complex(M):
    """Full residual R_V = [R_BL; R_ue], complex-safe."""
    N = M.foil.N
    Nsys = M.glob.Nsys
    M.glob.R = np.zeros(3 * Nsys, dtype=complex)
    M.glob.R_U = sparse.lil_matrix((3 * Nsys, 4 * Nsys), dtype=complex)
    M.glob.R_x = sparse.lil_matrix((3 * Nsys, Nsys), dtype=complex)
    M.glob.realloc = False
    mfoil.build_glob_sys(M)
    R = np.zeros(4 * Nsys, dtype=complex)
    R[:3 * Nsys] = M.glob.R
    ue = M.glob.U[3, :]
    ds = M.glob.U[1, :]
    phi = M.vsol.sigma_m @ (ds * ue)  # [Npanel]
    sgnue_ext = np.concatenate([M.isol.sgnue, np.ones(Nsys - N)])
    src = sgnue_ext * (M.vsol.ue_sigma @ phi)
    ueinv = mfoil.get_ueinv(M)
    R[3 * Nsys:] = ue - ueinv - src
    return R


def main():
    M = solve_base()
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw
    print(f'Base: cl={M.post.cl:.10f} cd={M.post.cd:.10f}')
    print(f'  N={N} Nw={Nw} Nsys={Nsys}')

    iw = M.vsol.Is[2][-1]

    # ----- ∂cd/∂Uw at end-of-wake station ----------------------------------
    # Promote arrays to complex
    M.glob.U = M.glob.U.astype(complex)
    M.isol.gam = M.isol.gam.astype(complex)
    M.isol.uewi = M.isol.uewi.astype(complex)
    M.isol.sgnue = M.isol.sgnue.astype(complex)

    cd_U_cs = np.zeros(4)
    for k in range(4):
        M.glob.U[k, iw] += 1j * CS_STEP
        with complex_step_mode():
            mfoil.calc_force(M)
        cd_U_cs[k] = float(np.imag(M.post.cd)) / CS_STEP
        M.glob.U[k, iw] -= 1j * CS_STEP
    # Restore real cd
    mfoil.calc_force(M)
    print(f'\ncd_U via CS at iw={iw}:')
    for k in range(4):
        print(f'  cd_U[{k}] = {cd_U_cs[k]:.10e}')

    # ----- ∂R_V/∂alpha at fixed glob.U -------------------------------------
    # Need ue m matrix to recompute residual at perturbed alpha. Hold U fixed.
    M2 = solve_base()
    M2.glob.U = M2.glob.U.astype(complex)
    M2.isol.gam = M2.isol.gam.astype(complex)
    M2.isol.uewi = M2.isol.uewi.astype(complex)
    M2.isol.uewiref = M2.isol.uewiref.astype(complex)
    M2.isol.gamref = M2.isol.gamref.astype(complex)
    M2.isol.sgnue = M2.isol.sgnue.astype(complex)
    M2.vsol.ue_m = M2.vsol.ue_m.astype(complex)
    M2.vsol.ue_sigma = M2.vsol.ue_sigma.astype(complex)
    M2.vsol.sigma_m = M2.vsol.sigma_m.astype(complex)
    M2.wake.x = M2.wake.x.astype(complex)
    M2.wake.t = M2.wake.t.astype(complex) if isinstance(M2.wake.t, np.ndarray) else M2.wake.t
    M2.foil.s = M2.foil.s.astype(complex)
    alpha0 = M2.oper.alpha
    # Sanity check: R at base should be ~0 (converged viscous state)
    R_base = compute_RV_complex(M2)
    print(f'Sanity: |R_base| = {np.linalg.norm(np.real(R_base)):.6e} '
          f'(should be ~0)')
    M2.oper.alpha = alpha0 + 1j * CS_STEP
    with complex_step_mode():
        # Re-run the alpha-dependent inviscid chain at fixed glob.U AND
        # fixed sgnue (matching cppfoil's RebuildIsol which updates gam from
        # gamref @ [cos α, sin α] and rebuilds wake/ue_m when redowake):
        M2.oper.redowake = True
        mfoil.rebuild_isol(M2)
        mfoil.stagpoint_move(M2)
        mfoil.set_wake_gap(M2)
        R_pert = compute_RV_complex(M2)
    dRV_dalpha = np.imag(R_pert) / CS_STEP
    print(f'\ndR_V/dalpha L2 norm = {np.linalg.norm(dRV_dalpha):.6e}')
    print(f'  top 5 entries:')
    idx = np.argsort(-np.abs(dRV_dalpha))[:5]
    for i in idx:
        print(f'    row {i:4d} value={dRV_dalpha[i]:.6e}')

    with h5py.File(OUTPUT_FILE, 'w') as f:
        f.create_dataset('N', data=N)
        f.create_dataset('Nw', data=Nw)
        f.create_dataset('Nsys', data=Nsys)
        f.create_dataset('iw', data=iw)
        f.create_dataset('alpha', data=ALPHA)
        f.create_dataset('Re', data=RE)
        f.create_dataset('Ma', data=MA)
        f.create_dataset('cd_U', data=cd_U_cs)
        f.create_dataset('dRV_dalpha', data=dRV_dalpha)
    print(f'\nWrote {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
