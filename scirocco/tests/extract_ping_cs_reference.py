#!/usr/bin/env python3
"""mfoil_cs complex-step Jacobian references for ping_test.cpp.

Replaces the in-test FD convergence-rate checks with machine-precision CS
references. One HDF5 file with subgroups per test, each containing the
CS-derived Jacobian/derivative at the exact same inputs the C++ test uses.

Tests covered:
  - BLPropertyDerivatives:   13 closures × 3 regimes × 2 Re values
  - ResidualStation:         3 regimes × 2 Re values
  - ResidualTransition:      Re=2e6, Ma=0, lam→turb
  - ResidualStation_Stagnation: simi=True at alpha=1, Re=1e6
  - BuildGlobSys (explicit): full Jacobian at converged state (NACA 0012,
                              130 panels, alpha=3, Re=1e6, Ma=0)
  - WakeSys:                 R_U at converged state (same as above)
  - StagnationState:         specific Ust/Ust_U/Ust_x at hardcoded U,x
  - CalcForce:               cl_ue, cl_alpha at converged state
"""
import os
import sys
import copy
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
                           'ping_reference.h5')


def make_param(Re, Ma, turb, wake):
    """Replicate scirocco's MakeParam(Re, Ma, turb, wake) using mfoil."""
    M = mfoil.mfoil(naca='0012', npanel=130)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=3.0, Re=Re, Ma=Ma)
    mfoil.init_thermo(M)
    M.param.turb = turb
    M.param.wake = wake
    M.param.simi = False
    return M


def cs_jacobian_4(func, U_base, M):
    """Generic CS Jacobian of a closure func(U, param) → scalar, w.r.t. U[0..3].
    Returns Vector4 of dv/dU[k]."""
    U_c = np.array(U_base, dtype=complex)
    J = np.zeros(4)
    for k in range(4):
        U_c[k] = U_base[k] + 1j * CS_STEP
        with complex_step_mode():
            v = func(U_c, M.param)
            val = v[0] if isinstance(v, tuple) else v
        J[k] = float(np.imag(val)) / CS_STEP
        U_c[k] = U_base[k]
    return J


def extract_bl_properties(grp):
    """A. BLPropertyDerivatives."""
    Ulam = np.array([0.01, 0.02, 8.4, 0.9])
    Uturb = np.array([0.023, 0.05, 0.031, 1.1])
    regimes = [
        ('lam', False, False, Ulam),
        ('turb', True, False, Uturb),
        ('wake', True, True, Uturb),
    ]
    Rev = [2e3, 1e5]

    funcs = {
        'GetHk': mfoil.get_Hk,
        'GetRet': mfoil.get_Ret,
        'GetCf': mfoil.get_cf,
        'GetCDi': mfoil.get_cDi,
        'GetHs': mfoil.get_Hs,
        'GetHss': mfoil.get_Hss,
        'GetUs': mfoil.get_Us,
        'GetCteq': mfoil.get_cteq,
        'GetCttr': mfoil.get_cttr,
        'GetDe': mfoil.get_de,
        'GetDamp': mfoil.get_damp,
        'GetMach2': mfoil.get_Mach2,
        'GetRho': mfoil.get_rho,
    }

    for Re in Rev:
        for name, turb, wake, U in regimes:
            M = make_param(Re, 0.4, turb, wake)
            sub = grp.require_group(f'{name}_Re{int(Re)}')
            for fname, ffunc in funcs.items():
                J = cs_jacobian_4(ffunc, U, M)
                sub.create_dataset(fname, data=J)
            sub.create_dataset('U', data=U)


def extract_residual_station(grp):
    """A2. ResidualStationDerivatives. Uses Umat=[U1,U2] of shape 4x2, xi
    [2], Aux [2]. We need the Jacobian R_U[3,8] and R_x[3,2] which the C++
    test combines with a dx step. Easier to just store the full Jacobians."""
    Ulam = np.array([0.01, 0.02, 8.4, 0.9])
    Uturb = np.array([0.023, 0.05, 0.031, 1.1])
    regimes = [
        ('lam', False, False, Ulam),
        ('turb', True, False, Uturb),
        ('wake', True, True, Uturb),
    ]
    Rev = [2e3, 1e5]
    xi_base = np.array([0.7, 0.8])
    Aux = np.array([0.002, 0.0018])

    for Re in Rev:
        for name, turb, wake, U1base in regimes:
            U2base = np.array([U1base[0]*1.1, U1base[1]*0.8,
                               U1base[2]*0.9, U1base[3]*1.2])
            M = make_param(Re, 0.4, turb, wake)
            # Build Umat. Compute full Jacobian via CS.
            Umat = np.column_stack([U1base, U2base])

            # CS each U[k][col] → R_U column
            R_U = np.zeros((3, 8))
            for k in range(8):
                col = k // 4
                comp = k % 4
                Uc = np.array(Umat, dtype=complex)
                Uc[comp, col] = Umat[comp, col] + 1j * CS_STEP
                with complex_step_mode():
                    R, *_ = mfoil.residual_station(M.param, xi_base, Uc, Aux)
                R_U[:, k] = np.imag(R) / CS_STEP

            # CS each xi[col] → R_x column
            R_x = np.zeros((3, 2))
            for col in range(2):
                xc = np.array(xi_base, dtype=complex)
                xc[col] = xi_base[col] + 1j * CS_STEP
                with complex_step_mode():
                    R, *_ = mfoil.residual_station(M.param, xc, Umat, Aux)
                R_x[:, col] = np.imag(R) / CS_STEP

            sub = grp.require_group(f'{name}_Re{int(Re)}')
            sub.create_dataset('U1', data=U1base)
            sub.create_dataset('U2', data=U2base)
            sub.create_dataset('xi', data=xi_base)
            sub.create_dataset('Aux', data=Aux)
            sub.create_dataset('R_U', data=R_U)
            sub.create_dataset('R_x', data=R_x)


def extract_residual_transition(grp):
    """B. ResidualTransition at Re=2e6, Ma=0.

    The test fixture solves NACA 0012 at alpha=3, Re=1e6, Ma=0, then
    overrides state.oper.Re=2e6 and rebuilds InitThermo / BuildParam.
    mfoil equivalent: solve at the same params then override and re-init."""
    M = mfoil.mfoil(naca='0012', npanel=130)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=3.0, Re=1e6, Ma=0.0)
    M.solve()
    M.oper.Re = 2e6
    M.oper.Ma = 0.0
    mfoil.init_thermo(M)
    M.param.turb = False  # mfoil's residual_transition sets these internally

    U = np.zeros((4, 2))
    U[:, 0] = [0.01, 0.02, 8.95, 0.9]
    U[:, 1] = [0.013, 0.023, 0.028, 0.85]
    x = np.array([0.7, 0.8])
    Aux = np.zeros(2)

    R_U = np.zeros((3, 8))
    for k in range(8):
        col = k // 4
        comp = k % 4
        Uc = np.array(U, dtype=complex)
        Uc[comp, col] = U[comp, col] + 1j * CS_STEP
        with complex_step_mode():
            R, *_ = mfoil.residual_transition(M, M.param, x, Uc, Aux)
        R_U[:, k] = np.imag(R) / CS_STEP
    R_x = np.zeros((3, 2))
    for col in range(2):
        xc = np.array(x, dtype=complex)
        xc[col] = x[col] + 1j * CS_STEP
        with complex_step_mode():
            R, *_ = mfoil.residual_transition(M, M.param, xc, U, Aux)
        R_x[:, col] = np.imag(R) / CS_STEP

    grp.create_dataset('U', data=U)
    grp.create_dataset('x', data=x)
    grp.create_dataset('Aux', data=Aux)
    grp.create_dataset('R_U', data=R_U)
    grp.create_dataset('R_x', data=R_x)


def extract_residual_stagnation(grp):
    """C. ResidualStation_Stagnation: simi=True, alpha=1, Re=1e6, Ma=0."""
    M = mfoil.mfoil(naca='0012', npanel=130)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=1.0, Re=1e6, Ma=0.0)
    mfoil.init_thermo(M)
    M.param.turb = False
    M.param.wake = False
    M.param.simi = True

    U = np.array([0.00004673616, 0.000104289, 0.0, 0.11977917547])
    xval = 4.590816441485401e-05
    Aux = np.zeros(2)
    Umat = np.column_stack([U, U])
    xvec = np.array([xval, xval])

    R_U = np.zeros((3, 8))
    for k in range(8):
        col = k // 4
        comp = k % 4
        Uc = np.array(Umat, dtype=complex)
        Uc[comp, col] = Umat[comp, col] + 1j * CS_STEP
        with complex_step_mode():
            R, *_ = mfoil.residual_station(M.param, xvec, Uc, Aux)
        R_U[:, k] = np.imag(R) / CS_STEP
    R_x = np.zeros((3, 2))
    for col in range(2):
        xc = np.array(xvec, dtype=complex)
        xc[col] = xvec[col] + 1j * CS_STEP
        with complex_step_mode():
            R, *_ = mfoil.residual_station(M.param, xc, Umat, Aux)
        R_x[:, col] = np.imag(R) / CS_STEP

    grp.create_dataset('U', data=U)
    grp.create_dataset('x', data=xval)
    grp.create_dataset('R_U', data=R_U)
    grp.create_dataset('R_x', data=R_x)


def extract_stagnation_state(grp):
    """F. StagnationState at hardcoded U, x. Ust = stagnation_state(U, x)."""
    U = np.zeros((4, 2))
    U[:, 0] = [5e-5, 1.1e-4, 0.0, 0.0348]
    U[:, 1] = [4.9e-5, 1.09e-4, 0.0, 0.07397]
    x = np.array([5.18e-4, 1.1e-3])

    Ust_U = np.zeros((4, 8))
    for k in range(8):
        col = k // 4
        comp = k % 4
        Uc = np.array(U, dtype=complex)
        Uc[comp, col] = U[comp, col] + 1j * CS_STEP
        with complex_step_mode():
            Ust, _Ust_U, _Ust_x, _xst = mfoil.stagnation_state(Uc, x)
        Ust_U[:, k] = np.imag(Ust) / CS_STEP
    Ust_x = np.zeros((4, 2))
    for col in range(2):
        xc = np.array(x, dtype=complex)
        xc[col] = x[col] + 1j * CS_STEP
        with complex_step_mode():
            Ust, _Ust_U2, _Ust_x2, _xst2 = mfoil.stagnation_state(U, xc)
        Ust_x[:, col] = np.imag(Ust) / CS_STEP

    grp.create_dataset('U', data=U)
    grp.create_dataset('x', data=x)
    grp.create_dataset('Ust_U', data=Ust_U)
    grp.create_dataset('Ust_x', data=Ust_x)


def base_solve_fixture():
    """Match scirocco's PingTest fixture: NACA 0012, 130 panels, alpha=3,
    Re=1e6, Ma=0. Returns the converged mfoil object."""
    M = mfoil.mfoil(naca='0012', npanel=130)
    M.param.verb = 0
    M.param.doplot = False
    M.param.niglob = 25
    M.param.rtol = 1e-6
    M.setoper(alpha=3.0, Re=1e6, Ma=0.0)
    M.solve()
    return M


def extract_build_glob_sys(grp):
    """D. BuildGlobSys Jacobians: store full R_U and R_x at converged state.

    Also stores converged U and xi so the C++ test loads them and rebuilds
    the Jacobian at the EXACT same state (both scirocco and mfoil converge
    only to rtol=1e-6, so their independent U's differ at 1e-6 — enough to
    make the Jacobian comparison fail at machine precision)."""
    M = base_solve_fixture()
    N = M.foil.N
    Nsys = M.glob.Nsys
    grp.create_dataset('Nsys', data=Nsys)
    grp.create_dataset('R_base', data=np.real(M.glob.R))
    grp.create_dataset('U', data=np.real(M.glob.U))
    grp.create_dataset('xi', data=np.real(M.isol.xi))
    R_U_dense = M.glob.R_U.toarray() if sparse.issparse(M.glob.R_U) else M.glob.R_U
    R_x_dense = M.glob.R_x.toarray() if sparse.issparse(M.glob.R_x) else M.glob.R_x
    grp.create_dataset('R_U', data=R_U_dense.astype(float))
    grp.create_dataset('R_x', data=R_x_dense.astype(float))


def extract_wake_sys(grp):
    """E. WakeSys derivative R_U at the converged state."""
    M = base_solve_fixture()
    il = M.vsol.Is[0][-1]
    iu = M.vsol.Is[1][-1]
    iw = M.vsol.Is[2][0]

    # Promote things for CS
    M.glob.U = M.glob.U.astype(complex)
    M.isol.gam = M.isol.gam.astype(complex)

    # mfoil wake_sys mutates param; mfoil exposes wake_sys with signature
    # wake_sys(M, param) returning (R, R_U, J).
    # We CS-perturb glob.U[k, J[j]] for each j in [il, iu, iw] and each k in
    # 0..3, then read R imag part.
    R_U = np.zeros((3, 12))
    J = (il, iu, iw)
    for j_idx, station in enumerate(J):
        for k in range(4):
            M.glob.U[k, station] += 1j * CS_STEP
            with complex_step_mode():
                p = copy.copy(M.param)
                R, _R_U, J_ret = mfoil.wake_sys(M, p)
            R_U[:, 4 * j_idx + k] = np.imag(R) / CS_STEP
            M.glob.U[k, station] -= 1j * CS_STEP

    grp.create_dataset('R_U', data=R_U)
    grp.create_dataset('J', data=np.array([il, iu, iw], dtype=np.int64))


def extract_calc_force(grp):
    """G. CalcForce: cl_ue (N) and cl_alpha (scalar) at converged state."""
    M = base_solve_fixture()
    N = M.foil.N
    Nsys = M.glob.Nsys

    # Promote
    M.glob.U = M.glob.U.astype(complex)
    M.isol.gam = M.isol.gam.astype(complex)
    M.isol.uewi = M.isol.uewi.astype(complex)
    M.isol.sgnue = M.isol.sgnue.astype(complex)

    cl_ue = np.zeros(N)
    for i in range(N):
        M.glob.U[3, i] += 1j * CS_STEP
        with complex_step_mode():
            mfoil.calc_force(M)
        cl_ue[i] = float(np.imag(M.post.cl)) / CS_STEP
        M.glob.U[3, i] -= 1j * CS_STEP

    alpha0 = M.oper.alpha
    M.oper.alpha = alpha0 + 1j * CS_STEP
    with complex_step_mode():
        mfoil.rebuild_isol(M)
        mfoil.calc_force(M)
    cl_alpha = float(np.imag(M.post.cl)) / CS_STEP

    grp.create_dataset('cl_ue', data=cl_ue)
    grp.create_dataset('cl_alpha', data=cl_alpha)


def main():
    with h5py.File(OUTPUT_FILE, 'w') as f:
        print('Extracting BLPropertyDerivatives ...')
        extract_bl_properties(f.create_group('bl_props'))
        print('Extracting ResidualStation ...')
        extract_residual_station(f.create_group('residual_station'))
        print('Extracting ResidualTransition ...')
        extract_residual_transition(f.create_group('residual_transition'))
        print('Extracting ResidualStation_Stagnation ...')
        extract_residual_stagnation(f.create_group('residual_stagnation'))
        print('Extracting StagnationState ...')
        extract_stagnation_state(f.create_group('stagnation_state'))
        print('Extracting BuildGlobSys ...')
        extract_build_glob_sys(f.create_group('build_glob_sys'))
        print('Extracting WakeSys ...')
        extract_wake_sys(f.create_group('wake_sys'))
        print('Extracting CalcForce ...')
        extract_calc_force(f.create_group('calc_force'))
    print(f'Wrote {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
