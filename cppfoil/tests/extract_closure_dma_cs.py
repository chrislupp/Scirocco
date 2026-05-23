#!/usr/bin/env python3
"""Extract mfoil_cs complex-step ∂closure/∂Ma at every BL station for the
ViscousMach test config (NACA 0012, alpha=3, Re=1e6, Ma=0.3). Writes one
HDF5 with arrays of CS derivatives per closure, indexed by global node.

Used by C++ probe to validate analytic _dMa chain at machine precision.
"""
import os
import sys
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'temp'))
import mfoil_cs  # noqa: E402 — patches mfoil for CS compat
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil  # noqa: E402

CS_STEP = 1e-30
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'closure_dma_cs_reference.h5')


def solve_base():
    """Solve viscous at Ma=0.3 and return the converged mfoil object."""
    M = mfoil.mfoil(naca='0012', npanel=100)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=3.0, Re=1e6, Ma=0.3)
    M.solve()
    return M


def extract():
    M = solve_base()
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw

    # Build per-station turb flags and surface assignment
    turb_flags = np.array(M.vsol.turb, dtype=bool)
    # Wake nodes are always turb
    wake_nodes = list(M.vsol.Is[2])

    # We capture U values from converged state (real-valued).
    U_real = M.glob.U[:, :Nsys].copy()
    foil_x_real = M.foil.x[:, :N].copy()

    closures = [
        'GetUk',     # GetUk(U[3], param)
        'GetMach2',  # GetMach2(U, param)
        'GetHk',     # GetHk(U, param)
        'GetHs_T',   # GetHs(U, param.turb=True)
        'GetHs_F',   # GetHs(U, param.turb=False)
        'GetRet',
        'GetUs_T',   # GetUs with turb=T at the param's surface
        'GetUs_F',   # GetUs with turb=F
        'GetCf_T',
        'GetCf_F',
        'GetCDi_T',  # full GetCDi with turb=T (wake => uses LamWake floor)
        'GetCDi_F',
        'GetCfxt_T',
        'GetCfxt_F',
        'GetCDixt_T',
        'GetCDixt_F',
        'GetDe_T',
        'GetDe_F',
        'GetCteq_T',
        'GetDamp',   # always lam in mfoil
    ]
    results = {name: np.zeros(Nsys, dtype=float) for name in closures}

    def perturb(M, Ma_complex):
        """Set Ma and re-init thermo. Returns a fresh copy at the perturbed
        Ma. Other state (U, foil.x, isol.gam, ...) stays REAL — we are
        computing ∂closure/∂Ma at fixed U."""
        # Build a perturbed param. mfoil's init_thermo modifies M.param
        # KTb, KTl, H0, cps, mu0, rho0, muinf, Minf, gam — using M.oper.Ma.
        # Strategy: set M.oper.Ma to the complex value, call init_thermo,
        # then call the closure. Restore afterward.
        save_oper_Ma = M.oper.Ma
        save_param = dict(M.param.__dict__)
        with complex_step_mode():
            M.oper.Ma = Ma_complex
            mfoil.init_thermo(M)
        return save_oper_Ma, save_param

    def restore(M, save_oper_Ma, save_param):
        M.oper.Ma = save_oper_Ma
        for k, v in save_param.items():
            setattr(M.param, k, v)

    Ma0 = M.oper.Ma

    for node in range(Nsys):
        U = U_real[:, node]
        is_wake = (node >= N)
        is_turb = (turb_flags[node] if node < N else True)

        # Determine xi for this node (used by GetCfxt/GetCDixt)
        xi = M.isol.xi[node]

        # Build per-node param view (turb + wake flags)
        def make_param(turb_flag, wake_flag):
            # mfoil.param doesn't have direct turb/wake; instead we pass
            # to closures. Use the standard build_param + station_param
            # logic but with turb override.
            # Approach: just set param.turb and param.wake directly.
            M.param.turb = bool(turb_flag)
            M.param.wake = bool(wake_flag)
            M.param.simi = False

        def cs(closure_fn, turb_flag, wake_flag):
            save = perturb(M, Ma0 + 1j * CS_STEP)
            make_param(turb_flag, wake_flag)
            try:
                with complex_step_mode():
                    v = closure_fn(U, M.param)
                    val = v[0] if isinstance(v, tuple) else v
                d = float(np.imag(val)) / CS_STEP
            finally:
                restore(M, *save)
                M.param.wake = False
                M.param.turb = False
            return d

        # GetUk takes (ue, param) — different signature
        def cs_uk():
            save = perturb(M, Ma0 + 1j * CS_STEP)
            try:
                with complex_step_mode():
                    v = mfoil.get_uk(U[3], M.param)
                    val = v[0]
                d = float(np.imag(val)) / CS_STEP
            finally:
                restore(M, *save)
            return d

        def cs_with_xi(closure_fn, turb_flag, wake_flag):
            # GetCfxt / GetCDixt take (U, x, param)
            save = perturb(M, Ma0 + 1j * CS_STEP)
            make_param(turb_flag, wake_flag)
            try:
                with complex_step_mode():
                    v = closure_fn(U, xi, M.param)
                    val = v[0]
                d = float(np.imag(val)) / CS_STEP
            finally:
                restore(M, *save)
                M.param.wake = False
                M.param.turb = False
            return d

        results['GetUk'][node]      = cs_uk()
        results['GetMach2'][node]   = cs(mfoil.get_Mach2,  is_turb, is_wake)
        results['GetHk'][node]      = cs(mfoil.get_Hk,     is_turb, is_wake)
        results['GetHs_T'][node]    = cs(mfoil.get_Hs,     True,    is_wake)
        results['GetHs_F'][node]    = cs(mfoil.get_Hs,     False,   is_wake)
        results['GetRet'][node]     = cs(mfoil.get_Ret,    is_turb, is_wake)
        results['GetUs_T'][node]    = cs(mfoil.get_Us,     True,    is_wake)
        results['GetUs_F'][node]    = cs(mfoil.get_Us,     False,   is_wake)
        results['GetCf_T'][node]    = cs(mfoil.get_cf,     True,    is_wake)
        results['GetCf_F'][node]    = cs(mfoil.get_cf,     False,   is_wake)
        results['GetCDi_T'][node]   = cs(mfoil.get_cDi,    True,    is_wake)
        results['GetCDi_F'][node]   = cs(mfoil.get_cDi,    False,   is_wake)
        results['GetCfxt_T'][node]  = cs_with_xi(mfoil.get_cfxt,  True,  is_wake)
        results['GetCfxt_F'][node]  = cs_with_xi(mfoil.get_cfxt,  False, is_wake)
        results['GetCDixt_T'][node] = cs_with_xi(mfoil.get_cDixt, True,  is_wake)
        results['GetCDixt_F'][node] = cs_with_xi(mfoil.get_cDixt, False, is_wake)
        results['GetDe_T'][node]    = cs(mfoil.get_de,     True,    is_wake)
        results['GetDe_F'][node]    = cs(mfoil.get_de,     False,   is_wake)
        results['GetCteq_T'][node]  = cs(mfoil.get_cteq,   True,    is_wake)
        results['GetDamp'][node]    = cs(mfoil.get_damp,   False,   False)

        if node % 10 == 0 or node == Nsys - 1:
            print(f'  node {node}/{Nsys-1}  '
                  f'GetCDi_T={results["GetCDi_T"][node]:.4e}  '
                  f'GetUs_T={results["GetUs_T"][node]:.4e}  '
                  f'GetCfxt_T={results["GetCfxt_T"][node]:.4e}')

    with h5py.File(OUTPUT_FILE, 'w') as f:
        f.create_dataset('N', data=N)
        f.create_dataset('Nw', data=Nw)
        f.create_dataset('Nsys', data=Nsys)
        f.create_dataset('Ma', data=float(Ma0))
        f.create_dataset('U', data=U_real)
        f.create_dataset('xi', data=M.isol.xi[:Nsys])
        f.create_dataset('turb', data=turb_flags.astype(np.int8))
        for name, arr in results.items():
            f.create_dataset(name, data=arr)
    print(f'Wrote {OUTPUT_FILE}')


if __name__ == '__main__':
    extract()
