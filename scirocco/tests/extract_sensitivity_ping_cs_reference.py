#!/usr/bin/env python3
"""mfoil_cs complex-step references for sensitivity_ping_test.cpp.

Panel-influence derivative tests use a fixed Xj, xi geometry. We CS-perturb
each Xj or xi component through mfoil's matching panel function and store
the resulting Jacobian. Force partial test uses NACA 2412 at α=5.
"""
import os
import sys
import copy
import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'mfoil'))
import mfoil_cs  # noqa: E402
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil  # noqa: E402

CS_STEP = 1e-30
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'sensitivity_ping_reference.h5')


def cs_perturb(f, base, idx, step=CS_STEP):
    """CS-derivative of f(perturbed_base) at base, perturbing flat index idx."""
    arr = np.array(base, dtype=complex)
    arr.flat[idx] += 1j * step
    with complex_step_mode():
        out = f(arr)
    return out


def extract_panel(grp):
    """Panel-influence derivatives at fixed Xj, xi."""
    Xj = np.array([[0.3, 0.6], [0.05, -0.02]])
    xi = np.array([0.45, 0.15])
    grp.create_dataset('Xj', data=Xj)
    grp.create_dataset('xi', data=xi)

    # mfoil's panel routines have different return signatures.
    # panel_linvortex_stream(Xj, xi) → (a, b)
    # panel_linvortex_velocity(Xj, xi, vdir=None, onmid=False) → (a, b) where a,b are 2-vectors
    # panel_constsource_stream(Xj, xi) → scalar
    # panel_constsource_velocity(Xj, xi, vdir=None) → 2-vector
    # panel_linsource_stream(Xj, xi) → (a, b)
    # panel_linsource_velocity(Xj, xi, vdir=None) → (a, b) 2-vectors each

    # ---- LinVortexStream ----
    def lvs(Xj_, xi_):
        return mfoil.panel_linvortex_stream(Xj_, xi_)
    a_base, b_base = lvs(Xj, xi)
    a_Xj = np.zeros(4)
    b_Xj = np.zeros(4)
    for k in range(4):
        col = k // 2  # k=0,1 → col 0 (Xj[:,0]); k=2,3 → col 1 (Xj[:,1])
        row = k % 2   # k=0,2 → row 0 (x); k=1,3 → row 1 (y)
        # Actually the C++ test maps: if (k<2) Xj(k,0) else Xj(k-2,1).
        # So k=0: Xj(0,0)=x of pt1; k=1: Xj(1,0)=y of pt1;
        #    k=2: Xj(0,1)=x of pt2; k=3: Xj(1,1)=y of pt2.
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lvs(Xjc, xi)
        a_Xj[k] = np.imag(a) / CS_STEP
        b_Xj[k] = np.imag(b) / CS_STEP
    a_xi = np.zeros(2)
    b_xi = np.zeros(2)
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lvs(Xj, xic)
        a_xi[k] = np.imag(a) / CS_STEP
        b_xi[k] = np.imag(b) / CS_STEP
    sub = grp.create_group('LinVortexStream')
    sub.create_dataset('a_Xj', data=a_Xj)
    sub.create_dataset('b_Xj', data=b_Xj)
    sub.create_dataset('a_xi', data=a_xi)
    sub.create_dataset('b_xi', data=b_xi)

    # ---- ConstSourceStream (scalar) ----
    def css(Xj_, xi_):
        return mfoil.panel_constsource_stream(Xj_, xi_)
    val_Xj = np.zeros(4)
    for k in range(4):
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            v = css(Xjc, xi)
        val_Xj[k] = np.imag(v) / CS_STEP
    val_xi = np.zeros(2)
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            v = css(Xj, xic)
        val_xi[k] = np.imag(v) / CS_STEP
    sub = grp.create_group('ConstSourceStream')
    sub.create_dataset('val_Xj', data=val_Xj)
    sub.create_dataset('val_xi', data=val_xi)

    # ---- LinVortexVelocity (2-vector a, b) ----
    def lvv(Xj_, xi_):
        return mfoil.panel_linvortex_velocity(Xj_, xi_, None, False)
    a_Xj = np.zeros((2, 4))
    b_Xj = np.zeros((2, 4))
    for k in range(4):
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lvv(Xjc, xi)
        a_Xj[:, k] = np.imag(a) / CS_STEP
        b_Xj[:, k] = np.imag(b) / CS_STEP
    a_xi = np.zeros((2, 2))
    b_xi = np.zeros((2, 2))
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lvv(Xj, xic)
        a_xi[:, k] = np.imag(a) / CS_STEP
        b_xi[:, k] = np.imag(b) / CS_STEP
    sub = grp.create_group('LinVortexVelocity')
    sub.create_dataset('a_Xj', data=a_Xj)
    sub.create_dataset('b_Xj', data=b_Xj)
    sub.create_dataset('a_xi', data=a_xi)
    sub.create_dataset('b_xi', data=b_xi)

    # ---- ConstSourceVelocity (2-vector) ----
    def csv(Xj_, xi_):
        return mfoil.panel_constsource_velocity(Xj_, xi_, None)
    val_Xj = np.zeros((2, 4))
    for k in range(4):
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            v = csv(Xjc, xi)
        val_Xj[:, k] = np.imag(v) / CS_STEP
    val_xi = np.zeros((2, 2))
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            v = csv(Xj, xic)
        val_xi[:, k] = np.imag(v) / CS_STEP
    sub = grp.create_group('ConstSourceVelocity')
    sub.create_dataset('val_Xj', data=val_Xj)
    sub.create_dataset('val_xi', data=val_xi)

    # ---- LinSourceStream ----
    def lss(Xj_, xi_):
        return mfoil.panel_linsource_stream(Xj_, xi_)
    a_Xj = np.zeros(4)
    b_Xj = np.zeros(4)
    for k in range(4):
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lss(Xjc, xi)
        a_Xj[k] = np.imag(a) / CS_STEP
        b_Xj[k] = np.imag(b) / CS_STEP
    a_xi = np.zeros(2)
    b_xi = np.zeros(2)
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lss(Xj, xic)
        a_xi[k] = np.imag(a) / CS_STEP
        b_xi[k] = np.imag(b) / CS_STEP
    sub = grp.create_group('LinSourceStream')
    sub.create_dataset('a_Xj', data=a_Xj)
    sub.create_dataset('b_Xj', data=b_Xj)
    sub.create_dataset('a_xi', data=a_xi)
    sub.create_dataset('b_xi', data=b_xi)

    # ---- LinSourceVelocity ----
    def lsv(Xj_, xi_):
        return mfoil.panel_linsource_velocity(Xj_, xi_, None)
    a_Xj = np.zeros((2, 4))
    b_Xj = np.zeros((2, 4))
    for k in range(4):
        Xjc = np.array(Xj, dtype=complex)
        if k < 2:
            Xjc[k, 0] += 1j * CS_STEP
        else:
            Xjc[k - 2, 1] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lsv(Xjc, xi)
        a_Xj[:, k] = np.imag(a) / CS_STEP
        b_Xj[:, k] = np.imag(b) / CS_STEP
    a_xi = np.zeros((2, 2))
    b_xi = np.zeros((2, 2))
    for k in range(2):
        xic = np.array(xi, dtype=complex)
        xic[k] += 1j * CS_STEP
        with complex_step_mode():
            a, b = lsv(Xj, xic)
        a_xi[:, k] = np.imag(a) / CS_STEP
        b_xi[:, k] = np.imag(b) / CS_STEP
    sub = grp.create_group('LinSourceVelocity')
    sub.create_dataset('a_Xj', data=a_Xj)
    sub.create_dataset('b_Xj', data=b_Xj)
    sub.create_dataset('a_xi', data=a_xi)
    sub.create_dataset('b_xi', data=b_xi)


def extract_force_partial(grp):
    """dcl/dx_partial at fixed gamma, NACA 2412 α=5. CS-perturb each foil
    coordinate component, compute Im(cl)/eps at fixed isol.gam.

    Matches C++ test: foil.Solve() then perturbs foil.x at fixed gamma.
    """
    M = mfoil.mfoil(naca='2412', npanel=80)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=5.0)
    M.solve()
    N = M.foil.N

    # We want partial dcl/dx at FIXED gamma. So we promote foil.x to complex,
    # perturb, but don't re-solve gamma. We DO need to recompute foil.s and
    # then calc_force.
    dcl_dx = np.zeros((2, N))
    base_x = np.array(M.foil.x, dtype=float)

    def recompute_geom_complex(M_):
        x = M_.foil.x
        Np = M_.foil.N
        ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
        s = np.zeros(Np, dtype=x.dtype)
        s[1:] = np.cumsum(ds)
        M_.foil.s = s

    for i in range(N):
        for d in range(2):
            Mp = copy.deepcopy(M)
            Mp.foil.x = Mp.foil.x.astype(complex)
            Mp.isol.gam = Mp.isol.gam.astype(complex)
            Mp.isol.sgnue = Mp.isol.sgnue.astype(complex)
            Mp.foil.x[d, i] += 1j * CS_STEP
            recompute_geom_complex(Mp)
            with complex_step_mode():
                # Fixed gamma. calc_force uses gam, foil.x, alpha → cl.
                mfoil.calc_force(Mp)
            dcl_dx[d, i] = float(np.imag(Mp.post.cl)) / CS_STEP
    grp.create_dataset('foil_x', data=base_x)
    grp.create_dataset('dcl_dx_partial', data=dcl_dx)
    grp.create_dataset('cl_base', data=float(M.post.cl))


def main():
    with h5py.File(OUTPUT_FILE, 'w') as f:
        print('Extracting panel influences ...')
        extract_panel(f.create_group('panels'))
        print('Extracting force partials ...')
        extract_force_partial(f.create_group('force_partial'))
    print(f'Wrote {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
