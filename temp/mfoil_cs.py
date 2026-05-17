"""
Complex-step compatible wrapper for mfoil.

Monkey-patches ~12 functions in mfoil that break with complex arithmetic
(abs() on complex, np.arctan2, np.sign, complex-rejection checks, comparisons).

Usage:
    import mfoil_cs as mfoil
    M = mfoil.mfoil(naca='0012', npanel=100)
    # ... proceed with complex-step differentiation
"""

import sys
import os
import numpy as np
from scipy import sparse

# Ensure mfoil is importable from the same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mfoil as _mfoil_module
# Keep a reference to the module before `from mfoil import *` shadows the name
import mfoil

# ---------------------------------------------------------------------------
# Helper: complex-safe atan2
#   For real inputs, behaves like np.arctan2.
#   For complex inputs, uses first-order complex-step linearization:
#     atan2(y, x) ≈ atan2(Re(y), Re(x)) + j * Im(y*conj(x)) / (Re(x)^2 + Re(y)^2)
# ---------------------------------------------------------------------------
def _cs_atan2(y, x):
    yr, xr = np.real(y), np.real(x)
    base = np.arctan2(yr, xr)
    if np.isrealobj(y) and np.isrealobj(x):
        return base
    yi, xi = np.imag(y), np.imag(x)
    denom = xr**2 + yr**2
    if np.ndim(denom) == 0:
        if denom == 0:
            return base
    return base + 1j * (yi * xr - xi * yr) / denom


# ---------------------------------------------------------------------------
# Patched: TE_info  (line 877: abs() on cross product)
# ---------------------------------------------------------------------------
_orig_TE_info = _mfoil_module.TE_info

def _patched_TE_info(X):
    t1 = X[:, 0] - X[:, 1]; t1 = t1 / _mfoil_module.norm2(t1)
    t2 = X[:, -1] - X[:, -2]; t2 = t2 / _mfoil_module.norm2(t2)
    t = 0.5 * (t1 + t2); t = t / _mfoil_module.norm2(t)
    s = X[:, -1] - X[:, 0]
    hTE = -s[0] * t[1] + s[1] * t[0]
    dtdx = t1[0] * t2[1] - t2[0] * t1[1]
    p = s / _mfoil_module.norm2(s)
    # Use np.real() instead of abs() — sign is fixed by geometry convention
    tcp = np.real(t[0] * p[1] - t[1] * p[0])
    tdp = np.dot(t, p)
    return t, hTE, dtdx, tcp, tdp


# ---------------------------------------------------------------------------
# Patched: panel_info  (lines 914-915: np.arctan2 with complex args)
# ---------------------------------------------------------------------------
_orig_panel_info = _mfoil_module.panel_info

def _patched_panel_info(Xj, xi):
    xj1, zj1 = Xj[0, 0], Xj[1, 0]
    xj2, zj2 = Xj[0, 1], Xj[1, 1]

    t = np.array([xj2 - xj1, zj2 - zj1]); t = t / _mfoil_module.norm2(t)
    n = np.array([-t[1], t[0]])

    xz = np.array([(xi[0] - xj1), (xi[1] - zj1)])
    x = np.dot(xz, t)
    z = np.dot(xz, n)

    d = _mfoil_module.dist(xj2 - xj1, zj2 - zj1)
    r1 = _mfoil_module.dist(x, z)
    r2 = _mfoil_module.dist(x - d, z)
    theta1 = _cs_atan2(z, x)
    theta2 = _cs_atan2(z, x - d)

    return t, n, x, z, d, r1, r2, theta1, theta2


# ---------------------------------------------------------------------------
# Patched panel functions: use .real for r < ep comparisons and theta < 0
# ---------------------------------------------------------------------------
def _patched_panel_linvortex_stream(Xj, xi):
    t, n, x, z, d, r1, r2, theta1, theta2 = _patched_panel_info(Xj, xi)
    ep = 1e-9
    logr1 = 0. if (np.real(r1) < ep) else np.log(r1)
    logr2 = 0. if (np.real(r2) < ep) else np.log(r2)
    P1 = (0.5 / np.pi) * (z * (theta2 - theta1) - d + x * logr1 - (x - d) * logr2)
    P2 = x * P1 + (0.5 / np.pi) * (0.5 * r2**2 * logr2 - 0.5 * r1**2 * logr1 - r2**2 / 4 + r1**2 / 4)
    a = P1 - P2 / d
    b = P2 / d
    return a, b


def _patched_panel_constsource_velocity(Xj, xi, vdir=None):
    t, n, x, z, d, r1, r2, theta1, theta2 = _patched_panel_info(Xj, xi)
    ep = 1e-9
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (np.real(r1) < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0, 0, 0) if (np.real(r2) < ep) else (np.log(r2), theta1, theta2)
    u = (0.5 / np.pi) * (logr1 - logr2)
    w = (0.5 / np.pi) * (theta2 - theta1)
    a = np.array([u * t[0] + w * n[0], u * t[1] + w * n[1]])
    if vdir is not None:
        a = np.dot(a, vdir)
    return a


def _patched_panel_constsource_stream(Xj, xi):
    t, n, x, z, d, r1, r2, theta1, theta2 = _patched_panel_info(Xj, xi)
    ep = 1e-9
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (np.real(r1) < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0, 0, 0) if (np.real(r2) < ep) else (np.log(r2), theta1, theta2)
    P = (x * (theta1 - theta2) + d * theta2 + z * logr1 - z * logr2) / (2 * np.pi)
    dP = d
    P = (P - 0.25 * dP) if (np.real(theta1 + theta2) > np.pi) else (P + 0.75 * dP)
    a = P
    return a


def _patched_panel_linsource_stream(Xj, xi):
    t, n, x, z, d, r1, r2, theta1, theta2 = _patched_panel_info(Xj, xi)
    if np.real(theta1) < 0:
        theta1 = theta1 + 2 * np.pi
    if np.real(theta2) < 0:
        theta2 = theta2 + 2 * np.pi
    ep = 1e-9
    logr1, theta1, theta2 = (0, np.pi, np.pi) if (np.real(r1) < ep) else (np.log(r1), theta1, theta2)
    logr2, theta1, theta2 = (0, 0, 0) if (np.real(r2) < ep) else (np.log(r2), theta1, theta2)
    P1 = (0.5 / np.pi) * (x * (theta1 - theta2) + theta2 * d + z * logr1 - z * logr2)
    P2 = x * P1 + (0.5 / np.pi) * (0.5 * r2**2 * theta2 - 0.5 * r1**2 * theta1 - 0.5 * z * d)
    a = P1 - P2 / d
    b = P2 / d
    return a, b


# ---------------------------------------------------------------------------
# Patched: build_gamma  (line 627: abs(hTE) comparison)
# ---------------------------------------------------------------------------
_orig_build_gamma = _mfoil_module.build_gamma

def _patched_build_gamma(M, alpha):
    N = M.foil.N
    # Use complex-safe dtype for arrays
    is_complex = np.iscomplexobj(M.foil.x)
    dt = complex if is_complex else float

    A = np.zeros([N + 1, N + 1], dtype=dt)
    rhs = np.zeros([N + 1, 2], dtype=dt)
    t, hTE, dtdx, tcp, tdp = _patched_TE_info(M.foil.x)
    nogap = (np.real(np.abs(hTE)) < 1e-10 * M.geom.chord)

    _mfoil_module.vprint(M.param, 1, '\n <<< Solving the inviscid problem >>> \n')

    for i in range(N):
        xi = M.foil.x[:, i]
        for j in range(N - 1):
            aij, bij = _patched_panel_linvortex_stream(M.foil.x[:, [j, j + 1]], xi)
            A[i, j] += aij
            A[i, j + 1] += bij
        A[i, N] = -1

        rhs[i, :] = [-xi[1], xi[0]]
        a = _patched_panel_constsource_stream(M.foil.x[:, [N - 1, 0]], xi)
        A[i, 0] += -a * (0.5 * tcp)
        A[i, N - 1] += a * (0.5 * tcp)
        a2, b2 = _patched_panel_linvortex_stream(M.foil.x[:, [N - 1, 0]], xi)
        A[i, 0] += -(a2 + b2) * (-0.5 * tdp)
        A[i, N - 1] += (a2 + b2) * (-0.5 * tdp)

    if nogap:
        A[N - 1, :] = 0
        A[N - 1, [0, 1, 2, N - 3, N - 2, N - 1]] = [1, -2, 1, -1, 2, -1]

    A[N, 0] = 1
    A[N, N - 1] = 1

    M.isol.AIC = A
    g = np.linalg.solve(M.isol.AIC, rhs)

    M.isol.gamref = g[0:N, :]
    M.isol.gam = M.isol.gamref[:, 0] * _mfoil_module.cosd(alpha) + M.isol.gamref[:, 1] * _mfoil_module.sind(alpha)


# ---------------------------------------------------------------------------
# Patched: update_state  (lines 1845-1904: complex checks, abs(), comparisons)
# ---------------------------------------------------------------------------
_orig_update_state = _mfoil_module.update_state

def _patched_update_state(M):
    # Skip complex rejection checks entirely
    It = np.nonzero(M.vsol.turb == True)[0]
    ctmax = np.real(max(M.glob.U[2, It]))

    omega = 1.0

    for k in range(2):
        Uk = M.glob.U[k, :]; dUk = M.glob.dU[k, :]
        fmin = np.real(min(dUk / Uk))
        om = np.real(np.abs(0.5 / fmin)) if (fmin < -0.5) else 1.
        if om < omega:
            omega = om
            _mfoil_module.vprint(M.param, 3, '  th/ds decrease: omega = %.5f' % (omega))

    Uk = M.glob.U[2, :]; dUk = M.glob.dU[2, :]
    for i in range(len(Uk)):
        if (not M.vsol.turb[i]) and (np.real(Uk[i]) < .2):
            continue
        if (M.vsol.turb[i]) and (np.real(Uk[i]) < 0.1 * ctmax):
            continue
        if (np.real(Uk[i]) == 0.) or (np.real(dUk[i]) == 0.):
            continue
        if (np.real(Uk[i] + dUk[i]) < 0):
            om = np.real(0.8 * np.abs(Uk[i] / dUk[i]))
            if om < omega:
                omega = om
                _mfoil_module.vprint(M.param, 3, '  neg sa: omega = %.5f' % (omega))

    I = np.nonzero(M.vsol.turb == False)[0]
    # Skip complex check
    dumax = np.real(max(np.abs(dUk[I])))
    om = np.real(np.abs(2. / dumax)) if (dumax > 0) else 1.
    if om < omega:
        omega = om
        _mfoil_module.vprint(M.param, 3, '  amp: omega = %.5f' % (omega))

    I = np.nonzero(M.vsol.turb == True)[0]
    dumax = np.real(max(np.abs(dUk[I])))
    om = np.real(np.abs(.05 / dumax)) if (dumax > 0) else 1.
    if om < omega:
        omega = om
        _mfoil_module.vprint(M.param, 3, '  ctau: omega = %.5f' % (omega))

    dUk = M.glob.dU[3, :]
    fmax = np.real(max(np.abs(dUk) / M.oper.Vinf))
    om = 0.2 / fmax if (fmax > 0) else 1.
    if om < omega:
        omega = om
        _mfoil_module.vprint(M.param, 3, '  ue: omega = %.5f' % (omega))

    if (np.real(np.abs(M.glob.dalpha)) > 2):
        omega = min(omega, np.real(np.abs(2 / M.glob.dalpha)))

    _mfoil_module.vprint(M.param, 2, '  state update: under-relaxation = %.5f' % (omega))
    M.glob.U += omega * M.glob.dU
    M.oper.alpha += omega * M.glob.dalpha


# ---------------------------------------------------------------------------
# Patched: solve_glob  (line 1950: max(abs(ue)), sparse dtype)
# ---------------------------------------------------------------------------
_orig_solve_glob = _mfoil_module.solve_glob

def _patched_solve_glob(M):
    Nsys = M.glob.Nsys
    docl = M.oper.givencl

    ue = M.glob.U[3, :]; ds = M.glob.U[1, :]
    uemax = np.real(max(np.abs(ue)))
    for i in range(len(ue)):
        ue[i] = ue[i] if (np.real(ue[i]) > 1e-10 * uemax) else 1e-10 * uemax

    ueinv = _mfoil_module.get_ueinv(M)

    NN = 4 * Nsys + docl
    is_complex = np.iscomplexobj(M.glob.U)
    sp_dtype = complex if is_complex else float

    if (M.glob.realloc) or (type(M.glob.R_V) == list) or not (M.glob.R_V.shape == (NN, NN)):
        alloc_R_V = True
        M.glob.R_V = sparse.lil_matrix((NN, NN), dtype=sp_dtype)
    else:
        alloc_R_V = False
        M.glob.R_V *= 0.

    Ids = slice(1, 4 * Nsys, 4)
    Iue = slice(3, 4 * Nsys, 4)

    R = np.concatenate((M.glob.R, ue - (ueinv + M.vsol.ue_m @ (ds * ue))))

    M.glob.R_V[0:3 * Nsys, 0:4 * Nsys] = M.glob.R_U
    I = slice(3 * Nsys, 4 * Nsys, 1)
    M.glob.R_V[I, Iue] = sparse.identity(Nsys, dtype=sp_dtype) - M.vsol.ue_m @ np.diag(ds)
    M.glob.R_V[I, Ids] = -M.vsol.ue_m @ np.diag(ue)

    if docl:
        Rcla, Ru_alpha, Rcla_U = _mfoil_module.clalpha_residual(M)
        R = np.append(R, Rcla)
        M.glob.R_V[I, 4 * Nsys] = Ru_alpha
        M.glob.R_V[4 * Nsys, :] = Rcla_U

    if alloc_R_V:
        M.glob.R_V = M.glob.R_V.tocsr()
    dV = -sparse.linalg.spsolve(M.glob.R_V, R)

    M.glob.dU = np.reshape(dV[0:4 * Nsys], (4, Nsys), order='F')
    if docl:
        M.glob.dalpha = dV[-1]


# ---------------------------------------------------------------------------
# Patched: march_amplification  (lines 2596, 2600, 2612: abs(), complex check)
# ---------------------------------------------------------------------------
_orig_march_amplification = _mfoil_module.march_amplification

def _patched_march_amplification(M, si):
    Is = M.vsol.Is[si]
    N = len(Is)
    param = _mfoil_module.build_param(M, si)
    U = M.glob.U[:, Is]
    turb = M.vsol.turb[Is]

    U[2, 0] = 0.
    param.turb, param.wake = False, False
    i = 1
    while i < N:
        U1, U2 = U[:, i - 1], U[:, i].copy()
        if turb[i]:
            U2[2] = U1[2] * 1.01
        dx = M.isol.xi[Is[i]] - M.isol.xi[Is[i - 1]]

        nNewton = 20
        for iNewton in range(nNewton):
            damp1, damp1_U1 = _mfoil_module.get_damp(U1, param)
            damp2, damp2_U2 = _mfoil_module.get_damp(U2, param)
            damp, damp_U = _mfoil_module.upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2)

            Ramp = U2[2] - U1[2] - damp * dx

            if iNewton > 11:
                _mfoil_module.vprint(param, 3, 'i=%d, iNewton=%d, sa = [%.5e, %.5e], damp = %.5e, Ramp = %.5e' % (
                    i, iNewton, np.real(U1[2]), np.real(U2[2]), np.real(damp), np.real(Ramp)))

            if np.real(np.abs(Ramp)) < 1e-12:
                break
            Ramp_U = np.r_[0, 0, -1, 0, 0, 0, 1, 0] - damp_U * dx
            dU = -Ramp / Ramp_U[6]
            omega = 1; dmax = 0.5 * (1.01 - iNewton / nNewton)
            if np.real(np.abs(dU)) > dmax:
                omega = dmax / np.abs(dU)
            U2[2] += omega * dU

        if iNewton >= nNewton - 1:
            _mfoil_module.vprint(param, 1, 'march amp Newton unconverged!')

        if np.real(U2[2]) > param.ncrit:
            _mfoil_module.vprint(param, 2, '  march_amplification (si,i=%d,%d): %.5e is above critical.' % (si, i, np.real(U2[2])))
            break
        else:
            M.glob.U[2, Is[i]] = U2[2]
            U[2, i] = U2[2]
            # Skip complex check

        i += 1

    return (i - 1)


# ---------------------------------------------------------------------------
# Patched: residual_transition  (line 2711: complex check; lines 2661-2662: abs())
# ---------------------------------------------------------------------------
_orig_residual_transition = _mfoil_module.residual_transition

def _patched_residual_transition(M, param, x, U, Aux):
    # This function is large but only needs a few lines changed.
    # We patch the inner Newton loop's abs() and the final complex check.
    # Rather than duplicating the entire function, we temporarily override
    # the module-level abs and np.sign, then call the original.
    #
    # Actually, the complex rejection is a simple check at line 2711.
    # The abs() calls at 2661-2662 are for Newton damping of xt.
    # Let's do a full replacement since the function is self-contained.

    I1 = slice(0, 4)
    I2 = slice(4, 8)
    x1, x2 = x[0], x[1]
    U1, U2 = U[:, 0], U[:, 1]
    dx = x2 - x1
    ncrit = param.ncrit
    sa = U[:, :][2, :]

    Z = np.zeros(4)

    # determine transition location xt
    xt = x1 + 0.5 * dx  # initial guess
    nNewton = 20
    for iNewton in range(nNewton):
        w2 = (xt - x1) / dx; w1 = 1. - w2
        Ut = w1 * U1 + w2 * U2
        Ut_xt = (U2 - U1) / dx
        Ut[2] = ncrit; Ut_xt[2] = 0.

        param_save = param.turb
        param.turb = False
        damp1, damp1_U1 = _mfoil_module.get_damp(U1, param)
        dampt, dampt_Ut = _mfoil_module.get_damp(Ut, param); dampt_Ut[2] = 0.
        param.turb = param_save

        Rxt = ncrit - sa[0] - 0.5 * (xt - x1) * (damp1 + dampt)
        Rxt_xt = -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * np.dot(dampt_Ut, Ut_xt)
        dxt = -Rxt / Rxt_xt
        _mfoil_module.vprint(param, 4, '   Transition: iNewton,Rxt,xt = %d,%.5e,%.5e' % (iNewton, np.real(Rxt), np.real(xt)))
        dmax = 0.2 * dx * (1.1 - iNewton / nNewton)
        if np.real(np.abs(dxt)) > np.real(dmax):
            dxt = dxt * dmax / np.abs(dxt)
        if np.real(np.abs(Rxt)) < 1e-10:
            break
        if iNewton < nNewton:
            xt += dxt

    if iNewton >= nNewton:
        _mfoil_module.vprint(param, 1, 'Transition location calculation failed.')
    M.vsol.xt = xt

    # prepare for xt linearizations
    Rxt_U = -0.5 * (xt - x1) * np.concatenate((damp1_U1 + dampt_Ut * w1, dampt_Ut * w2)); Rxt_U[2] -= 1.
    Ut_x1 = (U2 - U1) * (w2 - 1) / dx; Ut_x2 = (U2 - U1) * (-w2) / dx
    Ut_x1[2] = 0; Ut_x2[2] = 0
    Rxt_x1 = 0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * np.dot(dampt_Ut, Ut_x1)
    Rxt_x2 = -0.5 * (xt - x1) * np.dot(dampt_Ut, Ut_x2)

    xt_U = -Rxt_U / Rxt_xt; xt_U1 = xt_U[I1]; xt_U2 = xt_U[I2]
    xt_x1 = -Rxt_x1 / Rxt_xt; xt_x2 = -Rxt_x2 / Rxt_xt

    Ut_x1 += Ut_xt * xt_x1
    Ut_x2 += Ut_xt * xt_x2

    Ut_U1 = w1 * np.eye(4) + np.outer((U2 - U1), xt_U1) / dx
    Ut_U2 = w2 * np.eye(4) + np.outer((U2 - U1), xt_U2) / dx

    Utl = Ut.copy(); Utl_U1 = Ut_U1.copy(); Utl_U2 = Ut_U2.copy(); Utl_x1 = Ut_x1.copy(); Utl_x2 = Ut_x2.copy()
    Utl[2] = ncrit; Utl_U1[2, :] = Z; Utl_U2[2, :] = Z; Utl_x1[2] = 0; Utl_x2[2] = 0
    Utt = Ut.copy(); Utt_U1 = Ut_U1.copy(); Utt_U2 = Ut_U2.copy(); Utt_x1 = Ut_x1.copy(); Utt_x2 = Ut_x2.copy()

    param2 = _mfoil_module.build_param(M, 0)

    param2.turb = True
    cttr, cttr_Ut = _mfoil_module.get_cttr(Ut, param2)
    Utt[2] = cttr; Utt_U1[2, :] = np.dot(cttr_Ut, Ut_U1); Utt_U2[2, :] = np.dot(cttr_Ut, Ut_U2)
    Utt_x1[2] = np.dot(cttr_Ut, Ut_x1); Utt_x2[2] = np.dot(cttr_Ut, Ut_x2)

    param2.turb = False
    Rl, Rl_U, Rl_x = _mfoil_module.residual_station(param2, np.r_[x1, xt], np.stack((U1, Utl), axis=-1), Aux)
    Rl_U1 = Rl_U[:, I1]; Rl_Utl = Rl_U[:, I2]
    param2.turb = True
    Rt, Rt_U, Rt_x = _mfoil_module.residual_station(param2, np.r_[xt, x2], np.stack((Utt, U2), axis=-1), Aux)
    Rt_Utt = Rt_U[:, I1]; Rt_U2 = Rt_U[:, I2]

    R = Rl + Rt
    # Skip complex rejection check (was: if any(np.imag(R) != 0): raise)
    R_U1 = Rl_U1 + np.dot(Rl_Utl, Utl_U1) + np.outer(Rl_x[:, 1], xt_U1) + np.dot(Rt_Utt, Utt_U1) + np.outer(Rt_x[:, 0], xt_U1)
    R_U2 = np.dot(Rl_Utl, Utl_U2) + np.outer(Rl_x[:, 1], xt_U2) + np.dot(Rt_Utt, Utt_U2) + Rt_U2 + np.outer(Rt_x[:, 0], xt_U2)
    R_U = np.hstack((R_U1, R_U2))
    R_x = np.stack((Rl_x[:, 0] + Rl_x[:, 1] * xt_x1 + Rt_x[:, 0] * xt_x1 + np.dot(Rl_Utl, Utl_x1) + np.dot(Rt_Utt, Utt_x1),
                     Rl_x[:, 1] * xt_x2 + Rt_x[:, 1] + Rt_x[:, 0] * xt_x2 + np.dot(Rl_Utl, Utl_x2) + np.dot(Rt_Utt, Utt_x2)), axis=-1)

    return R, R_U, R_x


# ---------------------------------------------------------------------------
# Patched: get_upw  (line 2944: np.sign(aa))
# ---------------------------------------------------------------------------
_orig_get_upw = _mfoil_module.get_upw

def _patched_get_upw(U1, U2, param):
    Hk1, Hk1_U1 = _mfoil_module.get_Hk(U1, param)
    Hk2, Hk2_U2 = _mfoil_module.get_Hk(U2, param)
    Z = np.zeros(Hk1_U1.shape)
    Hut = 1.
    C = 1. if param.wake else 5.
    Huc = C * Hut / Hk2**2
    Huc_U = np.concatenate((Z, -2 * Huc / Hk2 * Hk2_U2))
    aa = (Hk2 - 1.) / (Hk1 - 1.)
    sga = np.sign(np.real(aa))  # complex-safe: use real part for sign
    la = np.log(sga * aa)
    la_U = np.concatenate((-1. / (Hk1 - 1.) * Hk1_U1, 1. / (Hk2 - 1.) * Hk2_U2))
    Hls = la**2; Hls_U = 2 * la * la_U
    if np.real(Hls) > 15:
        Hls, Hls_U = 15, Hls_U * 0.
    upw = 1. - 0.5 * np.exp(-Hls * Huc)
    upw_U = -0.5 * np.exp(-Hls * Huc) * (-Hls_U * Huc - Hls * Huc_U)
    return upw, upw_U


# ---------------------------------------------------------------------------
# Patched: BL init marching  (lines 2338-2430: abs(), max() on complex)
# ---------------------------------------------------------------------------
_orig_init_boundary_layer = _mfoil_module.init_boundary_layer

def _patched_init_boundary_layer(M):
    # Delegates to original. All builtins (abs, max, min) and numpy array creation
    # are handled by the complex_step_mode() context manager at the call site.
    _orig_init_boundary_layer(M)


# ---------------------------------------------------------------------------
# Patched: calc_force  (line 422: np.zeros(N) must be complex-compatible)
# ---------------------------------------------------------------------------
_orig_calc_force = _mfoil_module.calc_force

def _patched_calc_force(M):
    # The only issue is that cl_ue = np.zeros(N) is float but cp_ue can be complex.
    # Patch: force U row 3 dtype to propagate through by ensuring float arrays
    # are allocated with complex dtype when foil.x is complex.
    is_complex = np.iscomplexobj(M.foil.x) or np.iscomplexobj(M.glob.U)

    chord = M.geom.chord; xref = M.geom.xref
    Vinf = M.param.Vinf; rho = M.oper.rho; alpha = M.oper.alpha
    qinf = 0.5 * rho * Vinf**2
    N = M.foil.N

    ue = M.glob.U[3, :] if M.oper.viscous else _mfoil_module.get_ueinv(M)
    cp, cp_ue = _mfoil_module.get_cp(ue, M.param); M.post.cp = cp
    M.post.cpi, cpi_ue = _mfoil_module.get_cp(_mfoil_module.get_ueinv(M), M.param)

    # Recheck: ue or cp may be complex even if foil.x/U are real (e.g. complex alpha)
    is_complex = is_complex or np.iscomplexobj(ue) or np.iscomplexobj(cp)
    dt = complex if is_complex else float
    cl, cl_ue, cl_alpha, cm, cdpi = dt(0), np.zeros(N, dtype=dt), dt(0), dt(0), dt(0)
    for i0 in range(1, N + 1):
        i, ip = (0, N - 1) if (i0 == N) else (i0, i0 - 1)
        x1, x2 = M.foil.x[:, ip], M.foil.x[:, i]
        dxv = x2 - x1; dx1 = x1 - xref; dx2 = x2 - xref
        dx1nds = dxv[0] * dx1[0] + dxv[1] * dx1[1]
        dx2nds = dxv[0] * dx2[0] + dxv[1] * dx2[1]
        dx = -dxv[0] * _mfoil_module.cosd(alpha) - dxv[1] * _mfoil_module.sind(alpha)
        dz = dxv[1] * _mfoil_module.cosd(alpha) - dxv[0] * _mfoil_module.sind(alpha)
        cp1, cp2 = cp[ip], cp[i]; cpbar = 0.5 * (cp1 + cp2)
        cl = cl + dx * cpbar
        I = [ip, i]; cl_ue[I] += dx * 0.5 * cp_ue[I]
        cl_alpha += cpbar * (_mfoil_module.sind(alpha) * dxv[0] - _mfoil_module.cosd(alpha) * dxv[1]) * np.pi / 180
        cm += cp1 * dx1nds / 3 + cp1 * dx2nds / 6 + cp2 * dx1nds / 6 + cp2 * dx2nds / 3
        cdpi = cdpi + dz * cpbar
        cl /= chord; cm /= chord**2; cdpi /= chord
        M.post.cl, M.post.cl_ue, M.post.cl_alpha = cl, cl_ue, cl_alpha
        M.post.cm, M.post.cdpi = cm, cdpi

    cd, cdf = dt(0), dt(0)
    if M.oper.viscous:
        iw = M.vsol.Is[2][-1]
        U = M.glob.U[:, iw]
        H, H_U = _mfoil_module.get_H(U)
        uk, uk_ue = _mfoil_module.get_uk(U[3], M.param)
        cd = 2.0 * U[0] * (uk / Vinf)**((5 + H) / 2.)
        M.post.cd_U = 2.0 * U[0] * (uk / Vinf)**((5 + H) / 2.) * np.log(uk / Vinf) * 0.5 * H_U
        M.post.cd_U[0] += 2.0 * (uk / Vinf)**((5 + H) / 2.)
        M.post.cd_U[3] += 2.0 * U[0] * (5 + H) / 2. * (uk / Vinf)**((3 + H) / 2.) * (1. / Vinf) * uk_ue

        Df = dt(0)
        for si in range(2):
            Is = M.vsol.Is[si]
            param = _mfoil_module.build_param(M, si)
            _mfoil_module.station_param(M, param, Is[0])
            cf1 = dt(0)
            ue1 = dt(0)
            rho1 = rho
            x1 = M.isol.xstag
            for ii in range(len(Is)):
                _mfoil_module.station_param(M, param, Is[ii])
                cf2, cf2_U = _mfoil_module.get_cf(M.glob.U[:, Is[ii]], param)
                ue2, ue2_ue = _mfoil_module.get_uk(M.glob.U[3, Is[ii]], param)
                rho2, rho2_U = _mfoil_module.get_rho(M.glob.U[:, Is[ii]], param)
                x2 = M.foil.x[:, Is[ii]]; dxv = x2 - x1
                dxx = dxv[0] * _mfoil_module.cosd(alpha) + dxv[1] * _mfoil_module.sind(alpha)
                Df += 0.25 * (rho1 * cf1 * ue1**2 + rho2 * cf2 * ue2**2) * dxx
                cf1 = cf2; ue1 = ue2; x1 = x2; rho1 = rho2
        cdf = Df / (qinf * chord)

    M.post.cd, M.post.cdf, M.post.cdp = cd, cdf, cd - cdf

    s = '  alpha=%.2fdeg, cl=%.6f, cm=%.6f, cdpi=%.6f, cd=%.6f, cdf=%.6f, cdp=%.6f' % (
        np.real(M.oper.alpha), np.real(M.post.cl), np.real(M.post.cm),
        np.real(M.post.cdpi), np.real(M.post.cd), np.real(M.post.cdf), np.real(M.post.cdp))
    _mfoil_module.vprint(M.param, 1, s)


# ---------------------------------------------------------------------------
# Patched: inviscid_velocity (line 688: np.zeros(2) must be complex)
# ---------------------------------------------------------------------------
_orig_inviscid_velocity = _mfoil_module.inviscid_velocity

def _patched_inviscid_velocity(X, G, Vinf, alpha, x, dolin):
    N = X.shape[1]
    is_complex = np.iscomplexobj(G) or np.iscomplexobj(X) or np.iscomplexobj(x)
    dt = complex if is_complex else float
    V = np.zeros(2, dtype=dt)
    if dolin:
        V_G = np.zeros([2, N], dtype=dt)
    t, hTE, dtdx, tcp, tdp = _mfoil_module.TE_info(X)
    for j in range(N - 1):
        a, b = _mfoil_module.panel_linvortex_velocity(X[:, [j, j + 1]], x, None, False)
        V += a * G[j] + b * G[j + 1]
        if dolin:
            V_G[:, j] += a
            V_G[:, j + 1] += b
    a = _mfoil_module.panel_constsource_velocity(X[:, [N - 1, 0]], x, None)
    f1, f2 = a * (-0.5 * tcp), a * 0.5 * tcp
    V += f1 * G[0] + f2 * G[N - 1]
    if dolin:
        V_G[:, 0] += f1
        V_G[:, N - 1] += f2
    a2, b2 = _mfoil_module.panel_linvortex_velocity(X[:, [N - 1, 0]], x, None, False)
    f1, f2 = (a2 + b2) * (0.5 * tdp), (a2 + b2) * (-0.5 * tdp)
    V += f1 * G[0] + f2 * G[N - 1]
    if dolin:
        V_G[:, 0] += f1
        V_G[:, N - 1] += f2
    V += Vinf * np.array([_mfoil_module.cosd(alpha), _mfoil_module.sind(alpha)])
    if dolin:
        return V, V_G
    return V


# ---------------------------------------------------------------------------
# Patched: get_ueinv (returns complex when gam is complex)
# ---------------------------------------------------------------------------
_orig_get_ueinv = _mfoil_module.get_ueinv

def _patched_get_ueinv(M):
    # get_ueinv calls inviscid_velocity in a loop; with our patched version it works.
    # But it also creates np.zeros arrays. Let's just call the original since
    # inviscid_velocity is already patched.
    return _orig_get_ueinv(M)


# ---------------------------------------------------------------------------
# Apply all monkey-patches
# ---------------------------------------------------------------------------
_mfoil_module.TE_info = _patched_TE_info
_mfoil_module.panel_info = _patched_panel_info
_mfoil_module.panel_linvortex_stream = _patched_panel_linvortex_stream
_mfoil_module.panel_constsource_velocity = _patched_panel_constsource_velocity
_mfoil_module.panel_constsource_stream = _patched_panel_constsource_stream
_mfoil_module.panel_linsource_stream = _patched_panel_linsource_stream
_mfoil_module.build_gamma = _patched_build_gamma
_mfoil_module.update_state = _patched_update_state
_mfoil_module.solve_glob = _patched_solve_glob
_mfoil_module.march_amplification = _patched_march_amplification
_mfoil_module.residual_transition = _patched_residual_transition
_mfoil_module.get_upw = _patched_get_upw
_mfoil_module.init_boundary_layer = _patched_init_boundary_layer
_mfoil_module.calc_force = _patched_calc_force
_mfoil_module.inviscid_velocity = _patched_inviscid_velocity

# ---------------------------------------------------------------------------
# Patched: stagpoint_find  (line 806: abs() destroys complex imaginary part)
# ---------------------------------------------------------------------------
_orig_stagpoint_find = _mfoil_module.stagpoint_find

def _patched_stagpoint_find(M):
    """Complex-step safe stagpoint_find.

    The original uses abs(M.foil.s - sstag) to compute xi, which destroys
    imaginary parts when foil.s is complex. Replace with sign-aware computation
    using the known stagnation panel index j.
    """
    assert len(M.isol.gam) > 0, 'No inviscid solution'
    N = M.foil.N
    for j in range(N):
        if np.real(M.isol.gam[j]) > 0:
            break
    assert j < N - 1, 'no stagnation point'
    I = [j - 1, j]
    G = M.isol.gam[I]
    S = M.foil.s[I]
    M.isol.Istag = I
    den = G[1] - G[0]
    w1 = G[1] / den
    w2 = -G[0] / den
    M.isol.sstag = w1 * S[0] + w2 * S[1]
    M.isol.xstag = M.foil.x[:, j - 1] * w1 + M.foil.x[:, j] * w2
    M.isol.sstag_g = np.array([G[1] * (S[0] - S[1]) / (den * den), -G[0] * (S[0] - S[1]) / (den * den)])
    sgnue = -1 * np.ones(N)
    for i in range(j, N):
        sgnue[i] = 1
    M.isol.sgnue = sgnue

    # Complex-step safe xi: avoid abs() which strips imaginary part.
    # For nodes before stag (i < j): s < sstag, so xi = sstag - s = -(s - sstag)
    # For nodes at/after stag (i >= j): s >= sstag, so xi = s - sstag
    diff = M.foil.s - M.isol.sstag
    xi_foil = diff.copy()
    xi_foil[:j] = -diff[:j]
    M.isol.xi = np.concatenate((xi_foil, M.wake.s - M.isol.sstag))

_mfoil_module.stagpoint_find = _patched_stagpoint_find

# Patch init_thermo: comparison Minf > 0 with complex Ma
_orig_init_thermo = _mfoil_module.init_thermo

def _patched_init_thermo(M):
    g = M.param.gam; gmi = g - 1
    rhoinf = M.oper.rho
    Vinf = M.oper.Vinf; M.param.Vinf = Vinf
    M.param.muinf = rhoinf * Vinf * M.geom.chord / M.oper.Re
    Minf = M.oper.Ma; M.param.Minf = Minf
    if np.real(Minf) > 0:
        M.param.KTb = np.sqrt(1 - Minf**2)
        M.param.KTl = Minf**2 / (1 + M.param.KTb)**2
        M.param.H0 = (1 + 0.5 * gmi * Minf**2) * Vinf**2 / (gmi * Minf**2)
        Tr = 1 - 0.5 * Vinf**2 / M.param.H0
        finf = Tr**1.5 * (1 + M.param.Tsrat) / (Tr + M.param.Tsrat)
        M.param.cps = 2 / (g * Minf**2) * (((1 + 0.5 * gmi * Minf**2) / (1 + 0.5 * gmi))**(g / gmi) - 1)
    else:
        finf = 1
    M.param.mu0 = M.param.muinf / finf
    M.param.rho0 = rhoinf * (1 + 0.5 * gmi * Minf**2)**(1 / gmi)

_mfoil_module.init_thermo = _patched_init_thermo

# Patch BL property functions that compare param.Minf > 0
# Strategy: wrap each function to temporarily set param.Minf to its real part
# for comparison, then restore the complex value for computation.
# Actually simpler: each function only uses `if (param.Minf > 0):` as a branch.
# We can wrap the param object to make Minf comparison use real part.
# Most efficient: just patch the specific functions.

def _make_minf_safe(orig_func):
    """Wrap a function that does `if (param.Minf > 0):` to use np.real() for the check."""
    import functools
    @functools.wraps(orig_func)
    def wrapper(*args, **kwargs):
        # Find param argument (usually last positional arg)
        param = None
        for a in args:
            if hasattr(a, 'Minf'):
                param = a
                break
        if param is not None and isinstance(param.Minf, (complex, np.complexfloating)):
            # Temporarily replace Minf with a special wrapper
            orig_minf = param.Minf
            class _CmpFloat(complex):
                """Complex that compares using real part."""
                def __gt__(self, other): return self.real > other
                def __lt__(self, other): return self.real < other
                def __ge__(self, other): return self.real >= other
                def __le__(self, other): return self.real <= other
            param.Minf = _CmpFloat(orig_minf)
            try:
                return orig_func(*args, **kwargs)
            finally:
                param.Minf = orig_minf
        return orig_func(*args, **kwargs)
    return wrapper

for _fname in ['get_cp', 'get_uk', 'get_Mach2', 'get_Hss', 'get_rho', 'get_uq', 'get_Hk', 'get_Ret']:
    _orig = getattr(_mfoil_module, _fname)
    setattr(_mfoil_module, _fname, _make_minf_safe(_orig))

# Also patch the check at line 289 in solve_inviscid (M.oper.Ma > 0)
# This is inside solve_inviscid, but it's a check_cpstar_stall function
_orig_solve_inviscid = _mfoil_module.solve_inviscid

def _patched_solve_inviscid(M):
    # Temporarily make M.oper.Ma comparable
    orig_ma = M.oper.Ma
    if isinstance(orig_ma, (complex, np.complexfloating)):
        class _CmpFloat(complex):
            def __gt__(self, other): return self.real > other
            def __lt__(self, other): return self.real < other
            def __ge__(self, other): return self.real >= other
            def __le__(self, other): return self.real <= other
        M.oper.Ma = _CmpFloat(orig_ma)
    try:
        _orig_solve_inviscid(M)
    finally:
        M.oper.Ma = orig_ma

_mfoil_module.solve_inviscid = _patched_solve_inviscid

# Patch norm2 to use analytic norm (sum of squares) instead of |x_i|^2
_orig_norm2 = _mfoil_module.norm2

def _patched_norm2(x):
    """Analytic vector 2-norm: sqrt(sum(x_i^2)) instead of sqrt(sum(|x_i|^2))."""
    if np.iscomplexobj(x):
        return np.sqrt(np.dot(x, x))
    return _orig_norm2(x)

_mfoil_module.norm2 = _patched_norm2

# Patch dist similarly: dist(a,b) = sqrt(a^2 + b^2) — already analytic as written,
# but let's verify it uses the correct formula
_orig_dist = _mfoil_module.dist

def _patched_dist(a, b):
    """Analytic distance: sqrt(a^2 + b^2)."""
    return np.sqrt(a**2 + b**2)

_mfoil_module.dist = _patched_dist

# Patch vprint to handle complex values in format strings
_orig_vprint = _mfoil_module.vprint

def _patched_vprint(param, level, msg, *args):
    try:
        _orig_vprint(param, level, msg, *args)
    except (TypeError, np.ComplexWarning):
        # If formatting fails due to complex, just skip the print
        pass

_mfoil_module.vprint = _patched_vprint

# Also patch the module-level atan2 used by other panel functions we didn't replace
_mfoil_module.atan2 = _cs_atan2


# ---------------------------------------------------------------------------
# Complex-step context manager: forces np.zeros to create complex arrays
# inside the mfoil module when solving with complex perturbations.
# ---------------------------------------------------------------------------
import contextlib

_np_zeros_orig = np.zeros
_np_ones_orig = np.ones
_np_empty_orig = np.empty

def _complex_zeros(shape, dtype=None, order='C'):
    if dtype is None or dtype == float or dtype == np.float64:
        dtype = complex
    return _np_zeros_orig(shape, dtype=dtype, order=order)

def _complex_ones(shape, dtype=None, order='C'):
    if dtype is None or dtype == float or dtype == np.float64:
        dtype = complex
    return _np_ones_orig(shape, dtype=dtype, order=order)

def _complex_empty(shape, dtype=None, order='C'):
    if dtype is None or dtype == float or dtype == np.float64:
        dtype = complex
    return _np_empty_orig(shape, dtype=dtype, order=order)

@contextlib.contextmanager
def complex_step_mode():
    """Context manager that makes numpy array creation in mfoil produce complex arrays."""
    import builtins
    orig_abs = builtins.abs
    orig_max = builtins.max
    orig_min = builtins.min

    def _cs_abs(x):
        if isinstance(x, (complex, np.complexfloating)):
            return np.abs(x)
        try:
            if np.iscomplexobj(x):
                return np.abs(x)
        except (TypeError, ValueError):
            pass
        return orig_abs(x)

    def _cs_max(*args, **kwargs):
        if len(args) == 1:
            it = args[0]
            try:
                arr = np.asarray(it) if not isinstance(it, np.ndarray) else it
                if np.iscomplexobj(arr):
                    idx = np.argmax(np.real(arr))
                    return arr.flat[idx]
            except (TypeError, ValueError):
                pass
            return orig_max(*args, **kwargs)
        elif len(args) == 2:
            a, b = args
            ar = np.real(a) if isinstance(a, (complex, np.complexfloating)) else a
            br = np.real(b) if isinstance(b, (complex, np.complexfloating)) else b
            return a if ar >= br else b
        return orig_max(*args, **kwargs)

    def _cs_min(*args, **kwargs):
        if len(args) == 1:
            it = args[0]
            try:
                arr = np.asarray(it) if not isinstance(it, np.ndarray) else it
                if np.iscomplexobj(arr):
                    idx = np.argmin(np.real(arr))
                    return arr.flat[idx]
            except (TypeError, ValueError):
                pass
            return orig_min(*args, **kwargs)
        elif len(args) == 2:
            a, b = args
            ar = np.real(a) if isinstance(a, (complex, np.complexfloating)) else a
            br = np.real(b) if isinstance(b, (complex, np.complexfloating)) else b
            return a if ar <= br else b
        return orig_min(*args, **kwargs)

    import warnings
    # Suppress ComplexWarning from format strings (harmless: only affects diagnostic prints)
    _warn_ctx = warnings.catch_warnings()
    _warn_ctx.__enter__()
    # np.ComplexWarning was removed in numpy 2.x; ignore by name fallback.
    _cw = getattr(np, 'ComplexWarning', None)
    if _cw is None:
        try:
            from numpy.exceptions import ComplexWarning as _cw  # numpy >= 1.25
        except ImportError:
            _cw = Warning  # very loose fallback
    warnings.simplefilter('ignore', _cw)

    # Patch numpy globally (mfoil uses `import numpy as np`, which references the module)
    np.zeros = _complex_zeros
    np.ones = _complex_ones
    np.empty = _complex_empty
    builtins.abs = _cs_abs
    builtins.max = _cs_max
    builtins.min = _cs_min

    # Patch scipy.sparse to create complex sparse matrices
    from scipy import sparse as _sp
    _orig_lil_matrix = _sp.lil_matrix
    _orig_identity = _sp.identity

    class _complex_lil_matrix(_orig_lil_matrix):
        def __init__(self, *args, dtype=None, **kwargs):
            if dtype is None or dtype == float or dtype == np.float64:
                dtype = complex
            super().__init__(*args, dtype=dtype, **kwargs)

    def _complex_identity(n, dtype=None, format=None):
        if dtype is None or dtype == float or dtype == np.float64:
            dtype = complex
        return _orig_identity(n, dtype=dtype, format=format)

    _sp.lil_matrix = _complex_lil_matrix
    _sp.identity = _complex_identity

    try:
        yield
    finally:
        np.zeros = _np_zeros_orig
        np.ones = _np_ones_orig
        np.empty = _np_empty_orig
        builtins.abs = orig_abs
        builtins.max = orig_max
        builtins.min = orig_min
        _sp.lil_matrix = _orig_lil_matrix
        _sp.identity = _orig_identity
        _warn_ctx.__exit__(None, None, None)


# Re-export everything from mfoil
from mfoil import *
