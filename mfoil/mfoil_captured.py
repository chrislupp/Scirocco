"""Captured laminar-turbulent transition for mfoil (prototype).

The standard mfoil boundary layer stores one third state per node: the
amplification factor n at laminar stations, or the root shear-stress
coefficient ct at turbulent ones. Transition is a cut cell [x1, xt, x2], so
the discrete equations switch type whenever xt crosses a node. Outputs are
continuous, but their derivatives jump: on NACA 2412 at Re 1e6, dCd/dalpha
saw-tooths by +-40% each time transition moves one panel.

Here both variables live at every node. The state is [th, ds, ct, ue, n]
(ct stays in row 2, n is row 4), and every interval carries the same four
equations, discretized exactly as mfoil.residual_station:

  momentum, shape  with closures Hs, cf, cDi blended per station,
                   f = (1 - g) f_lam + g f_turb,  g = S((n - nc) / wc)
  amplification    dn/dx = (1 - gs) damp + gs smax(damp, C / th),
                   gs = S((n - ncrit) / ws)
  shear lag        mfoil's lag equation for ct, with equilibrium target
                   (1 - r) cttr + r cteq and pressure-gradient term times r,
                   r = S((n - ncrit + 2.5 ws) / ws)

S(t) = (1 + tanh t) / 2; in the wake g = r = 1 and n2 = n1.

- Up to ncrit the equations are the standard laminar ones, so onset is where
  the cut cell puts transition. ct tracks cttr (XFOIL's initial value) and
  is released into the lag equation just before onset.
- After onset n grows at least at C / th, whatever the profile. Otherwise
  the blend would lower Hk, the laminar rate would collapse, and n would
  stall; that feedback produced multiple solutions in earlier versions. The
  smooth max keeps the faster laminar rate where it is larger (laminar
  separation bubbles); a hard floor alone over-predicted bubble drag by up
  to 22%.
- The closures blend as n crosses the band of g: 1% at ncrit + b ws, and
  10-90% over about L momentum thicknesses (C = 2.2 wc / L). b = 0 starts
  intermittency at onset, as in Dhawan & Narasimha (1958).
- Blending closures, not residuals, keeps each equation a single ODE whose
  discretization telescopes. Blending residuals left a node-dependent source
  term and was rougher than the cut cell.
- At the stagnation similarity station the lag equation degenerates
  (dx ~ 1e-6), so ct is pinned there algebraically, to the first node's ct
  rather than an extrapolation (which can go negative).

Parameters (on M.param): capt_L (default 60, momentum thicknesses; calibrated
on Eppler 387 drag, see scirocco/validation/captured_transition),
capt_b (0), capt_ws (1), capt_wc (1). capt_Lmodel = 'narasimha' takes the
length from the Dhawan-Narasimha spot-formation-rate correlation instead
(see transition_rate); it lowers Cd by 3-11% at Re >= 1e6, consistent with
double-counting the ramp the shear-lag equation already models, so it is
not the default. The zone should span several panels: at Re 3e6 the
dCd/dalpha ripple is 4.3% at 160 panels and 1.2% at 320.

Forced transition (M.param.xtrip = [lower, upper], located by
mfoil_cs.find_trip): in the interval that holds the trip, a fraction f along
it, n jumps to at least ncrit,

  npre  = n1 + q1 f dx                    laminar rate up to the trip
  npost = npre + ws softplus((ncrit - npre) / ws)
  n2    = npost + (q1* + q2) / 2 (1 - f) dx

with q1* the rate at station 1 with its growth switch at npost. The trip
acts as onset does: the closures blend over capt_L momentum thicknesses
behind it, and where free transition came first the jump vanishes.

Usage:
    M = mfoil.mfoil(naca='2412', npanel=160)
    M.setoper(alpha=4, Re=1e6)
    solve_captured(M)          # standard solve, then captured Newton
"""
import copy

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as splinalg

import mfoil as mf

# State rows
TH, DS, CT, UE, AMP = 0, 1, 2, 3, 4
LAM_ROWS = [TH, DS, AMP, UE]  # 4-state view for laminar residuals
TURB_ROWS = [TH, DS, CT, UE]  # 4-state view for turbulent residuals


def _cfg(param):
    """Switch layout: centers are offsets from ncrit, widths in n.

    G : onset switch gs; the n rate goes from the laminar rate to Gfull.
    cl: closure intermittency g, 1% at ncrit + b ws (param.capt_b, default 0).
    ct: shear-lag release, ahead of onset.
    The post-onset rate Gfull = 2.2 wc / (L th) makes n cross the 10-90%
    band of g in L momentum thicknesses (param.capt_L).
    """
    L = getattr(param, 'capt_L', 60.0)
    ws = getattr(param, 'capt_ws', 1.0)
    wc = getattr(param, 'capt_wc', 1.0)
    b = getattr(param, 'capt_b', 0.0)
    d = dict(cl=(b*ws + 2.3*wc, wc), G=(0.0, ws), ct=(-2.5*ws, ws),
             C=2.2*wc/L)
    d.update(getattr(param, 'capt_cfg', {}))
    return d


def _switch(n, center, w):
    t = np.tanh((n - center)/w)
    return 0.5*(1 + t), 0.5*(1 - t**2)/w


def gamma(n, param):
    """g(n): intermittency of the closures, and d/dn."""
    off, w = _cfg(param)['cl']
    return _switch(n, param.ncrit + off, w)


def growth_switch(n, param):
    off, w = _cfg(param)['G']
    return _switch(n, param.ncrit + off, w)


def release(n, param):
    off, w = _cfg(param)['ct']
    return _switch(n, param.ncrit + off, w)


def _C(param):
    return _cfg(param)['C']


def _smax(a, a_U, b, b_U, eps, eps_U):
    """Smooth max(a, b) with rounding scale eps, and its linearization."""
    r = np.sqrt((a - b)**2 + eps**2)
    m = 0.5*(a + b + r)
    m_U = 0.5*(a_U + b_U + ((a - b)*(a_U - b_U) + eps*eps_U)/r)
    return m, m_U


def _rate(d, d_U, s, s_n, C, th, k, E):
    """Post-onset n rate at one station: (1 - s) d + s smax(d, C / th).

    The floor C / th carries n through the closure blend without feedback
    from the profile; the max keeps fast growth where the laminar rate is
    larger (laminar separation bubbles). d_U is 5-state; k is unused.
    """
    c = C/th; c_U = -C/th**2*E[TH]
    ke = 0.1
    m, m_U = _smax(d, d_U, c, c_U, ke*c, ke*c_U)
    q = (1 - s)*d + s*m
    q_U = (1 - s)*d_U + s*m_U + (m - d)*s_n*E[AMP]
    r = np.sqrt((d - c)**2 + (ke*c)**2)
    q_C = s*0.5*(1 - (d - c)/r + ke*ke*c/r)/th
    return q, q_U, q_C


def _trip_jump(npre, param):
    """J = ws softplus((ncrit - npre) / ws) and -dJ/dnpre."""
    ws = getattr(param, 'capt_ws', 1.0)
    t = (param.ncrit - npre)/ws
    if np.real(t) > 30:
        return ws*t, 1.0
    e = np.exp(t)
    return ws*np.log1p(e), e/(1 + e)


def _find_trip(M, si):
    """(station, frac) of the trip on surface si; station < 0 if none."""
    if si > 1 or getattr(M.param, 'xtrip', None) is None:
        return -1, 1.
    import mfoil_cs
    i, frac, _ = mfoil_cs.find_trip(M, si)
    return i, frac


def _dN(L):
    """Increase of the spot formation rate N over its zero-gradient value.

    Dey & Narasimha (1989), as summarized in Narasimha (1990, NASA CR-187487,
    Fig. 3), against the Thwaites parameter L = th^2/nu dU/dx: 0.25 L^2 for
    0 <= L <= 0.12; for adverse gradients a fit to the curve through
    Gostelow's (1989) data, 5.25 L^2 + 25.2 L^3 (-0.07 <= L <= 0).
    """
    if np.real(L) >= 0:
        return 0.25*L**2, 0.5*L
    return 5.25*L**2 + 25.2*L**3, 10.5*L + 75.6*L**2


def transition_rate(param, U1, U2, x):
    """C such that n crosses the closure band in the transition length.

    Constant (param.capt_L momentum thicknesses) unless
    param.capt_Lmodel == 'narasimha': then the Dhawan-Narasimha 10-90%
    length 1.863 lambda, lambda / th = sqrt(0.41 Re_th / N) (Dey et al.
    1995, eqs. 3, 10), with N = N0 + dN(L) from local Re_th, th and dU/dx;
    N0 = param.capt_N0 (default 0.7e-3, Narasimha's value for
    turbulence-driven transition). Returns C, dC/dU (10), dC/dx (2).
    """
    cfg = _cfg(param)
    if getattr(param, 'capt_Lmodel', 'constant') != 'narasimha':
        return cfg['C'], np.zeros(10), np.zeros(2)
    wc = cfg['cl'][1]
    pt = copy.copy(param); pt.turb = False; pt.wake = False
    R1, R1_U = mf.get_Ret(U1[TURB_ROWS], pt); R2, R2_U = mf.get_Ret(U2[TURB_ROWS], pt)
    Rt = 0.5*(R1 + R2)
    Rt_U = 0.5*np.concatenate((_e5(R1_U, TURB_ROWS), _e5(R2_U, TURB_ROWS)))
    E = np.eye(5)
    th = 0.5*(U1[TH] + U2[TH]); th_U = 0.5*np.concatenate((E[TH], E[TH]))
    ue = 0.5*(U1[UE] + U2[UE]); ue_U = 0.5*np.concatenate((E[UE], E[UE]))
    dx = x[1] - x[0]
    du = (U2[UE] - U1[UE])/dx
    du_U = np.concatenate((-E[UE], E[UE]))/dx
    du_x = -du/dx*np.r_[-1., 1.]
    # Thwaites parameter th^2/nu dU/dx = Re_th (th/ue) dU/dx, smoothly bounded
    Lt = Rt*th/ue*du
    Lt_U = (Rt_U*th/ue + Rt*th_U/ue - Rt*th/ue**2*ue_U)*du + Rt*th/ue*du_U
    Lt_x = Rt*th/ue*du_x
    a = 0.1
    t = np.tanh(Lt/a); Lc = a*t; Lc_L = 1 - t**2
    dN, dN_L = _dN(Lc)
    N = getattr(param, 'capt_N0', 0.7e-3) + dN
    N_U = dN_L*Lc_L*Lt_U; N_x = dN_L*Lc_L*Lt_x
    Lth = 1.863*np.sqrt(0.41*Rt/N)
    Lth_U = 0.5*Lth*(Rt_U/Rt - N_U/N); Lth_x = -0.5*Lth*N_x/N
    C = 2.2*wc/Lth
    return C, -C/Lth*Lth_U, -C/Lth*Lth_x


def _e5(v4, rows):
    """Embed a 4-state linearization in the 5-state layout."""
    out = np.zeros(5, dtype=np.result_type(v4, float))
    out[rows] += v4
    return out


def _e10(v8, rows):
    """Embed an 8-entry (two-station) linearization in 10 entries."""
    out = np.zeros(10, dtype=np.result_type(v8, float))
    out[np.array(rows)] += v8[0:4]
    out[5 + np.array(rows)] += v8[4:8]
    return out


def _up(upw, upw_U, f1, f1_U, f2, f2_U):
    """mfoil.upwind for 5-state linearizations (upw_U has 10 entries)."""
    return ((1 - upw)*f1 + upw*f2,
            -upw_U*f1 + upw_U*f2 + np.concatenate(((1 - upw)*f1_U, upw*f2_U)))


def _blended(name, U5, param):
    """Hs, cf or cDi: (1 - g) laminar + g turbulent at one station."""
    fn = {'Hs': mf.get_Hs, 'cf': mf.get_cf, 'cDi': mf.get_cDi}[name]
    pt = copy.copy(param); pt.turb = True
    ft, ft_U = fn(U5[TURB_ROWS], pt)
    if param.wake:
        return ft, _e5(ft_U, TURB_ROWS)
    pl = copy.copy(param); pl.turb = False
    fl, fl_U = fn(U5[LAM_ROWS], pl)
    g, g_n = gamma(U5[AMP], param)
    f_U = (1 - g)*_e5(fl_U, LAM_ROWS) + g*_e5(ft_U, TURB_ROWS)
    f_U[AMP] += (ft - fl)*g_n
    return (1 - g)*fl + g*ft, f_U


def _xt(name, U5, x, param):
    """cf x / th or cDi x / th with the blended closure."""
    f, f_U = _blended(name, U5, param)
    v = f*x/U5[TH]
    v_U = f_U*x/U5[TH]; v_U[TH] -= v/U5[TH]
    return v, v_U, f/U5[TH]


def residual_captured(param, x, Uin, Aux, trip=None):
    """Four residuals of one interval; U is 5x2. Returns R, R_U (4x10), R_x.

    trip: fraction of the interval at which a trip sits, or None.

    Same discretization as mfoil.residual_station, with the closures Hs, cf
    and cDi blended per station by g(n), one shear-lag equation whose
    equilibrium target is released from cttr to cteq by r(n), and one
    amplification equation (laminar rate plus the post-onset growth).
    """
    U = Uin.copy()
    U[DS, :] -= Aux
    U1, U2 = U[:, 0], U[:, 1]
    Um = 0.5*(U1 + U2)
    Ut1, Ut2 = U1[TURB_ROWS], U2[TURB_ROWS]
    pt = copy.copy(param); pt.turb = True
    pl = copy.copy(param); pl.turb = False; pl.wake = False
    E = np.eye(5)
    Z10 = np.zeros(10)

    def two(a, b):
        return np.concatenate((a, b))

    th, ct, n = U[TH, :], U[CT, :], U[AMP, :]
    uk1, uk1_u = mf.get_uk(U1[UE], param)
    uk2, uk2_u = mf.get_uk(U2[UE], param)
    thlog = np.log(th[1]/th[0])
    thlog_U = two(-E[TH]/th[0], E[TH]/th[1])
    uelog = np.log(uk2/uk1)
    uelog_U = two(-uk1_u/uk1*E[UE], uk2_u/uk2*E[UE])
    xlog = np.log(x[1]/x[0]); xlog_x = np.r_[-1./x[0], 1./x[1]]
    dx = x[1] - x[0]; dx_x = np.r_[-1., 1.]

    upw, upw8 = mf.get_upw(Ut1, Ut2, pt)
    upw_U = _e10(upw8, TURB_ROWS)

    H1, H1_U = mf.get_H(Ut1); H2, H2_U = mf.get_H(Ut2)
    H = 0.5*(H1 + H2); H_U = 0.5*two(_e5(H1_U, TURB_ROWS), _e5(H2_U, TURB_ROWS))

    Hs1, Hs1_U = _blended('Hs', U1, param)
    Hs2, Hs2_U = _blended('Hs', U2, param)
    Hs, Hs_U = _up(0.5, Z10, Hs1, Hs1_U, Hs2, Hs2_U)
    Hslog = np.log(Hs2/Hs1)
    Hslog_U = two(-Hs1_U/Hs1, Hs2_U/Hs2)

    if param.simi:
        thlog = 0.; thlog_U = 0*thlog_U
        Hslog = 0.; Hslog_U = 0*Hslog_U
        uelog = 1.; uelog_U = 0*uelog_U
        xlog = 1.; xlog_x = np.r_[0., 0.]
        dx = 0.5*(x[0] + x[1]); dx_x = np.r_[0.5, 0.5]

    Hw1, Hw1_U = mf.get_Hw(Ut1, Aux[0]); Hw2, Hw2_U = mf.get_Hw(Ut2, Aux[1])
    Hw = 0.5*(Hw1 + Hw2)
    Hw_U = 0.5*two(_e5(Hw1_U, TURB_ROWS), _e5(Hw2_U, TURB_ROWS))

    # ---- shear-lag equation (turbulent closures; target released by r(n))
    salog = np.log(ct[1]/ct[0])
    salog_U = two(-E[CT]/ct[0], E[CT]/ct[1])
    de1, de1_U = mf.get_de(Ut1, pt); de2, de2_U = mf.get_de(Ut2, pt)
    de, de_U = _up(0.5, Z10, de1, _e5(de1_U, TURB_ROWS), de2, _e5(de2_U, TURB_ROWS))
    Us1, Us1_U = mf.get_Us(Ut1, pt); Us2, Us2_U = mf.get_Us(Ut2, pt)
    Us, Us_U = _up(0.5, Z10, Us1, _e5(Us1_U, TURB_ROWS), Us2, _e5(Us2_U, TURB_ROWS))
    Hk1, Hk1_U = mf.get_Hk(Ut1, pt); Hk2, Hk2_U = mf.get_Hk(Ut2, pt)
    Hk, Hk_U = _up(upw, upw_U, Hk1, _e5(Hk1_U, TURB_ROWS), Hk2, _e5(Hk2_U, TURB_ROWS))
    Rt1, Rt1_U = mf.get_Ret(Ut1, pt); Rt2, Rt2_U = mf.get_Ret(Ut2, pt)
    Ret, Ret_U = _up(0.5, Z10, Rt1, _e5(Rt1_U, TURB_ROWS), Rt2, _e5(Rt2_U, TURB_ROWS))
    cft1, cft1_U = mf.get_cf(Ut1, pt); cft2, cft2_U = mf.get_cf(Ut2, pt)
    cft, cft_U = _up(upw, upw_U, cft1, _e5(cft1_U, TURB_ROWS), cft2, _e5(cft2_U, TURB_ROWS))
    dsa = 0.5*(U[DS, 0] + U[DS, 1]); dsa_U = 0.5*two(E[DS], E[DS])
    uq, uq_U = mf.get_uq(dsa, dsa_U, cft, cft_U, Hk, Hk_U, Ret, Ret_U, pt)
    tg, tg_U = [], []
    for Uk in (U1, U2):
        ceq, ceq_U = mf.get_cteq(Uk[TURB_ROWS], copy.copy(pt))
        if param.wake:
            tg.append(ceq); tg_U.append(_e5(ceq_U, TURB_ROWS))
            continue
        ctr, ctr_U = mf.get_cttr(Uk[TURB_ROWS], copy.copy(pt))
        r, r_n = release(Uk[AMP], param)
        t_U = (1 - r)*_e5(ctr_U, TURB_ROWS) + r*_e5(ceq_U, TURB_ROWS)
        t_U[AMP] += (ceq - ctr)*r_n
        tg.append((1 - r)*ctr + r*ceq); tg_U.append(t_U)
    tgt, tgt_U = _up(upw, upw_U, tg[0], tg_U[0], tg[1], tg_U[1])
    saa, saa_U = _up(upw, upw_U, ct[0], E[CT], ct[1], E[CT])
    if param.wake:
        rb, rb_U = 1.0, Z10
    else:
        ra, ra_n = release(n[0], param); rc, rc_n = release(n[1], param)
        rb = 0.5*(ra + rc); rb_U = 0.5*two(ra_n*E[AMP], rc_n*E[AMP])
    Clag = param.SlagK/param.GB/(1. + Us)
    Clag_U = -Clag/(1. + Us)*Us_U
    ald = param.Dlr if param.wake else 1.0
    pg = 2*de*(uq*dx - uelog)*param.Cuq
    pg_U = 2*de_U*(uq*dx - uelog)*param.Cuq + 2*de*(uq_U*dx - uelog_U)*param.Cuq
    Rlag = Clag*(tgt - ald*saa)*dx - 2*de*salog + rb*pg
    Rlag_U = (Clag_U*(tgt - ald*saa)*dx + Clag*(tgt_U - ald*saa_U)*dx
              - 2*de_U*salog - 2*de*salog_U + rb_U*pg + rb*pg_U)
    Rlag_x = Clag*(tgt - ald*saa)*dx_x + rb*2*de*uq*dx_x*param.Cuq
    if param.simi:
        # the ODE degenerates at the stagnation point (dx ~ 1e-6); pin ct
        # algebraically, with the same 2 de scale as the transport term
        Rlag = 2*de*np.log(ct[1]/tg[1])
        Rlag_U = (2*de_U*np.log(ct[1]/tg[1])
                  + 2*de*np.concatenate((np.zeros(5), E[CT]/ct[1] - tg_U[1]/tg[1])))
        Rlag_x = np.zeros(2)

    # ---- amplification equation
    if param.wake:
        Ramp = n[1] - n[0]; Ramp_U = two(-E[AMP], E[AMP]); Ramp_x = np.zeros(2)
    elif param.simi:
        Ramp = n[0] + n[1]; Ramp_U = two(E[AMP], E[AMP]); Ramp_x = np.zeros(2)
    else:
        d1, d1_U = mf.get_damp(U1[LAM_ROWS], copy.copy(pl))
        d2, d2_U = mf.get_damp(U2[LAM_ROWS], copy.copy(pl))
        damp = 0.5*(d1 + d2)
        damp_U = 0.5*two(_e5(d1_U, LAM_ROWS), _e5(d2_U, LAM_ROWS))
        # rate (1 - gs) damp + gs C / th: profile independent after onset
        C, C_U, C_x = transition_rate(param, U1, U2, x)
        s1, s1_n = growth_switch(n[0], param); s2, s2_n = growth_switch(n[1], param)
        q1, q1_U, q1_C = _rate(d1, _e5(d1_U, LAM_ROWS), s1, s1_n, C, th[0], 0, E)
        q2, q2_U, q2_C = _rate(d2, _e5(d2_U, LAM_ROWS), s2, s2_n, C, th[1], 1, E)
        damp = 0.5*(q1 + q2)
        damp_U = 0.5*two(q1_U, q2_U) + 0.5*(q1_C + q2_C)*C_U
        Ramp = n[1] - n[0] - damp*dx
        Ramp_U = two(-E[AMP], E[AMP]) - damp_U*dx
        Ramp_x = -damp*dx_x - 0.5*(q1_C + q2_C)*C_x*dx
        if trip is not None:
            # forced transition (constant C only): see the module docstring
            f = trip; g = (1 - f)*dx
            npre = n[0] + q1*f*dx
            J, sig = _trip_jump(npre, param)
            a = 1 - sig
            npost = npre + J
            ss, ss_n = growth_switch(npost, param)
            d1_5 = _e5(d1_U, LAM_ROWS)
            qs, qs_U, _ = _rate(d1, d1_5, ss, 0.0, C, th[0], 0, E)
            qs_n = (qs - d1)/ss*ss_n if np.real(ss) > 1e-300 else 0.0
            k = 1 + 0.5*g*qs_n
            npost_U = a*(E[AMP] + q1_U*f*dx)
            Ramp = n[1] - npost - 0.5*(qs + q2)*g
            Ramp_U = two(-k*npost_U - 0.5*g*qs_U, E[AMP] - 0.5*g*q2_U)
            Ramp_x = -(k*a*q1*f + 0.5*(qs + q2)*(1 - f))*dx_x

    # ---- momentum equation (blended cf)
    Ms1, Ms1_U = mf.get_Mach2(Ut1, param); Ms2, Ms2_U = mf.get_Mach2(Ut2, param)
    Ms = 0.5*(Ms1 + Ms2); Ms_U = 0.5*two(_e5(Ms1_U, TURB_ROWS), _e5(Ms2_U, TURB_ROWS))
    c1, c1_U, c1_x = _xt('cf', U1, x[0], param)
    c2, c2_U, c2_x = _xt('cf', U2, x[1], param)
    cm, cm_U, cm_x = _xt('cf', Um, 0.5*(x[0] + x[1]), param)
    cfxt = 0.25*c1 + 0.5*cm + 0.25*c2
    cfxt_U = 0.25*two(c1_U + cm_U, cm_U + c2_U)
    cfxt_x = 0.25*np.r_[c1_x + cm_x, cm_x + c2_x]
    Rmom = thlog + (2 + H + Hw - Ms)*uelog - 0.5*xlog*cfxt
    Rmom_U = thlog_U + (H_U + Hw_U - Ms_U)*uelog + (2 + H + Hw - Ms)*uelog_U - 0.5*xlog*cfxt_U
    Rmom_x = -0.5*xlog_x*cfxt - 0.5*xlog*cfxt_x

    # ---- shape equation (blended Hs, cf, cDi)
    D1, D1_U, D1_x = _xt('cDi', U1, x[0], param)
    D2, D2_U, D2_x = _xt('cDi', U2, x[1], param)
    cDixt, cDixt_U = _up(upw, upw_U, D1, D1_U, D2, D2_U)
    cDixt_x = np.r_[(1 - upw)*D1_x, upw*D2_x]
    cfxtu, cfxtu_U = _up(upw, upw_U, c1, c1_U, c2, c2_U)
    cfxtu_x = np.r_[(1 - upw)*c1_x, upw*c2_x]
    Hss1, Hss1_U = mf.get_Hss(Ut1, param); Hss2, Hss2_U = mf.get_Hss(Ut2, param)
    Hss, Hss_U = _up(0.5, Z10, Hss1, _e5(Hss1_U, TURB_ROWS), Hss2, _e5(Hss2_U, TURB_ROWS))
    Rshape = Hslog + (2*Hss/Hs + 1 - H - Hw)*uelog + xlog*(0.5*cfxtu - cDixt)
    Rshape_U = (Hslog_U + (2*Hss_U/Hs - 2*Hss/Hs**2*Hs_U - H_U - Hw_U)*uelog
                + (2*Hss/Hs + 1 - H - Hw)*uelog_U + xlog*(0.5*cfxtu_U - cDixt_U))
    Rshape_x = xlog_x*(0.5*cfxtu - cDixt) + xlog*(0.5*cfxtu_x - cDixt_x)

    R = np.array([Rmom, Rshape, Ramp, Rlag])
    R_U = np.vstack((Rmom_U, Rshape_U, Ramp_U, Rlag_U))
    R_x = np.vstack((Rmom_x, Rshape_x, Ramp_x, Rlag_x))
    return R, R_U, R_x


def stagnation_state(U, x):
    """5-state version of mfoil.stagnation_state (ue row quadratic)."""
    U1, U2, x1, x2 = U[:, 0], U[:, 1], x[0], x[1]
    dx = x2 - x1; dx_x = np.r_[-1, 1]
    rx = x2/x1; rx_x = np.r_[-rx, 1]/x1
    w1 = x2/dx; w1_x = -w1/dx*dx_x + np.r_[0, 1]/dx
    w2 = -x1/dx; w2_x = -w2/dx*dx_x + np.r_[-1, 0]/dx
    Ust = U1*w1 + U2*w2
    wk1 = rx/dx; wk1_x = rx_x/dx - wk1/dx*dx_x
    wk2 = -1/(rx*dx); wk2_x = -wk2*(rx_x/rx + dx_x/dx)
    K = wk1*U1[UE] + wk2*U2[UE]
    K_x = U1[UE]*wk1_x + U2[UE]*wk2_x
    xst = 1e-6
    Ust[UE] = K*xst
    E = np.eye(5); E[UE, UE] = 0
    Ust_U = np.hstack((w1*E, w2*E))
    Ust_U[UE, UE] = wk1*xst; Ust_U[UE, 5 + UE] = wk2*xst
    Ust_x = np.outer(U1, w1_x) + np.outer(U2, w2_x)
    Ust_x[UE] = K_x*xst
    # ct: take the first node's value. Extrapolating can turn it negative
    # (cteq changes fast at Re_theta ~ 10), and ct is not used there anyway.
    Ust[CT] = U1[CT]
    Ust_U[CT] = 0; Ust_U[CT, CT] = 1
    Ust_x[CT] = 0
    return Ust, Ust_U, Ust_x, xst


def wake_sys(M, param):
    """Wake initialization: th, ds sums; theta-weighted ct; n = 0."""
    il = M.vsol.Is[0][-1]; Ul = M.glob.U[:, il]
    iu = M.vsol.Is[1][-1]; Uu = M.glob.U[:, iu]
    iw = M.vsol.Is[2][0]; Uw = M.glob.U[:, iw]
    hTE = mf.TE_info(M.foil.x)[1]
    thsum = Ul[TH] + Uu[TH]
    ctw = (Ul[CT]*Ul[TH] + Uu[CT]*Uu[TH])/thsum
    R = np.r_[Uw[TH] - thsum, Uw[DS] - (Ul[DS] + Uu[DS] + hTE),
              Uw[AMP], Uw[CT] - ctw]
    R_U = np.zeros((4, 15), dtype=np.result_type(R, ctw))
    for b in (0, 5):
        R_U[0, b + TH] = -1; R_U[1, b + DS] = -1
    R_U[0, 10 + TH] = 1; R_U[1, 10 + DS] = 1; R_U[2, 10 + AMP] = 1
    R_U[3, 10 + CT] = 1
    R_U[3, 0 + CT] = -Ul[TH]/thsum; R_U[3, 0 + TH] = -(Ul[CT] - ctw)/thsum
    R_U[3, 5 + CT] = -Uu[TH]/thsum; R_U[3, 5 + TH] = -(Uu[CT] - ctw)/thsum
    return R, R_U, [il, iu, iw]


def build_glob_sys(M):
    """Assemble R (4 Nsys), R_U (4 Nsys x 5 Nsys) and R_x (4 Nsys x Nsys)."""
    Nsys = M.glob.Nsys
    dt = M.glob.U.dtype
    R = np.zeros(4*Nsys, dtype=dt)
    R_U = sparse.lil_matrix((4*Nsys, 5*Nsys), dtype=dt)
    R_x = sparse.lil_matrix((4*Nsys, Nsys), dtype=dt)

    def put(i, Ri, Ri_U, J, Ri_x=None):
        Ig = slice(4*i, 4*i + 4)
        R[Ig] += Ri
        for j, nj in enumerate(J):
            R_U[Ig, 5*nj:5*nj + 5] += Ri_U[:, 5*j:5*j + 5]
            if Ri_x is not None:
                R_x[Ig, nj] += Ri_x[:, j:j + 1]

    for si in range(3):
        Is = M.vsol.Is[si]
        xi = M.isol.xi[Is]
        N = len(Is)
        U = M.glob.U[:, Is]
        Aux = np.zeros(N, dtype=dt)
        param = mf.build_param(M, si)
        if si == 2:
            Aux[:] = M.vsol.wgap
        i0 = 1 if (si < 2) and (np.real(xi[0]) < 1e-8*np.real(xi[-1])) else 0
        if si < 2:
            Ip = [i0, i0 + 1]
            Ust, Ust_U, Ust_x, xst = stagnation_state(U[:, Ip], xi[Ip])
            param.turb, param.simi = False, True
            R1, R1_Ut, R1_xs = residual_captured(param, np.r_[xst, xst],
                                                 np.stack((Ust, Ust), axis=-1),
                                                 Aux[[i0, i0]])
            param.simi = False
            R1_Ust = R1_Ut[:, 0:5] + R1_Ut[:, 5:10]
            J = [Is[i0], Is[i0 + 1]]
            if i0 == 1:  # first node sits on the stagnation point
                rows = [TH, DS, CT, AMP]
                Rs = U[rows, 0] - Ust[rows]
                Rs_U = np.zeros((4, 15), dtype=dt)
                Rs_U[:, 0:5] = np.eye(5)[rows]
                Rs_U[:, 5:15] = -Ust_U[rows]
                put(Is[0], Rs, Rs_U, [Is[0]] + J)
                R_x[4*Is[0]:4*Is[0] + 4, J] = -Ust_x[rows]
            put(Is[i0], R1, R1_Ust @ Ust_U, J, R1_Ust @ Ust_x)
        else:
            R1, R1_U, J = wake_sys(M, param)
            put(Is[i0], R1, R1_U, J)
            param.turb, param.wake = True, True
        itrip, ftrip = _find_trip(M, si)
        for i in range(i0 + 1, N):
            Ip = [i - 1, i]
            Ri, Ri_U, Ri_x = residual_captured(
                param, xi[Ip], U[:, Ip], Aux[Ip],
                ftrip if i == itrip else None)
            put(Is[i], Ri, Ri_U, [Is[i - 1], Is[i]], Ri_x)

    # stagnation-point movement: R_ue += R_x x_st st_ue (as in mfoil)
    x_st = np.concatenate((-M.isol.sgnue, -np.ones(M.wake.N)))
    R_x = R_x.tocsr()
    R_st = R_x @ x_st
    R_U = R_U.tocsc()
    Ist, st_ue = M.isol.Istag, M.isol.sstag_ue
    col = sparse.lil_matrix((4*Nsys, 5*Nsys), dtype=np.result_type(R_st, st_ue))
    col[:, 5*Ist[0] + UE] = (R_st*st_ue[0])[:, None]
    col[:, 5*Ist[1] + UE] = (R_st*st_ue[1])[:, None]
    M.glob.R, M.glob.R_U, M.glob.R_x = R, (R_U + col.tocsc()).tocsr(), R_x


def solve_glob(M):
    """Newton update for [th, ds, ct, ue, n] with the ue coupling rows."""
    assert not M.oper.givencl, 'captured prototype: alpha must be prescribed'
    Nsys = M.glob.Nsys
    U = M.glob.U
    ue = U[UE, :].copy(); ds = U[DS, :]
    uemax = np.max(np.abs(np.real(ue)))
    ue = np.where(np.real(ue) > 1e-10*uemax, ue, 1e-10*uemax)
    ueinv = mf.get_ueinv(M)
    R = np.concatenate((M.glob.R, ue - (ueinv + M.vsol.ue_m @ (ds*ue))))
    Ids = np.arange(Nsys)*5 + DS
    Iue = np.arange(Nsys)*5 + UE
    Bue = (sparse.identity(Nsys) - sparse.csr_matrix(M.vsol.ue_m @ np.diag(ds)))
    Bds = sparse.csr_matrix(-M.vsol.ue_m @ np.diag(ue))
    Sue = sparse.csr_matrix((np.ones(Nsys), (Iue, np.arange(Nsys))), shape=(5*Nsys, Nsys))
    Sds = sparse.csr_matrix((np.ones(Nsys), (Ids, np.arange(Nsys))), shape=(5*Nsys, Nsys))
    coupling = Bue @ Sue.T + Bds @ Sds.T
    R_V = sparse.vstack((M.glob.R_U, coupling)).tocsc()
    dV = -splinalg.spsolve(R_V, R)
    M.glob.dU = np.reshape(dV, (5, Nsys), order='F')


def update_state(M):
    """Damped update keeping th, ds, ct positive and changes bounded."""
    U, dU = M.glob.U, M.glob.dU
    omega = 1.0
    for k in (TH, DS):
        fmin = np.min(np.real(dU[k]/U[k]))
        if fmin < -0.5:
            omega = min(omega, abs(0.5/fmin))
    # limit ct decreases only where ct matters (released or blended);
    # upstream it just tracks cttr and is kept positive after the update
    nr = np.real(U[AMP])
    act = np.maximum(np.real(release(nr, M.param)[0]),
                     np.real(gamma(nr, M.param)[0])) > 1e-3
    for si in (2,):
        act[M.vsol.Is[si]] = True
    r = np.real(dU[CT]/U[CT])
    if np.any(act) and np.min(r[act]) < -0.8:
        omega = min(omega, abs(0.8/np.min(r[act])))
    # limit n changes only where n still matters (closure blend unsaturated)
    live = np.real(gamma(U[AMP], M.param)[0]) < 0.999
    if np.any(live):
        dn = np.max(np.abs(np.real(dU[AMP, live])))
        if dn > 2:
            omega = min(omega, 2/dn)
    fmax = np.max(np.abs(np.real(dU[UE])))/M.oper.Vinf
    if fmax > 0.2:
        omega = min(omega, 0.2/fmax)
    mf.vprint(M.param, 2, '  captured update: omega = %.5f' % omega)
    M.glob.U = U + omega*dU
    low = np.real(M.glob.U[CT]) < 0.2*np.real(U[CT])
    M.glob.U[CT, low] = 0.2*U[CT, low]
    for si in range(3):  # fix bad Hk, as mfoil does
        Hkmin = 1.00005 if si == 2 else 1.02
        param = mf.build_param(M, si)
        for j in M.vsol.Is[si]:
            Hk = mf.get_Hk(M.glob.U[TURB_ROWS, j], param)[0]
            if np.real(Hk) < Hkmin:
                M.glob.U[DS, j] += 2*(Hkmin - Hk)*M.glob.U[DS, j]


def convert_state(M):
    """Standard [th, ds, sa, ue] solution -> [th, ds, ct, ue, n].

    th, ds, ue are kept. n is marched from the stagnation point with the
    captured amplification equation (laminar rate + gs C / th) on the given
    th, ds, ue, so it is consistent with the residual. ct is cttr upstream
    of the cut-cell transition and the standard ct state downstream.
    """
    U4 = M.glob.U
    Nsys = U4.shape[1]
    U = np.zeros((5, Nsys), dtype=U4.dtype)
    U[0:4] = U4
    C = _C(M.param)
    for si in range(3):
        Is = M.vsol.Is[si]
        param = mf.build_param(M, si)
        itrip, ftrip = _find_trip(M, si)
        for k, j in enumerate(Is):
            if si == 2 or M.vsol.turb[j]:
                U[CT, j] = U4[2, j]
            else:
                p = copy.copy(param); p.turb = True
                U[CT, j] = mf.get_cttr(U4[:, j], p)[0]
            if si == 2 or k == 0:
                U[AMP, j] = 0.0
                continue
            p = copy.copy(param); p.turb = False
            jm = Is[k - 1]
            Um = U[LAM_ROWS, jm]; Uj = U[LAM_ROWS, j].copy()
            dx = M.isol.xi[j] - M.isol.xi[jm]
            Uj[2] = Um[2]
            U5m = U[:, jm]; U5j = U[:, j].copy()
            C = transition_rate(param, U5m, U5j, M.isol.xi[[jm, j]])[0]
            sm = growth_switch(Um[2], M.param)[0]
            E = np.eye(5)
            d1 = mf.get_damp(Um, p)[0]
            q1 = _rate(d1, np.zeros(5), sm, 0.0, C, Um[0], 0, E)[0]
            for _ in range(60):
                d2, d2_U = mf.get_damp(Uj, p)
                s2, s2_n = growth_switch(Uj[2], M.param)
                q2, q2_U, _ = _rate(d2, _e5(d2_U, LAM_ROWS), s2, s2_n, C, Uj[0], 1, E)
                q2_n = q2_U[AMP]
                Ra = Uj[2] - Um[2] - 0.5*(q1 + q2)*dx
                Ra_n = 1 - 0.5*q2_n*dx
                if k == itrip:  # as residual_captured
                    npre = Um[2] + q1*ftrip*dx
                    npost = npre + _trip_jump(npre, M.param)[0]
                    qs = _rate(d1, np.zeros(5), growth_switch(npost, M.param)[0],
                               0.0, C, Um[0], 0, E)[0]
                    Ra = Uj[2] - npost - 0.5*(qs + q2)*(1 - ftrip)*dx
                    Ra_n = 1 - 0.5*q2_n*(1 - ftrip)*dx
                Uj[2] -= Ra/Ra_n
                if abs(Ra) < 1e-12:
                    break
            U[AMP, j] = Uj[2]
    M.glob.U = U


def solve_captured(M, nNewton=None):
    """Standard viscous solve, then converge the captured system."""
    M.solve()
    convert_state(M)
    M.glob.realloc = True
    M.glob.conv = False
    nNewton = nNewton or M.param.niglob
    for it in range(nNewton):
        build_glob_sys(M)
        Rn = np.sqrt(np.sum(np.abs(np.real(M.glob.R))**2))
        mf.vprint(M.param, 1, 'captured Newton %d, Rnorm = %.10e' % (it, Rn))
        if Rn < max(M.param.rtol, 1e-10):
            M.glob.conv = True
            break
        solve_glob(M)
        update_state(M)
        mf.stagpoint_move(M)
    # report transition where gamma crosses 1/2; forces from the wake
    g = gamma(np.real(M.glob.U[AMP]), M.param)[0]
    M.vsol.turb = g > 0.5
    for j in M.vsol.Is[2]:
        M.vsol.turb[j] = True
    U5 = M.glob.U
    M.glob.U = U5[0:4]
    mf.calc_force(M)
    M.glob.U = U5
    return M.glob.conv
