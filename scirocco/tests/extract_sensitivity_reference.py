#!/usr/bin/env python3
"""
Generate machine-precision sensitivity reference data using complex-step
differentiation on the mfoil solver.

Outputs: sensitivity_reference.h5

HDF5 structure:
  viscous_shape/     NACA 0012, 100 panels, alpha=3, Re=1e6
  viscous_alpha/     same base case
  viscous_mach/      NACA 0012, 100 panels, alpha=3, Re=1e6, Ma=0.3
  viscous_re/        NACA 0012, 100 panels, alpha=3, Re=1e6
  viscous_shape_onset/  NACA 2412, 160 panels, alpha=4, Re=1e6, upper-surface
                     nodes 85-93 around the amplification onset (issue #22)
  inviscid_shape/    NACA 2412, 80 panels, alpha=5
  inviscid_alpha/    NACA 2412, 80 panels, alpha=5
  inviscid_mach/     NACA 0012, 100 panels, alpha=3, Ma=0.3
  envelope_shape/    NACA 2412, 160 panels, alpha=4, Re=1e6, onset='envelope'
  envelope_scalar/   same case: d/dalpha and d/dRe
  envelope_mach/     same case at Ma=0.3: d/dMa
  viscous_re_laminar_te/  NACA 2412, 160 panels, alpha=4, Re=1e6, default
                     onset: d/dRe with the lower surface laminar at the TE
  captured_primal/   NACA 2412, 160 panels, cumulative arc length:
                     TransitionModel::Captured cl/cd/cm (mfoil_captured)
  captured_shape/    same, alpha=4, Re=1e6: CS dcl/dcd/dcm/dx around transition
  captured_scalar/   same: CS d/dalpha, d/dRe (Ma=0); d/dMa at Ma=0.3
  trip_shape/        NACA 2412, 160 panels, alpha=4, Re=1e6, forced transition
                     at x/c = 0.2 (upper) and 0.6 (lower): CS dcl/dcd/dcm/dx
                     at the nodes around both trips
  trip_scalar/       same: CS d/dalpha, d/dRe (Ma=0); d/dMa at Ma=0.3
  cf_floor_shape/    same case with the laminar skin-friction floor
                     (Param::cf_laminar_floor), lower trip at x/c = 0.01 and
                     free transition above: CS dcl/dcd/dcm/dx around both
  cf_floor_scalar/   same: CS d/dalpha, d/dRe (Ma=0); d/dMa at Ma=0.3

Usage:
  cd mfoil
  python ../scirocco/tests/extract_sensitivity_reference.py
  python ../scirocco/tests/extract_sensitivity_reference.py --onset-only
  python ../scirocco/tests/extract_sensitivity_reference.py --envelope-only

--onset-only (re)writes just viscous_shape_onset/ and --envelope-only just the
envelope_* groups into the existing file; the two can be combined.
--captured-only (re)writes just captured_primal/.
--laminar-te-only (re)writes just viscous_re_laminar_te/.
--captured-cs-only (re)writes just captured_shape/ and captured_scalar/.

--captured-trip-only (re)writes just captured_trip_shape/ and
captured_trip_scalar/: the captured model with the trips of trip_shape/.
--trip-only (re)writes just trip_shape/ and trip_scalar/.
--cf-floor-only (re)writes just cf_floor_shape/ and cf_floor_scalar/.
"""

import sys
import hashlib
import os
import time
import numpy as np
import h5py

# Add mfoil/ to path for mfoil_cs
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'mfoil'))
import mfoil_cs  # noqa: E402 — patches mfoil for complex-step compatibility
from mfoil_cs import complex_step_mode  # noqa: E402
import mfoil     # noqa: E402

CS_STEP = 1e-20

OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'sensitivity_reference.h5')


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def solve_inviscid(naca, npanel, alpha, Ma=0.0):
    """Solve inviscid case and return mfoil object."""
    M = mfoil.mfoil(naca=naca, npanel=npanel)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=alpha, Ma=Ma)
    M.solve()
    return M


def solve_viscous(naca, npanel, alpha, Re, Ma=0.0):
    """Solve viscous case and return mfoil object."""
    M = mfoil.mfoil(naca=naca, npanel=npanel)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=alpha, Re=Re, Ma=Ma)
    M.solve()
    return M


def recompute_geom(M):
    """Recompute arc-length distribution from (possibly complex) foil.x.

    Does NOT recompute chord — chord is a geometric property set during
    paneling and should not change when perturbing node positions.
    """
    x = M.foil.x
    N = M.foil.N
    ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
    s = np.zeros(N, dtype=x.dtype)
    s[1:] = np.cumsum(ds)
    M.foil.s = s


def cs_inviscid_scalar(naca, npanel, alpha, Ma, var, step):
    """
    Complex-step derivative of (cl, cm) w.r.t. a scalar operating condition.
    var: 'alpha' or 'Ma'
    """
    kwargs = dict(alpha=alpha, Ma=Ma)
    if var == 'alpha':
        kwargs['alpha'] = alpha + 1j * step
    elif var == 'Ma':
        kwargs['Ma'] = Ma + 1j * step
    else:
        raise ValueError(f"Unknown var: {var}")

    # Create mfoil object OUTSIDE complex mode (scipy spline needs real arrays)
    M = mfoil.mfoil(naca=naca, npanel=npanel)
    M.param.verb = 0
    M.param.doplot = False
    # Now enter complex mode for setoper + solve
    with complex_step_mode():
        M.setoper(**kwargs)
        M.solve()

    dcl = np.imag(M.post.cl) / step
    dcm = np.imag(M.post.cm) / step
    return dcl, dcm


def cs_viscous_scalar(naca, npanel, alpha, Re, Ma, var, step, onset='cubic',
                      xtrip=None, cf_floor=False):
    """
    Complex-step derivative of (cl, cd, cm) w.r.t. a scalar operating condition.
    var: 'alpha', 'Ma', or 'Re'
    """
    kwargs = dict(alpha=alpha, Re=Re, Ma=Ma)
    if var == 'alpha':
        kwargs['alpha'] = alpha + 1j * step
    elif var == 'Ma':
        kwargs['Ma'] = Ma + 1j * step
    elif var == 'Re':
        kwargs['Re'] = Re + 1j * step
    else:
        raise ValueError(f"Unknown var: {var}")

    # Create mfoil object OUTSIDE complex mode (scipy spline needs real arrays)
    M = mfoil.mfoil(naca=naca, npanel=npanel)
    M.param.verb = 0
    M.param.doplot = False
    M.param.onset = onset
    M.param.xtrip = xtrip
    M.param.cf_laminar_floor = cf_floor
    # Now enter complex mode for setoper + solve
    with complex_step_mode():
        M.setoper(**kwargs)
        M.solve()

    dcl = np.imag(M.post.cl) / step
    dcd = np.imag(M.post.cd) / step
    dcm = np.imag(M.post.cm) / step
    return dcl, dcd, dcm


def cs_inviscid_shape(naca, npanel, alpha, Ma, test_nodes):
    """
    Complex-step shape sensitivities for inviscid case.
    Returns dcl_dx [2, N+1], dcm_dx [2, N+1] (nonzero only at test_nodes).
    """
    # M_base must use cumulative foil.s to match perturbation calls (see
    # cs_viscous_shape for rationale).
    M_base = mfoil.mfoil(naca=naca, npanel=npanel)
    M_base.param.verb = 0
    M_base.param.doplot = False
    recompute_geom(M_base)
    M_base.setoper(alpha=alpha, Ma=Ma)
    M_base.solve()
    N = M_base.foil.N
    dcl_dx = np.zeros((2, N))
    dcm_dx = np.zeros((2, N))

    for node in test_nodes:
        for dim in range(2):
            # Create object OUTSIDE complex mode (scipy needs real arrays)
            M = mfoil.mfoil(naca=naca, npanel=npanel)
            M.param.verb = 0
            M.param.doplot = False
            M.foil.x = M.foil.x.astype(complex)
            M.foil.x[dim, node] += 1j * CS_STEP
            recompute_geom(M)
            with complex_step_mode():
                M.setoper(alpha=alpha, Ma=Ma)
                M.solve()

            dcl_dx[dim, node] = np.imag(M.post.cl) / CS_STEP
            dcm_dx[dim, node] = np.imag(M.post.cm) / CS_STEP

            print(f"    inviscid shape node={node} dim={dim}: "
                  f"dcl={dcl_dx[dim, node]:.8e}, dcm={dcm_dx[dim, node]:.8e}")

    return M_base, dcl_dx, dcm_dx


def cs_viscous_shape(naca, npanel, alpha, Re, Ma, test_nodes, onset='cubic',
                     xtrip=None, cf_floor=False):
    """
    Complex-step shape sensitivities for viscous case.
    Returns dcl_dx, dcd_dx, dcm_dx [2, N+1] (nonzero only at test_nodes).
    """
    # IMPORTANT: M_base must use the SAME foil.s convention as the perturbation
    # calls below (which apply recompute_geom for complex-step compatibility).
    # mfoil's default uses scipy spline arc-length which differs from the
    # cumulative-panel-norm by ~1e-5; if M_base uses spline while perturbations
    # use cumulative, the saved cl_base lives at a different fixed point than
    # what CS evaluates, and scirocco's analytic gradient (which lines up with
    # the CS path) can't be checked against a primal that doesn't match.
    M_base = mfoil.mfoil(naca=naca, npanel=npanel)
    M_base.param.verb = 0
    M_base.param.doplot = False
    M_base.param.onset = onset
    M_base.param.xtrip = xtrip
    M_base.param.cf_laminar_floor = cf_floor
    recompute_geom(M_base)  # cumulative foil.s, same path as perturbations
    M_base.setoper(alpha=alpha, Re=Re, Ma=Ma)
    M_base.solve()
    N = M_base.foil.N
    dcl_dx = np.zeros((2, N))
    dcd_dx = np.zeros((2, N))
    dcm_dx = np.zeros((2, N))

    for node in test_nodes:
        for dim in range(2):
            M = mfoil.mfoil(naca=naca, npanel=npanel)
            M.param.verb = 0
            M.param.doplot = False
            M.param.onset = onset
            M.param.xtrip = xtrip
            M.param.cf_laminar_floor = cf_floor
            M.foil.x = M.foil.x.astype(complex)
            M.foil.x[dim, node] += 1j * CS_STEP
            recompute_geom(M)
            with complex_step_mode():
                M.setoper(alpha=alpha, Re=Re, Ma=Ma)
                M.solve()


            dcl_dx[dim, node] = np.imag(M.post.cl) / CS_STEP
            dcd_dx[dim, node] = np.imag(M.post.cd) / CS_STEP
            dcm_dx[dim, node] = np.imag(M.post.cm) / CS_STEP

            print(f"    viscous shape node={node} dim={dim}: "
                  f"dcl={dcl_dx[dim, node]:.8e}, dcd={dcd_dx[dim, node]:.8e}, "
                  f"dcm={dcm_dx[dim, node]:.8e}")

    return M_base, dcl_dx, dcd_dx, dcm_dx


# ---------------------------------------------------------------------------
# Validation: verify complex-step against FD before generating reference data
# ---------------------------------------------------------------------------
def validate_inviscid():
    """Complex-step dcl/dalpha for NACA 2412 should match FD to ~13 digits."""
    print("Validation: inviscid dcl/dalpha (NACA 2412, 80 panels, alpha=5)")

    # Complex-step
    dcl_cs, dcm_cs = cs_inviscid_scalar('2412', 80, 5.0, 0.0, 'alpha', CS_STEP)

    # Finite difference
    ep = 1e-7
    M_p = solve_inviscid('2412', 80, 5.0 + ep)
    M_m = solve_inviscid('2412', 80, 5.0 - ep)
    dcl_fd = (M_p.post.cl - M_m.post.cl) / (2.0 * ep)
    dcm_fd = (M_p.post.cm - M_m.post.cm) / (2.0 * ep)

    err_cl = abs(dcl_cs - dcl_fd) / max(abs(dcl_fd), 1e-30)
    err_cm = abs(dcm_cs - dcm_fd) / max(abs(dcm_fd), 1e-30)
    print(f"  dcl/dalpha: CS={dcl_cs:.15e}, FD={dcl_fd:.15e}, rel_err={err_cl:.2e}")
    print(f"  dcm/dalpha: CS={dcm_cs:.15e}, FD={dcm_fd:.15e}, rel_err={err_cm:.2e}")

    if err_cl > 1e-6:
        raise RuntimeError(f"Inviscid validation FAILED: dcl rel_err={err_cl:.2e}")
    print("  PASSED")
    return True


def validate_viscous():
    """Complex-step dcl/dalpha for NACA 0012 Re=1e6 should match FD to ~10 digits."""
    print("Validation: viscous dcl/dalpha (NACA 0012, 100 panels, alpha=3, Re=1e6)")

    # Complex-step
    dcl_cs, dcd_cs, dcm_cs = cs_viscous_scalar('0012', 100, 3.0, 1e6, 0.0, 'alpha', CS_STEP)

    # Finite difference
    ep = 5e-5
    M_p = solve_viscous('0012', 100, 3.0 + ep, 1e6)
    M_m = solve_viscous('0012', 100, 3.0 - ep, 1e6)
    dcl_fd = (M_p.post.cl - M_m.post.cl) / (2.0 * ep)
    dcd_fd = (M_p.post.cd - M_m.post.cd) / (2.0 * ep)

    err_cl = abs(dcl_cs - dcl_fd) / max(abs(dcl_fd), 1e-30)
    err_cd = abs(dcd_cs - dcd_fd) / max(abs(dcd_fd), 1e-30)
    print(f"  dcl/dalpha: CS={dcl_cs:.15e}, FD={dcl_fd:.15e}, rel_err={err_cl:.2e}")
    print(f"  dcd/dalpha: CS={dcd_cs:.15e}, FD={dcd_fd:.15e}, rel_err={err_cd:.2e}")

    # FD has O(ep^2) error (~1e-5 to 1e-3 depending on smoothness).
    # CS should be more accurate. Allow 5% difference (FD is the less accurate one).
    if err_cl > 0.05:
        raise RuntimeError(f"Viscous validation FAILED: dcl rel_err={err_cl:.2e}")
    print("  PASSED")
    return True


# ---------------------------------------------------------------------------
# Main: generate reference data
# ---------------------------------------------------------------------------
def write_viscous_shape_group(f, name, M, test_nodes, dcl_dx, dcd_dx, dcm_dx):
    """Write one viscous shape-sensitivity group (base state + CS partials)."""
    g = f.create_group(name)
    g.create_dataset('x_base', data=np.real(M.foil.x[:, :M.foil.N]))
    g.create_dataset('cl_base', data=float(np.real(M.post.cl)))
    g.create_dataset('cd_base', data=float(np.real(M.post.cd)))
    g.create_dataset('cm_base', data=float(np.real(M.post.cm)))
    g.create_dataset('test_nodes', data=np.array(test_nodes, dtype=np.int64))
    g.create_dataset('dcl_dx', data=dcl_dx)
    g.create_dataset('dcd_dx', data=dcd_dx)
    g.create_dataset('dcm_dx', data=dcm_dx)


# Issue #22: on NACA 2412 at alpha=4, Re=1e6, 160 panels, the upper-surface
# amplification starts growing between nodes 89 and 90. The onset ramp in
# the damping rate is narrow in Re_theta, so the shape gradient swings in
# sign from node to node there; refining the paneling does not smooth it.
# CS shows the same spike, so it belongs to the
# discrete model; these nodes pin the adjoint to it.
ONSET_CASE = dict(naca='2412', npanel=160, alpha=4.0, Re=1e6, Ma=0.0)
ONSET_NODES = list(range(85, 94))


def write_onset_group(f):
    """CS shape sensitivities around the upper-surface amplification onset."""
    c = ONSET_CASE
    print(f"\n--- viscous_shape_onset (NACA {c['naca']}, {c['npanel']} panels, "
          f"alpha={c['alpha']}, Re={c['Re']:.0e}) ---")
    M, dcl_dx, dcd_dx, dcm_dx = cs_viscous_shape(
        c['naca'], c['npanel'], c['alpha'], c['Re'], c['Ma'], ONSET_NODES)
    write_viscous_shape_group(f, 'viscous_shape_onset', M, ONSET_NODES,
                              dcl_dx, dcd_dx, dcm_dx)


# Opt-in envelope onset model (OnsetModel::Envelope / param.onset='envelope').
# Shape nodes straddle the upper-surface amplification onset (85-93) plus one
# lower- and one upper-surface node away from it.
ENVELOPE_CASE = dict(naca='2412', npanel=160, alpha=4.0, Re=1e6)
ENVELOPE_NODES = [40] + list(range(85, 94)) + [120]


def write_envelope_groups(f):
    """CS references for the envelope onset model: shape, alpha, Re, Ma."""
    c = ENVELOPE_CASE
    print("\n--- envelope_shape (NACA 2412, 160 panels, alpha=4, Re=1e6) ---")
    M, dcl_dx, dcd_dx, dcm_dx = cs_viscous_shape(
        c['naca'], c['npanel'], c['alpha'], c['Re'], 0.0, ENVELOPE_NODES,
        onset='envelope')
    write_viscous_shape_group(f, 'envelope_shape', M, ENVELOPE_NODES,
                              dcl_dx, dcd_dx, dcm_dx)

    print("\n--- envelope_scalar (alpha, Re) ---")
    g = f.create_group('envelope_scalar')
    for var, suffix in (('alpha', 'dalpha'), ('Re', 'dRe')):
        dcl, dcd, dcm = cs_viscous_scalar(c['naca'], c['npanel'], c['alpha'],
                                          c['Re'], 0.0, var, CS_STEP,
                                          onset='envelope')
        print(f"    d/{var}: dcl={dcl:.10e}, dcd={dcd:.10e}, dcm={dcm:.10e}")
        g.create_dataset('dcl_' + suffix, data=dcl)
        g.create_dataset('dcd_' + suffix, data=dcd)
        g.create_dataset('dcm_' + suffix, data=dcm)

    print("\n--- envelope_mach (Ma=0.3) ---")
    dcl, dcd, dcm = cs_viscous_scalar(c['naca'], c['npanel'], c['alpha'],
                                      c['Re'], 0.3, 'Ma', CS_STEP,
                                      onset='envelope')
    print(f"    d/dMa: dcl={dcl:.10e}, dcd={dcd:.10e}, dcm={dcm:.10e}")
    g = f.create_group('envelope_mach')
    g.create_dataset('dcl_dMa', data=dcl)
    g.create_dataset('dcd_dMa', data=dcd)
    g.create_dataset('dcm_dMa', data=dcm)


CAPTURED_CASES = [(4.0, 1e6), (0.0, 1e6), (4.0, 2e5), (8.0, 2e5), (4.0, 3e6)]


def write_laminar_te_re_group(f):
    """CS d/dRe with one surface still laminar at the trailing edge, so the
    wake initialization takes its shear stress from the Re-dependent
    transition closure (get_cttr) rather than from the state."""
    c = ENVELOPE_CASE
    print("\n--- viscous_re_laminar_te (NACA 2412, 160 panels, alpha=4, Re=1e6) ---")
    M = solve_viscous(c['naca'], c['npanel'], c['alpha'], c['Re'])
    if all(M.vsol.turb[Is[-1]] for Is in M.vsol.Is[:2]):
        raise RuntimeError('expected a laminar trailing edge on one surface')
    dcl, dcd, dcm = cs_viscous_scalar(c['naca'], c['npanel'], c['alpha'],
                                      c['Re'], 0.0, 'Re', CS_STEP)
    print(f"    d/dRe: dcl={dcl:.10e}, dcd={dcd:.10e}, dcm={dcm:.10e}")
    g = f.create_group('viscous_re_laminar_te')
    g.create_dataset('dcl_dRe', data=dcl)
    g.create_dataset('dcd_dRe', data=dcd)
    g.create_dataset('dcm_dRe', data=dcm)


# Forced transition (mfoil_cs.find_trip; Scirocco Param::xtrip_lower/upper).
# xtrip is [lower, upper]. Free transition would be at x/c = 0.40 on the
# upper surface and absent on the lower one, so both trips act.
TRIP_CASE = dict(naca='2412', npanel=160, alpha=4.0, Re=1e6)
TRIP_X = [0.6, 0.2]


def trip_nodes(M):
    """The two nodes around each transition point (trip or free) and their
    neighbours, plus one node on each surface away from them."""
    nodes = set()
    for si in range(2):
        Is = M.vsol.Is[si]
        i, _, at_node = mfoil_cs.find_trip(M, si)
        if i < 0:
            turb = [k for k in range(len(Is)) if M.vsol.turb[Is[k]]]
            if not turb:
                continue
            i = turb[0]
        elif at_node or not M.vsol.turb[Is[i]]:
            raise RuntimeError('expected a trip inside an interval')
        nodes.update(int(Is[k]) for k in range(i - 2, i + 2))
    nodes.update((20, 140))
    return sorted(nodes)


def write_trip_groups(f, prefix='trip', xtrip=None, cf_floor=False):
    """CS references for forced transition: shape, alpha, Re, Ma."""
    c = TRIP_CASE
    xtrip = TRIP_X if xtrip is None else xtrip
    shape, scalar = prefix + '_shape', prefix + '_scalar'
    print(f"\n--- {shape} (NACA 2412, 160 panels, alpha=4, Re=1e6, "
          f"xtrip lower/upper = {xtrip}, cf floor = {cf_floor}) ---")
    M0 = mfoil.mfoil(naca=c['naca'], npanel=c['npanel'])
    M0.param.verb = 0
    M0.param.doplot = False
    M0.param.xtrip = xtrip
    M0.param.cf_laminar_floor = cf_floor
    recompute_geom(M0)
    M0.setoper(alpha=c['alpha'], Re=c['Re'], Ma=0.0)
    M0.solve()
    nodes = trip_nodes(M0)
    M, dcl_dx, dcd_dx, dcm_dx = cs_viscous_shape(
        c['naca'], c['npanel'], c['alpha'], c['Re'], 0.0, nodes, xtrip=xtrip,
        cf_floor=cf_floor)
    write_viscous_shape_group(f, shape, M, nodes, dcl_dx, dcd_dx, dcm_dx)
    f[shape].create_dataset('xtrip', data=np.array(xtrip))
    f[shape].create_dataset('xt', data=np.real(M.vsol.Xt[:, 1]))

    print(f"\n--- {scalar} (alpha, Re; Ma at Ma=0.3) ---")
    g = f.create_group(scalar)
    for var, suffix, Ma in (('alpha', 'dalpha', 0.0), ('Re', 'dRe', 0.0),
                            ('Ma', 'dMa', 0.3)):
        dcl, dcd, dcm = cs_viscous_scalar(c['naca'], c['npanel'], c['alpha'],
                                          c['Re'], Ma, var, CS_STEP,
                                          xtrip=xtrip, cf_floor=cf_floor)
        print(f"    d/{var}: dcl={dcl:.10e}, dcd={dcd:.10e}, dcm={dcm:.10e}")
        g.create_dataset('dcl_' + suffix, data=dcl)
        g.create_dataset('dcd_' + suffix, data=dcd)
        g.create_dataset('dcm_' + suffix, data=dcm)


# Laminar floor on the turbulent skin friction (Param::cf_laminar_floor,
# mfoil_cs._patched_get_cf). The lower trip is close behind the stagnation
# point, where Re_theta is small; the upper surface transitions freely
# behind a separation bubble. The floor acts in both places.
FLOOR_X = [0.01, 1.0]


def write_floor_groups(f):
    write_trip_groups(f, prefix='cf_floor', xtrip=FLOOR_X, cf_floor=True)


def write_captured_group(f):
    """Captured-transition primal (mfoil_captured) on cumulative-s geometry."""
    import mfoil_captured
    print("\n--- captured_primal (NACA 2412, 160 panels) ---")
    g = f.create_group('captured_primal')
    out = []
    for alpha, Re in CAPTURED_CASES:
        M = mfoil.mfoil(naca='2412', npanel=160)
        M.param.verb = 0
        M.param.doplot = False
        recompute_geom(M)
        M.setoper(alpha=alpha, Re=Re)
        if not mfoil_captured.solve_captured(M):
            raise RuntimeError('captured solve failed at %g, %g' % (alpha, Re))
        out.append([alpha, Re, M.post.cl, M.post.cd, M.post.cm])
        print(f"    alpha={alpha} Re={Re:.0e}: cl={M.post.cl:.12f} "
              f"cd={M.post.cd:.12f} cm={M.post.cm:.12f}")
        x_base = np.real(M.foil.x[:, :M.foil.N])
    g.create_dataset('x_base', data=x_base)
    g.create_dataset('cases', data=np.array(out))  # alpha, Re, cl, cd, cm


CAPTURED_NODES = [40, 95, 110, 113, 116, 118, 120, 123, 140]


def captured_cs_solve(alpha, Re, Ma, node=None, dim=None, var=None,
                      xtrip=None):
    """Captured solve in complex_step_mode with one complex perturbation.

    Converges the real part to 1e-11, then runs extra Newton steps so the
    imaginary (derivative) part converges as well. Returns the mfoil object.
    """
    import mfoil_captured as mcap
    M = mfoil.mfoil(naca='2412', npanel=160)
    M.param.verb = 0
    M.param.doplot = False
    M.param.xtrip = xtrip
    if node is not None:
        M.foil.x = M.foil.x.astype(complex)
        M.foil.x[dim, node] += 1j * CS_STEP
    recompute_geom(M)
    kw = dict(alpha=alpha, Re=Re, Ma=Ma)
    if var is not None:
        kw[var] = kw[var] + 1j * CS_STEP
    with complex_step_mode():
        M.setoper(**kw)
        M.solve()
        mcap.convert_state(M)
        extra = None
        for it in range(60):
            mcap.build_glob_sys(M)
            if extra is None and np.linalg.norm(np.real(M.glob.R)) < 1e-11:
                extra = 3
            if extra is not None:
                if extra == 0:
                    break
                extra -= 1
            mcap.solve_glob(M)
            mcap.update_state(M)
            mfoil.stagpoint_move(M)
        if extra is None or extra > 0:
            raise RuntimeError('captured CS solve did not converge')
        U5 = M.glob.U
        M.glob.U = U5[0:4]
        mfoil.calc_force(M)
        M.glob.U = U5
    return M


def write_captured_cs_groups(f, prefix='captured', xtrip=None):
    """CS references for the captured model: shape, alpha, Re, Ma.

    With xtrip ([lower, upper]) the trips act, and the shape nodes are those
    around the trips (trip_nodes)."""
    print(f"\n--- {prefix}_shape (NACA 2412, 160, alpha=4, Re=1e6, "
          f"xtrip={xtrip}) ---")
    M0 = captured_cs_solve(4.0, 1e6, 0.0, xtrip=xtrip)
    N = M0.foil.N
    nodes = CAPTURED_NODES if xtrip is None else trip_nodes(M0)
    dcl = np.zeros((2, N)); dcd = np.zeros((2, N)); dcm = np.zeros((2, N))
    for node in nodes:
        for dim in range(2):
            M = captured_cs_solve(4.0, 1e6, 0.0, node, dim, xtrip=xtrip)
            dcl[dim, node] = np.imag(M.post.cl) / CS_STEP
            dcd[dim, node] = np.imag(M.post.cd) / CS_STEP
            dcm[dim, node] = np.imag(M.post.cm) / CS_STEP
            print(f"    node={node} dim={dim}: dcl={dcl[dim, node]:.8e} "
                  f"dcd={dcd[dim, node]:.8e} dcm={dcm[dim, node]:.8e}")
    g = f.create_group(prefix + '_shape')
    g.create_dataset('x_base', data=np.real(M0.foil.x[:, :N]))
    g.create_dataset('cl_base', data=float(np.real(M0.post.cl)))
    g.create_dataset('cd_base', data=float(np.real(M0.post.cd)))
    g.create_dataset('test_nodes', data=np.array(nodes, dtype=np.int64))
    if xtrip is not None:
        g.create_dataset('xtrip', data=np.array(xtrip))
    g.create_dataset('dcl_dx', data=dcl)
    g.create_dataset('dcd_dx', data=dcd)
    g.create_dataset('dcm_dx', data=dcm)

    print(f"\n--- {prefix}_scalar ---")
    g = f.create_group(prefix + '_scalar')
    for var, Ma, suffix in (('alpha', 0.0, 'dalpha'), ('Re', 0.0, 'dRe'),
                            ('Ma', 0.3, 'dMa')):
        M = captured_cs_solve(4.0, 1e6, Ma, var=var, xtrip=xtrip)
        vals = [np.imag(getattr(M.post, k)) / CS_STEP for k in ('cl', 'cd', 'cm')]
        print(f"    d/{var}: dcl={vals[0]:.10e} dcd={vals[1]:.10e} dcm={vals[2]:.10e}")
        for k, v in zip(('dcl_', 'dcd_', 'dcm_'), vals):
            g.create_dataset(k + suffix, data=v)


def main():
    if '--laminar-te-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            if 'viscous_re_laminar_te' in f:
                del f['viscous_re_laminar_te']
            write_laminar_te_re_group(f)
        print(f"Wrote viscous_re_laminar_te/ to {OUTPUT_FILE}")
        return

    if '--trip-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            for name in ('trip_shape', 'trip_scalar'):
                if name in f:
                    del f[name]
            write_trip_groups(f)
        print(f"Wrote trip_shape/ and trip_scalar/ to {OUTPUT_FILE}")
        return

    if '--captured-trip-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            for name in ('captured_trip_shape', 'captured_trip_scalar'):
                if name in f:
                    del f[name]
            write_captured_cs_groups(f, 'captured_trip', TRIP_X)
        print(f"Wrote captured_trip_* groups to {OUTPUT_FILE}")
        return

    if '--cf-floor-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            for name in ('cf_floor_shape', 'cf_floor_scalar'):
                if name in f:
                    del f[name]
            write_floor_groups(f)
        print(f"Wrote cf_floor_shape/ and cf_floor_scalar/ to {OUTPUT_FILE}")
        return

    if '--captured-cs-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            for name in ('captured_shape', 'captured_scalar'):
                if name in f:
                    del f[name]
            write_captured_cs_groups(f)
        print(f"Wrote captured_shape/ and captured_scalar/ to {OUTPUT_FILE}")
        return

    if '--captured-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            if 'captured_primal' in f:
                del f['captured_primal']
            write_captured_group(f)
        print(f"Wrote captured_primal/ to {OUTPUT_FILE}")
        return

    args = sys.argv[1:]
    if '--onset-only' in args or '--envelope-only' in args:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            if '--onset-only' in args:
                if 'viscous_shape_onset' in f:
                    del f['viscous_shape_onset']
                write_onset_group(f)
                print(f"Wrote viscous_shape_onset/ to {OUTPUT_FILE}")
            if '--envelope-only' in args:
                for name in ('envelope_shape', 'envelope_scalar',
                             'envelope_mach'):
                    if name in f:
                        del f[name]
                write_envelope_groups(f)
                print(f"Wrote envelope_* groups to {OUTPUT_FILE}")
        return

    t0 = time.time()

    # --- Validate first ---
    print("=" * 60)
    print("Phase 1: Validation")
    print("=" * 60)
    validate_inviscid()
    validate_viscous()

    print()
    print("=" * 60)
    print("Phase 2: Generating reference data")
    print("=" * 60)

    with h5py.File(OUTPUT_FILE, 'w') as f:
        for source in ('mfoil.py', 'mfoil_cs.py', 'complex_spline.py'):
            source_path = os.path.join(os.path.dirname(__file__), '../../mfoil', source)
            with open(source_path, 'rb') as source_file:
                f.attrs[source+'_sha256'] = hashlib.sha256(source_file.read()).hexdigest()

        # ----- inviscid_shape -----
        print("\n--- inviscid_shape (NACA 2412, 80 panels, alpha=5) ---")
        inv_test_nodes = [10, 30, 50, 70]
        M_inv, dcl_dx_inv, dcm_dx_inv = cs_inviscid_shape(
            '2412', 80, 5.0, 0.0, inv_test_nodes)

        g = f.create_group('inviscid_shape')
        g.create_dataset('test_nodes', data=np.array(inv_test_nodes, dtype=np.int64))
        g.create_dataset('x_base', data=np.real(M_inv.foil.x[:, :M_inv.foil.N]))
        g.create_dataset('cl_base', data=float(M_inv.post.cl))
        g.create_dataset('cm_base', data=float(M_inv.post.cm))
        g.create_dataset('dcl_dx', data=dcl_dx_inv)
        g.create_dataset('dcm_dx', data=dcm_dx_inv)

        # ----- inviscid_alpha -----
        print("\n--- inviscid_alpha (NACA 2412, 80 panels, alpha=5) ---")
        dcl_da_inv, dcm_da_inv = cs_inviscid_scalar('2412', 80, 5.0, 0.0, 'alpha', CS_STEP)
        print(f"    dcl/dalpha={dcl_da_inv:.10e}, dcm/dalpha={dcm_da_inv:.10e}")

        g = f.create_group('inviscid_alpha')
        g.create_dataset('dcl_dalpha', data=dcl_da_inv)
        g.create_dataset('dcm_dalpha', data=dcm_da_inv)

        # ----- inviscid_mach -----
        print("\n--- inviscid_mach (NACA 0012, 100 panels, alpha=3, Ma=0.3) ---")
        dcl_dMa_inv, dcm_dMa_inv = cs_inviscid_scalar('0012', 100, 3.0, 0.3, 'Ma', CS_STEP)
        M_inv_ma = solve_inviscid('0012', 100, 3.0, 0.3)
        print(f"    dcl/dMa={dcl_dMa_inv:.10e}, dcm/dMa={dcm_dMa_inv:.10e}")

        g = f.create_group('inviscid_mach')
        g.create_dataset('cl_base', data=float(M_inv_ma.post.cl))
        g.create_dataset('dcl_dMa', data=dcl_dMa_inv)
        g.create_dataset('dcm_dMa', data=dcm_dMa_inv)

        # ----- viscous_shape -----
        print("\n--- viscous_shape (NACA 0012, 100 panels, alpha=3, Re=1e6) ---")
        visc_test_nodes = [20, 50, 80]
        M_visc, dcl_dx_visc, dcd_dx_visc, dcm_dx_visc = cs_viscous_shape(
            '0012', 100, 3.0, 1e6, 0.0, visc_test_nodes)

        write_viscous_shape_group(f, 'viscous_shape', M_visc, visc_test_nodes,
                                  dcl_dx_visc, dcd_dx_visc, dcm_dx_visc)

        write_onset_group(f)

        # ----- viscous_alpha -----
        print("\n--- viscous_alpha (NACA 0012, 100 panels, alpha=3, Re=1e6) ---")
        dcl_da_v, dcd_da_v, dcm_da_v = cs_viscous_scalar(
            '0012', 100, 3.0, 1e6, 0.0, 'alpha', CS_STEP)
        print(f"    dcl/dalpha={dcl_da_v:.10e}, dcd/dalpha={dcd_da_v:.10e}, dcm/dalpha={dcm_da_v:.10e}")

        g = f.create_group('viscous_alpha')
        g.create_dataset('dcl_dalpha', data=dcl_da_v)
        g.create_dataset('dcd_dalpha', data=dcd_da_v)
        g.create_dataset('dcm_dalpha', data=dcm_da_v)

        # ----- viscous_mach -----
        print("\n--- viscous_mach (NACA 0012, 100 panels, alpha=3, Re=1e6, Ma=0.3) ---")
        dcl_dMa_v, dcd_dMa_v, dcm_dMa_v = cs_viscous_scalar(
            '0012', 100, 3.0, 1e6, 0.3, 'Ma', CS_STEP)
        M_visc_ma = solve_viscous('0012', 100, 3.0, 1e6, 0.3)
        print(f"    dcl/dMa={dcl_dMa_v:.10e}, dcd/dMa={dcd_dMa_v:.10e}, dcm/dMa={dcm_dMa_v:.10e}")

        g = f.create_group('viscous_mach')
        g.create_dataset('cl_base', data=float(M_visc_ma.post.cl))
        g.create_dataset('dcl_dMa', data=dcl_dMa_v)
        g.create_dataset('dcd_dMa', data=dcd_dMa_v)
        g.create_dataset('dcm_dMa', data=dcm_dMa_v)

        # ----- viscous_re -----
        print("\n--- viscous_re (NACA 0012, 100 panels, alpha=3, Re=1e6) ---")
        dcl_dRe, dcd_dRe, dcm_dRe = cs_viscous_scalar(
            '0012', 100, 3.0, 1e6, 0.0, 'Re', CS_STEP)
        print(f"    dcl/dRe={dcl_dRe:.10e}, dcd/dRe={dcd_dRe:.10e}, dcm/dRe={dcm_dRe:.10e}")

        g = f.create_group('viscous_re')
        g.create_dataset('dcl_dRe', data=dcl_dRe)
        g.create_dataset('dcd_dRe', data=dcd_dRe)
        g.create_dataset('dcm_dRe', data=dcm_dRe)

        write_envelope_groups(f)
        write_laminar_te_re_group(f)
        write_captured_group(f)
        write_captured_cs_groups(f)
        write_trip_groups(f)
        write_floor_groups(f)

    elapsed = time.time() - t0
    print(f"\nDone. Wrote {OUTPUT_FILE}")
    print(f"Total time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
