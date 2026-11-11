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
  inviscid_shape/    NACA 2412, 80 panels, alpha=5
  inviscid_alpha/    NACA 2412, 80 panels, alpha=5
  inviscid_mach/     NACA 0012, 100 panels, alpha=3, Ma=0.3
  envelope_shape/    NACA 2412, 160 panels, alpha=4, Re=1e6, onset='envelope'
  envelope_scalar/   same case: d/dalpha and d/dRe
  envelope_mach/     same case at Ma=0.3: d/dMa

Usage:
  cd mfoil
  python ../scirocco/tests/extract_sensitivity_reference.py
  python ../scirocco/tests/extract_sensitivity_reference.py --envelope-only

--envelope-only (re)writes just the envelope_* groups into the existing file.
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


def cs_viscous_scalar(naca, npanel, alpha, Re, Ma, var, step, onset='cubic'):
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


def cs_viscous_shape(naca, npanel, alpha, Re, Ma, test_nodes, onset='cubic'):
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
    g = f.create_group('envelope_shape')
    g.create_dataset('x_base', data=np.real(M.foil.x[:, :M.foil.N]))
    g.create_dataset('cl_base', data=float(np.real(M.post.cl)))
    g.create_dataset('cd_base', data=float(np.real(M.post.cd)))
    g.create_dataset('cm_base', data=float(np.real(M.post.cm)))
    g.create_dataset('test_nodes', data=np.array(ENVELOPE_NODES, dtype=np.int64))
    g.create_dataset('dcl_dx', data=dcl_dx)
    g.create_dataset('dcd_dx', data=dcd_dx)
    g.create_dataset('dcm_dx', data=dcm_dx)

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


def main():
    if '--envelope-only' in sys.argv[1:]:
        with h5py.File(OUTPUT_FILE, 'a') as f:
            for name in ('envelope_shape', 'envelope_scalar', 'envelope_mach'):
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

        g = f.create_group('viscous_shape')
        g.create_dataset('x_base', data=np.real(M_visc.foil.x[:, :M_visc.foil.N]))
        g.create_dataset('cl_base', data=float(np.real(M_visc.post.cl)))
        g.create_dataset('cd_base', data=float(np.real(M_visc.post.cd)))
        g.create_dataset('cm_base', data=float(np.real(M_visc.post.cm)))
        g.create_dataset('test_nodes', data=np.array(visc_test_nodes, dtype=np.int64))
        g.create_dataset('dcl_dx', data=dcl_dx_visc)
        g.create_dataset('dcd_dx', data=dcd_dx_visc)
        g.create_dataset('dcm_dx', data=dcm_dx_visc)

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

    elapsed = time.time() - t0
    print(f"\nDone. Wrote {OUTPUT_FILE}")
    print(f"Total time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
