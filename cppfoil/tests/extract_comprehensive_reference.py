#!/usr/bin/env python3
"""
Generate comprehensive CS reference data for debugging viscous shape
sensitivity adjoint.

For EVERY foil node k (0..N-1) and dim (0,1):
  - CS-perturb foil.x[dim, k]
  - Rebuild everything: foil.s, build_gamma, stagpoint_find, build_wake,
    identify_surfaces, set_wake_gap, calc_ue_m
  - Evaluate full residual R = [R_BL; R_ue] at fixed BL state U
  - Store dR[:, dim, k] = Im(R_pert) / step

Also stores:
  - Force partials (dcl, dcd, dcm w.r.t. foil.x at fixed gamma)
  - Lambda vectors from mfoil adjoint solve (for independent verification)
  - Full BL state U and converged gamma

Usage:
  cd ProjectTrysil/temp
  python ../cppfoil/tests/extract_comprehensive_reference.py
"""

import sys
import os
import copy
import time
import numpy as np
import h5py
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'temp'))
import mfoil_cs
from mfoil_cs import complex_step_mode
import mfoil

CS_STEP = 1e-20
OUTPUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'comprehensive_reference.h5')

NACA = '0012'
NPANEL = 100
ALPHA = 3.0
RE = 1e6
MA = 0.0


def solve_base():
    """Solve viscous base case, then switch foil.s to cumulative panel norms
    so the base matches the perturbation-chain convention. UE-BASED xi
    (Newton-converged) — no post-solve stagpoint_find, matching cppfoil's
    production SolveViscous (commit 4051c2d removed the gamma-xi swap)."""
    M = mfoil.mfoil(naca=NACA, npanel=NPANEL)
    M.param.verb = 0
    M.param.doplot = False
    M.setoper(alpha=ALPHA, Re=RE, Ma=MA)
    M.solve()
    # Replace foil.s with cumulative-panel-norm value (matches recompute_geom
    # used during CS perturbations). Then re-run the geom chain + refresh
    # xi/sstag at the NEW cumulative foil.s by calling stagpoint_move at
    # the converged state. This matches cppfoil's natural state after
    # SolveViscous at cumulative foil.s (xi from last Newton stagpoint_move
    # uses the current foil.s).
    recompute_geom(M)
    mfoil.build_gamma(M, ALPHA)
    mfoil.build_wake(M)
    mfoil.identify_surfaces(M)
    mfoil.set_wake_gap(M)
    mfoil.calc_ue_m(M)
    # Refresh ue-based sstag and xi for the NEW cumulative foil.s
    mfoil.stagpoint_move(M)
    return M


def recompute_geom(M):
    """Recompute arc-length from (complex) foil.x."""
    x = M.foil.x
    N = M.foil.N
    ds = np.sqrt((x[0, 1:] - x[0, :-1])**2 + (x[1, 1:] - x[1, :-1])**2)
    s = np.zeros(N, dtype=x.dtype)
    s[1:] = np.cumsum(ds)
    M.foil.s = s


def compute_R_ue(M):
    """Compute R_ue = ue - ueinv - sgnue * ue_sigma * phi at fixed U."""
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw

    ue = M.glob.U[3, :]
    ds_bl = M.glob.U[1, :]
    phi = M.vsol.sigma_m @ (ds_bl * ue)
    ueinv = mfoil.get_ueinv(M)
    sgnue = np.concatenate((M.isol.sgnue, np.ones(Nw)))
    return ue - ueinv - sgnue * (M.vsol.ue_sigma @ phi)


def compute_R_BL(M):
    """Compute R_BL [3*Nsys] by calling build_glob_sys.

    Pre-promotes all arrays to complex dtype before calling, then extracts
    M.glob.R (the BL part of the residual).
    """
    Nsys = M.glob.Nsys

    # Pre-promote arrays that build_glob_sys allocates internally
    M.glob.R = np.zeros(3 * Nsys, dtype=complex)
    M.glob.R_U = sparse.lil_matrix((3 * Nsys, 4 * Nsys), dtype=complex)
    M.glob.R_x = sparse.lil_matrix((3 * Nsys, Nsys), dtype=complex)
    M.glob.realloc = False  # prevent re-allocation inside build_glob_sys

    mfoil.build_glob_sys(M)

    return M.glob.R.copy()


def compute_full_residual(M):
    """Compute full residual [R_BL; R_ue] at fixed U."""
    R_BL = compute_R_BL(M)
    R_ue = compute_R_ue(M)
    return np.concatenate((R_BL, R_ue))


def cs_force_partials(M_base):
    """CS force partials dcl/dx, dcd/dx, dcm/dx at fixed gamma.

    Returns dcl [2, N], dcd [2, N], dcm [2, N].
    """
    N = M_base.foil.N
    dcl = np.zeros((2, N))
    dcd = np.zeros((2, N))
    dcm = np.zeros((2, N))

    for k in range(N):
        for dim in range(2):
            M = copy.deepcopy(M_base)
            M.foil.x = M.foil.x.astype(complex)
            M.foil.x[dim, k] += 1j * CS_STEP

            with complex_step_mode():
                mfoil.calc_force(M)

            dcl[dim, k] = np.imag(M.post.cl) / CS_STEP
            dcd[dim, k] = np.imag(M.post.cd) / CS_STEP
            dcm[dim, k] = np.imag(M.post.cm) / CS_STEP

        if k % 20 == 0 or k == N - 1:
            print(f'  force partial: node={k}/{N-1}')

    return dcl, dcd, dcm


def cs_full_dRdx(M_base):
    """CS derivative of full residual w.r.t. ALL foil.x at fixed BL state.

    For each (k, dim): perturb foil.x[dim, k], rebuild everything
    (foil.s, gamma, wake, ue_m, xi, wgap), evaluate residual at fixed U.

    Returns dR_dx [4*Nsys, 2, N]
    """
    N = M_base.foil.N
    Nw = M_base.wake.N
    Nsys = N + Nw
    dR_dx = np.zeros((4 * Nsys, 2, N))

    t0 = time.time()

    for k in range(N):
        for dim in range(2):
            M = copy.deepcopy(M_base)
            M.foil.x = M.foil.x.astype(complex)
            M.foil.x[dim, k] += 1j * CS_STEP
            recompute_geom(M)

            with complex_step_mode():
                # Promote arrays to complex (base BL state stays fixed)
                M.glob.U = M.glob.U.astype(complex)
                M.isol.uewiref = M.isol.uewiref.astype(complex)
                M.vsol.sigma_m = M.vsol.sigma_m.astype(complex)

                # Rebuild everything that depends on foil.x. UE-BASED xi:
                # call stagpoint_move (not stagpoint_find) to refresh sstag
                # AND xi at the perturbed cumulative foil.s — matching
                # cppfoil's natural state after SolveViscous (commit 4051c2d
                # removed the post-Newton xi swap).
                mfoil.build_gamma(M, ALPHA)
                mfoil.build_wake(M)
                mfoil.identify_surfaces(M)
                mfoil.set_wake_gap(M)
                mfoil.calc_ue_m(M)
                mfoil.stagpoint_move(M)

                R_pert = compute_full_residual(M)
            dR_dx[:, dim, k] = np.imag(R_pert) / CS_STEP

        if k % 10 == 0 or k == N - 1:
            elapsed = time.time() - t0
            rate = (2 * (k + 1)) / max(elapsed, 0.01)
            remaining = (2 * (N - k - 1)) / max(rate, 0.01)
            print(f'  dR/dx: node={k}/{N-1}  '
                  f'|dR_BL|={np.linalg.norm(dR_dx[:3*Nsys, :, k]):.3e}  '
                  f'|dR_ue|={np.linalg.norm(dR_dx[3*Nsys:, :, k]):.3e}  '
                  f'({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)')

    return dR_dx


def compute_mfoil_lambda(M):
    """Solve for lambda vectors in mfoil using the R_V Jacobian.

    R_V^T * lambda = dJ/dV  for each objective.

    Returns dict with lambda_cl, lambda_cd, lambda_cm [4*Nsys].
    """
    Nsys = M.glob.Nsys

    # Need to build R_V. Use solve_glob's assembly but extract R_V.
    # Actually, we need the converged Jacobian. Let's call build_glob_sys
    # to get R_U and R_x, then assemble R_V as solve_glob does.

    # Rebuild BL system at converged state to get R_U
    mfoil.build_glob_sys(M)

    # Assemble R_V [4*Nsys x 4*Nsys]
    ue = M.glob.U[3, :]
    ds = M.glob.U[1, :]
    ueinv = mfoil.get_ueinv(M)
    Ids = slice(1, 4 * Nsys, 4)
    Iue = slice(3, 4 * Nsys, 4)
    docl = M.oper.givencl
    NN = 4 * Nsys + docl

    R_V = sparse.lil_matrix((NN, NN))
    R_V[0:3 * Nsys, 0:4 * Nsys] = M.glob.R_U
    I = slice(3 * Nsys, 4 * Nsys, 1)
    R_V[I, Iue] = sparse.identity(Nsys) - M.vsol.ue_m @ np.diag(ds)
    R_V[I, Ids] = -M.vsol.ue_m @ np.diag(ue)

    if docl:
        Rcla, Ru_alpha, Rcla_U = mfoil.clalpha_residual(M)
        R_V[I, 4 * Nsys] = Ru_alpha
        R_V[4 * Nsys, :] = Rcla_U

    R_V = R_V.tocsr()

    # Construct dJ/dV for cl, cd, cm
    N = M.foil.N
    # cl
    dcl_dV = np.zeros(NN)
    for i in range(N):
        dcl_dV[4 * i + 3] = M.post.cl_ue[i]

    # cd: Squire-Young formulation (dcd/dV stored in M.post.cd_U at wake end)
    dcd_dV = np.zeros(NN)
    iw = M.vsol.Is[2][-1]
    dcd_dV[4 * iw:4 * iw + 4] = M.post.cd_U

    # cm
    dcm_dV = np.zeros(NN)
    # cm depends on ue at airfoil nodes through cp
    # cm_ue is populated by calc_force if we store it
    # For now, we'll use the cl-like pattern. Actually mfoil doesn't store
    # cm_ue directly. We'll skip cm lambda for now and focus on cl and cd.

    # Solve R_V^T * lambda = dJ/dV
    lambda_cl = sparse.linalg.spsolve(R_V.T, dcl_dV)
    lambda_cd = sparse.linalg.spsolve(R_V.T, dcd_dV)

    # Also return the assembled R_V (dense) so cppfoil can diff entries.
    # Trim alpha row/col if present.
    R_V_dense = R_V.toarray()[:4 * Nsys, :4 * Nsys]

    return {
        'lambda_cl': lambda_cl[:4 * Nsys],  # trim docl row if present
        'lambda_cd': lambda_cd[:4 * Nsys],
        'R_V': R_V_dense,
    }


def fd_validate(M_base, dR_dx, validate_nodes, Nsys):
    """Validate a few CS columns against FD for sanity check."""
    print("\n  FD validation (should agree with CS to ~1e-5):")
    ep = 1e-7

    for k in validate_nodes:
        for dim in range(1):  # just dim=0 for speed
            M_p = copy.deepcopy(M_base)
            M_p.foil.x[dim, k] += ep
            recompute_geom(M_p)
            mfoil.build_gamma(M_p, ALPHA)
            mfoil.build_wake(M_p)
            mfoil.identify_surfaces(M_p)
            mfoil.set_wake_gap(M_p)
            mfoil.calc_ue_m(M_p)
            R_p = compute_full_residual(M_p)

            M_m = copy.deepcopy(M_base)
            M_m.foil.x[dim, k] -= ep
            recompute_geom(M_m)
            mfoil.build_gamma(M_m, ALPHA)
            mfoil.build_wake(M_m)
            mfoil.identify_surfaces(M_m)
            mfoil.set_wake_gap(M_m)
            mfoil.calc_ue_m(M_m)
            R_m = compute_full_residual(M_m)

            dR_fd = (R_p - R_m) / (2 * ep)
            dR_cs = dR_dx[:, dim, k]

            err = np.linalg.norm(dR_fd - dR_cs)
            ref = max(np.linalg.norm(dR_cs), 1e-30)
            print(f'    k={k} d={dim}: |CS|={ref:.3e} |FD-CS|={err:.3e} '
                  f'rel={err/ref:.3e}')


def main():
    t0 = time.time()
    print("Solving base case...")
    M = solve_base()
    N = M.foil.N
    Nw = M.wake.N
    Nsys = N + Nw
    print(f"  N={N}, Nw={Nw}, Nsys={Nsys}")
    print(f"  cl={M.post.cl:.10f}, cd={M.post.cd:.10f}, cm={M.post.cm:.10f}")

    # --- 1. Force partials ---
    print(f"\nComputing force partials via CS (all {N} nodes)...")
    dcl_partial, dcd_partial, dcm_partial = cs_force_partials(M)

    # --- 2. Full dR/dx ---
    print(f"\nComputing full dR/dx via CS (all {N} nodes, 2 dims = {2*N} evals)...")
    dR_dx = cs_full_dRdx(M)

    # --- 3. FD validation on a few nodes ---
    fd_validate(M, dR_dx, [5, 25, 50, 75, N - 5], Nsys)

    # --- 4. Lambda vectors from mfoil ---
    print("\nComputing mfoil lambda vectors...")
    try:
        lambdas = compute_mfoil_lambda(M)
        have_lambda = True
        for name, lam in lambdas.items():
            print(f"  {name}: |lambda|={np.linalg.norm(lam):.6e}")
    except Exception as e:
        print(f"  WARNING: lambda computation failed: {e}")
        have_lambda = False

    # --- 5. Cross-check: partial - lambda^T * dR/dx should match total sens ---
    if have_lambda:
        print("\n  Cross-check: partial - lambda^T @ dR/dx vs total sensitivity")
        lam_cl = lambdas['lambda_cl']
        lam_cd = lambdas['lambda_cd']
        for k in [20, 50, 80]:
            for dim in range(2):
                vjp_cl = lam_cl @ dR_dx[:, dim, k]
                vjp_cd = lam_cd @ dR_dx[:, dim, k]
                total_cl = dcl_partial[dim, k] - vjp_cl
                total_cd = dcd_partial[dim, k] - vjp_cd
                print(f'    k={k} d={dim}: dcl_total={total_cl:.6e}  '
                      f'dcd_total={total_cd:.6e}')

    # --- Write HDF5 ---
    print(f"\nWriting {OUTPUT_FILE}...")
    with h5py.File(OUTPUT_FILE, 'w') as f:
        g = f.create_group('comprehensive')
        g.attrs['naca'] = NACA
        g.attrs['npanel'] = NPANEL
        g.attrs['alpha'] = ALPHA
        g.attrs['Re'] = RE
        g.attrs['Ma'] = MA

        # Base state
        g.create_dataset('foil_x', data=np.real(M.foil.x[:, :N]))
        g.create_dataset('foil_s', data=np.real(M.foil.s))
        g.create_dataset('U', data=np.real(M.glob.U))
        g.create_dataset('gam', data=np.real(M.isol.gam))
        g.create_dataset('xi', data=np.real(M.isol.xi))
        g.create_dataset('wgap', data=np.real(M.vsol.wgap))

        # Full residual Jacobian w.r.t. foil.x [4*Nsys, 2, N]
        g.create_dataset('dR_dx', data=dR_dx)

        # BL part [3*Nsys, 2, N] and ue part [Nsys, 2, N] (split for convenience)
        g.create_dataset('dR_BL_dx', data=dR_dx[:3 * Nsys, :, :])
        g.create_dataset('dR_ue_dx', data=dR_dx[3 * Nsys:, :, :])

        # Force partials [2, N]
        g.create_dataset('dcl_partial', data=dcl_partial)
        g.create_dataset('dcd_partial', data=dcd_partial)
        g.create_dataset('dcm_partial', data=dcm_partial)

        # Base force values
        g.create_dataset('cl_base', data=float(M.post.cl))
        g.create_dataset('cd_base', data=float(M.post.cd))
        g.create_dataset('cm_base', data=float(M.post.cm))

        # Lambda vectors (if computed)
        if have_lambda:
            g.create_dataset('lambda_cl', data=lambdas['lambda_cl'])
            g.create_dataset('lambda_cd', data=lambdas['lambda_cd'])
            # Dense R_V matrix for cppfoil-vs-mfoil entry-by-entry comparison.
            # 4*Nsys x 4*Nsys, ~ 488^2 * 8 bytes = 1.9 MB; fine to store.
            g.create_dataset('R_V', data=lambdas['R_V'])

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s.")


if __name__ == '__main__':
    main()
