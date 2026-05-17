/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#include "cppfoil/core/viscous_residual_deriv.h"

#include <vector>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/inviscid_deriv.h"
#include "cppfoil/core/panel_influence_deriv.h"
#include "cppfoil/core/viscous_geom_deriv.h"
#include "cppfoil/properties/bl_properties.h"

namespace cppfoil {

// ============================================================================
// DRblDFoilX
//
//   dR_BL/d(foil.x)
//     = R_x · DXiDFoilX                       (closed-form xi chain)
//     + (wgap chain at non-initial wake stations from ResidualStation)
//     + (hTE chain at the first wake station from WakeSys)
//
// R_x is the [3*Nsys x Nsys] sparse linearisation already assembled by
// BuildGlobSys.
//
// wgap chain: ResidualStation modifies U via "U.row(1) -= Aux", so
//   ∂R_station/∂Aux(0) = -∂R_station/∂U[1] at the first column (lin_U(:, 1)),
//   ∂R_station/∂Aux(1) = -∂R_station/∂U[1] at the second column (lin_U(:, 5))
// These derivatives live in R_U at rows 3*Is[2][ip]..3*Is[2][ip]+2,
// columns 4*Is[2][ip-1]+1 and 4*Is[2][ip]+1 respectively. Aux at wake-station
// local index `i` equals wgap(i), so the chain through DWgapDFoilX is:
//   contribution to dR_BL row r = -R_U(r, 4*Is[2][ip-1]+1) * dwgap(ip-1)
//                                 -R_U(r, 4*Is[2][ip]+1)   * dwgap(ip)
// for each wake-station index ip = 1..Nw-1 (the first wake station is
// handled by WakeSys, not ResidualStation).
//
// hTE chain at WakeSys: residual is R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE),
// so ∂R(1)/∂hTE = -1, and ∂hTE/∂foil.x[k]_dim comes from GetTEInfoDeriv.
// Only the second component of the first-wake-station residual triple is
// affected.
// ============================================================================
MatrixXd DRblDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    if (state.glob.R_x.rows() != 3 * Nsys || state.glob.R_x.cols() != Nsys) {
        throw std::runtime_error(
            "DRblDFoilX: state.glob.R_x not assembled — call BuildGlobSys "
            "first");
    }

    MatrixXd dxi_dx = DXiDFoilX(state, dfoil_s_dx);    // [Nsys x 2N]
    MatrixXd dwgap_dx = DWgapDFoilX(state, dfoil_s_dx); // [Nw   x 2N]

    // ---- R_x · DXi ----
    MatrixXd dRbl = state.glob.R_x * dxi_dx;  // sparse * dense → dense

    // ---- wgap chain at wake stations ip = 1..Nw-1 ----
    //
    // ResidualStation modifies U via "U.row(1) -= Aux" AND uses Aux directly
    // in GetHw (Hw = Aux/U[0]). So Aux has TWO chain paths into R:
    //   1. U_mod[1] = U[1] - Aux → dU_mod[1]/dAux = -1, contributes -R_U[:, U1_col]
    //   2. Hw = Aux/U[0] → dHw/dAux = 1/U[0], contributes ∂R/∂Hw * 1/U[0]
    //
    // Hw appears in (Rmom, Rshape, Rlag) with coefficients (uelog, -uelog, 0).
    // Total ∂R/∂Aux = -R_U[:, U1_col] + [uelog, -uelog, 0] / U[0].
    //
    // For wake stations with cf=0, these two paths CANCEL for Rmom and Rshape
    // (the U_mod chain through H = (U[1]-Aux)/U[0] gives ±uelog/U[0] that
    // exactly cancels the Hw-direct path). Missing the Hw path produces a
    // spurious 0.5-0.9-magnitude chain term — caught at the Phase 2 test.
    const auto& Iw = state.vsol.Is[2];
    if (static_cast<int>(Iw.size()) != Nw) {
        throw std::runtime_error("DRblDFoilX: Is[wake].size() != Nw");
    }
    const SparseMatrixd& RU = state.glob.R_U;
    for (int ip = 1; ip < Nw; ++ip) {
        int Ig = 3 * Iw[ip];
        int node_prev = Iw[ip - 1];
        int node_curr = Iw[ip];
        int col_prev = 4 * node_prev + 1;
        int col_curr = 4 * node_curr + 1;

        // Compute uelog and U[0] at the two endpoints (Karman-Tsien-corrected
        // edge velocity uk via GetUk to match residual_station internals).
        auto [uk_prev, _uk_U_prev] = GetUk(state.glob.U(3, node_prev), state.param);
        auto [uk_curr, _uk_U_curr] = GetUk(state.glob.U(3, node_curr), state.param);
        double uelog = std::log(uk_curr / uk_prev);
        double U0_prev = state.glob.U(0, node_prev);
        double U0_curr = state.glob.U(0, node_curr);

        // Hw in residual is the AVERAGE: Hw_avg = 0.5*(Hw1 + Hw2). So the
        // per-endpoint Hw chain contributes 0.5*∂R/∂Hw_avg, not the full thing.
        // R coefficients for Hw_avg: +uelog (Rmom), -uelog (Rshape), 0 (Rlag).
        // Per-endpoint Hw_i: factor 0.5 from the average.
        Vector3d Hw_coef(0.5 * uelog, -0.5 * uelog, 0.0);

        for (int r = 0; r < 3; ++r) {
            double cprev = RU.coeff(Ig + r, col_prev);
            double ccurr = RU.coeff(Ig + r, col_curr);
            // Correct ∂R/∂Aux = (-R_U via U_mod chain) + (Hw chain via 1/U[0])
            double dRdAux_prev = -cprev + Hw_coef(r) / U0_prev;
            double dRdAux_curr = -ccurr + Hw_coef(r) / U0_curr;
            // Chain to foil.x: contribute (∂R/∂Aux) * dwgap_dx
            if (dRdAux_prev != 0.0) {
                dRbl.row(Ig + r) += dRdAux_prev * dwgap_dx.row(ip - 1);
            }
            if (dRdAux_curr != 0.0) {
                dRbl.row(Ig + r) += dRdAux_curr * dwgap_dx.row(ip);
            }
        }
    }

    // ---- hTE chain at the first wake station (WakeSys) ----
    // R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE), so dR(1)/dhTE = -1.
    // hTE depends on foil.x[k]_dim only for k ∈ {0, 1, N-2, N-1}.
    int Ig_wake_init = 3 * Iw[0] + 1;  // second residual of first wake station
    for (int k : {0, 1, N - 2, N - 1}) {
        for (int dim = 0; dim < 2; ++dim) {
            double dtcp, dtdp, dhTE, ddtdx;
            GetTEInfoDeriv(state.foil.x, k, dim, dtcp, dtdp, dhTE, ddtdx);
            if (dhTE != 0.0) {
                dRbl(Ig_wake_init, 2 * k + dim) += -dhTE;
            }
        }
    }

    return dRbl;
}

// ============================================================================
// ComputeRueVJP — reverse-mode VJP through R_ue chain.
//
// Returns lambda_ue^T * dR_ue/dfoil_x as a [2N] vector.
//
// R_ue[i] = U[3, i] - ueinv[i] - sgue[i] * (ue_sigma * phi)[i]
// phi = sigma_m * mass,  mass[k] = U[1, k] * U[3, k]  (mass held fixed at U)
// sgue = [sgnue (foil); ones (wake)]
//
// Implemented chain pieces (partial-analytical):
//   A1. Foil-side dueinv via dgam: -lambda^T * sgnue .* dgam/dfoil_x
//        = +AICResidualSensVJP(state, AIC^-T * lambda_g_ext, alpha)
//        where lambda_g_ext = [(sgnue .* lambda_ue[0:N]); 0] in the N+1 system.
//   C. sigma_m chain via dfoil_s_dx and dwake_s_dx:
//        beta = ue_sigma^T * (sgue .* lambda_ue)
//        alpha_panel[k] = beta[k] * phi[k] / ds[k]
//        contribution += sum_k alpha_panel[k] * d(ds[k])/d(foil.x)
//
// NOT YET IMPLEMENTED (will give wrong VJP for cases where these matter):
//   A2. Wake-side dueinv chain through wake panel velocities (duewi/dfoil_x).
//   B. dueisigma chain (panel velocity derivatives at airfoil/wake stations).
//
// These TODOs prevent full Phase 7 cl/cd/cm accuracy when wake contributions
// dominate. Foil-side adjoints with weak wake coupling should be mostly OK.
// ============================================================================
VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    int twoN = 2 * N;
    double alpha = state.oper.alpha;

    if (lambda_ue.size() != Nsys) {
        throw std::runtime_error("ComputeRueVJP: lambda_ue size != Nsys");
    }

    VectorXd vjp = VectorXd::Zero(twoN);

    // ---- A1. Foil-side dueinv chain via dgam ---------------------------------
    // lambda_g (N-vector) = sgnue[0..N) .* lambda_ue[0..N)
    // Extended to N+1 with 0 in the streamfunction slot.
    int AIC_dim = static_cast<int>(state.isol.AIC.rows());  // N + 1
    VectorXd lambda_g_ext = VectorXd::Zero(AIC_dim);
    for (int i = 0; i < N; ++i) {
        lambda_g_ext(i) = state.isol.sgnue(i) * lambda_ue(i);
    }
    // mu = AIC^{-T} * lambda_g_ext
    MatrixXd AIC_T = state.isol.AIC.transpose();
    auto lu = AIC_T.partialPivLu();
    VectorXd mu = lu.solve(lambda_g_ext);
    // Contribution: +AICResidualSensVJP(state, mu, alpha)
    //   (derivation: -lambda^T * dueinv_foil = -lambda_g^T * dgam,
    //    dgam = -(AIC)^{-1} * dR_AIC/dfoil_x, so -lambda_g^T * dgam =
    //    +(AIC^{-T} lambda_g)^T * dR_AIC = +AICResidualSensVJP.)
    vjp += AICResidualSensVJP(state, mu, alpha);

    // ---- C. dphi chain via dsigma_m ----------------------------------------
    // beta = ue_sigma^T * (sgue .* lambda_ue),  shape [Npanel = N+Nw-2]
    VectorXd sgue_lambda(Nsys);
    for (int i = 0; i < N; ++i) sgue_lambda(i) = state.isol.sgnue(i) * lambda_ue(i);
    for (int i = N; i < Nsys; ++i) sgue_lambda(i) = lambda_ue(i);
    VectorXd beta = state.vsol.ue_sigma.transpose() * sgue_lambda;  // [Npanel]

    // phi = sigma_m * mass,  mass[k] = U[1, k] * U[3, k]
    VectorXd mass(Nsys);
    for (int i = 0; i < Nsys; ++i) mass(i) = state.glob.U(1, i) * state.glob.U(3, i);
    VectorXd phi = state.vsol.sigma_m * mass;  // [Npanel]

    // Foil panels k = 0 .. N-2: ds_k = foil.s[k+1] - foil.s[k]
    // dsigma_m[k, i]/d(foil.x) = +sgnue[i]/ds_k^2 * d(ds_k)/d(foil.x)  for i = k or k+1 sign-flipped
    // dphi[k] = sum_i dsigma_m[k, i] * mass[i] = -phi[k]/ds_k * d(ds_k)
    // (Mass contribution multiplied through cancels via -1/ds^2 * (sgnue[k+1]*mass[k+1] - sgnue[k]*mass[k]) = -phi[k]/ds.)
    //
    // Total contribution to VJP for panel k:
    //   -lambda_sg^T * ue_sigma * dphi[k]
    //     = -beta[k] * dphi[k]
    //     = -beta[k] * (-phi[k]/ds_k * d(ds_k)/dfoil_x)
    //     = (beta[k] * phi[k] / ds_k) * d(ds_k)/dfoil_x
    for (int k = 0; k < N - 1; ++k) {
        double ds_k = state.foil.s(k + 1) - state.foil.s(k);
        if (ds_k < 1e-14) continue;
        double coef = beta(k) * phi(k) / ds_k;
        // d(ds_k)/dfoil_x = dfoil_s_dx.row(k+1) - dfoil_s_dx.row(k)
        vjp += coef * (dfoil_s_dx.row(k + 1).transpose() -
                       dfoil_s_dx.row(k).transpose());
    }
    // Wake panels k = N-1 .. Npanel-1: ds_w = wake.s[w+1] - wake.s[w]
    // where w = k - (N-1) is the wake panel index.
    // Same structure: need d(ds_w)/d(foil.x) = dwake_s_dx.row(w+1) - dwake_s_dx.row(w)
    MatrixXd dwake_s_dx = DWakeSDFoilX(state, dfoil_s_dx);
    for (int w = 0; w < Nw - 1; ++w) {
        int k = (N - 1) + w;
        double ds_w = state.wake.s(w + 1) - state.wake.s(w);
        if (ds_w < 1e-14) continue;
        double coef = beta(k) * phi(k) / ds_w;
        vjp += coef * (dwake_s_dx.row(w + 1).transpose() -
                       dwake_s_dx.row(w).transpose());
    }

    // ---- A2. Wake-side dueinv chain (uewi) ---------------------------------
    // uewi[j] = V(foil.x, gam, Vinf, alpha, wake.x[j]) . wake.t[j]
    //
    // For j < Nw-1, wake.t[j] = V0[j] / |V0[j]| (set during BuildWake from V
    // at the same wake.x[j]). So uewi[j] = |V0[j]| and
    //   duewi[j]/dfoil_x = d|V|/dfoil_x = wake.t[j]^T * dV[j]/dfoil_x
    // because V[j]^T * (I - t t^T)/|V| = t^T - (t^T t) t^T = 0 by orthogonality.
    //
    // For j == Nw-1, wake.t[Nw-1] is V/|V| at x_pred_last (predictor point
    // from BuildWake's last iteration), NOT at wake.x[Nw-1]. The chain through
    // the predictor (which depends on wake.x[Nw-2] and v1 from the prior
    // station) is more involved. We add the dV term and a placeholder for
    // dwake.t (TODO: implement the proper predictor chain). Magnitude of
    // contribution is one station out of Nw, so adjoint impact is bounded.
    {
        MatrixXd dgam_dx_full = DGammaDFoilX(state, alpha);
        MatrixXd dwake_x_dx = DWakeXDFoilX(state);  // [2*Nw x 2N]
        for (int j = 0; j < Nw; ++j) {
            VelDeriv vd = InviscidVelocityWithDeriv(
                state, state.wake.x.col(j), dgam_dx_full, alpha);
            // d(wake.x[:, j])/d(foil.x) extracted from dwake_x_dx layout
            // (row d*Nw + j = d/dfoil.x of wake.x(d, j)).
            MatrixXd dwxj(2, twoN);
            dwxj.row(0) = dwake_x_dx.row(0 * Nw + j);
            dwxj.row(1) = dwake_x_dx.row(1 * Nw + j);
            MatrixXd dV = vd.V_xfoil + vd.V_xi * dwxj;  // [2 x 2N]
            Vector2d t_wake = state.wake.t.col(j);
            // duewi[j]/dfoil_x first term: t^T * dV
            vjp -= lambda_ue(N + j) *
                   (t_wake.transpose() * dV).transpose();
            // For j == Nw-1: V is NOT aligned with t_wake (t was set from
            // predictor). Add an approximate dwake.t chain via the local
            // norm Jacobian — this is exact when V and t are aligned but
            // approximate otherwise. TODO: full predictor chain.
            if (j == Nw - 1) {
                double Vnorm = vd.V.norm();
                if (Vnorm > 1e-12) {
                    Matrix2d nj = (Matrix2d::Identity() -
                                   t_wake * t_wake.transpose()) /
                                  Vnorm;
                    MatrixXd dwt = nj * dV;
                    vjp -= lambda_ue(N + j) *
                           (vd.V.transpose() * dwt).transpose();
                }
            }
        }
    }

    // ---- B0: precompute pieces shared by A2 and B for wake-row chains.
    // (Hoist out of A2 block so B can reuse.)
    // (No-op here — A2 already did this work; just note that dgam_dx_full and
    // dwake_x_dx are needed below too. We recompute in B for clarity, since
    // the A2 block has its own scope. TODO: refactor to share.)

    // ---- B. dueisigma chain (panel velocity derivative chain) -------------
    //
    // ue_sigma = [Bp; Dw] where Bp = -inv(AIC) * B is the airfoil source-to-ue
    // matrix and Dw = Cgam * Bp + Csig is the wake source-to-ue matrix (with
    // row 0 overwritten by Bp[N-1, :] for continuity).
    //
    // For the VJP we need -lambda_sg^T * dueisigma * phi.
    //
    // FOIL-SIDE PART (B1): -lambda_f^T * dBp * phi
    //   dBp = -A^{-1} (dA * Bp + dB), so:
    //   -lambda_f^T * dBp * phi = mu_f^T * (dA * Bp * phi + dB * phi)
    //   where mu_f = A^{-T} * lambda_f_ext  (lambda_f extended with 0).
    //
    //   - mu_f^T * dA * v (with v = [Bp*phi; 0]) → AICMatVecVJP
    //   - mu_f^T * dB * phi → panel-source-streamfunction derivative chain
    //
    // FIRST WAKE ROW (Dw[0, :] = Bp[N-1, :]) is folded into the foil-side
    // chain by adding lambda_w[0] into the effective lambda at foil node N-1.
    //
    // OTHER WAKE ROWS (Cgam*Bp + Csig with j > 0): TODO. Requires:
    //   - mu_w^T * dCgam * (Bp * phi) chain (wake velocity from gamma)
    //   - mu_w^T * Cgam * (dBp * phi) chain (combine with B1)
    //   - mu_w^T * dCsig * phi chain (wake source velocity)
    //   - chain through wake.x for both Csig and Cgam.
    {
        // Extract Bp from state.vsol.ue_sigma (top N rows).
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            // Build effective foil lambda including first wake row.
            // lambda_f_eff[i] = sgnue[i] * lambda_ue[i] (foil-side as before)
            // Plus lambda_w[0] * Bp[N-1, :] contribution: same as adding
            // lambda_w[0] * 1 to row N-1 of the foil lambda (since Dw[0]=Bp[N-1]).
            // Wait — Dw[0] uses lambda_w[0] directly (no sgnue scaling since
            // sgue[N+0]=1). So we add lambda_w[0]*Bp[N-1, :] which means we
            // include +lambda_w[0] at index N-1 of the effective lambda (with
            // sgue=1, no sgnue scaling here).
            VectorXd lambda_f_eff(N);
            for (int i = 0; i < N; ++i) {
                lambda_f_eff(i) = state.isol.sgnue(i) * lambda_ue(i);
            }
            // Add Dw[0] = Bp[N-1] contribution: lambda_w[0] * 1 * Bp[N-1]
            lambda_f_eff(N - 1) += lambda_ue(N);  // sgue at wake = 1

            // For j > 0 wake rows, Dw[j, :] = Cgam[j, :] * Bp + Csig[j, :].
            // The Cgam*Bp part folds into the foil-side dBp chain via
            // effective lambda for foil index k:
            //   lambda_f_eff[k] += sum_{j>0} lambda_ue[N+j] * Cgam[j, k]
            // (We need to rebuild Cgam here since CalcUeM doesn't save it.)
            MatrixXd Cgam(Nw, N);
            for (int j = 0; j < Nw; ++j) {
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d ti_w = state.wake.t.col(j);
                MatrixXd V_G;
                InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0,
                                 xi_w, true, &V_G);
                for (int k = 0; k < N; ++k) {
                    Cgam(j, k) = V_G(0, k) * ti_w(0) + V_G(1, k) * ti_w(1);
                }
            }
            // Cgam[0, :] contribution is replaced by Bp[N-1, :] in Dw[0], so
            // only j > 0 rows of Cgam enter the dBp chain.
            for (int j = 1; j < Nw; ++j) {
                double lw = lambda_ue(N + j);
                if (std::abs(lw) < 1e-20) continue;
                for (int k = 0; k < N; ++k) {
                    lambda_f_eff(k) += lw * Cgam(j, k);
                }
            }

            // mu_f = A^{-T} * [lambda_f_eff; 0]
            VectorXd lambda_f_ext = VectorXd::Zero(AIC_dim);
            lambda_f_ext.head(N) = lambda_f_eff;
            VectorXd mu_f = lu.solve(lambda_f_ext);

            // v = Bp * phi (foil columns of ue_sigma * phi gives v)
            // ue_sigma's top N rows = Bp, so v = ue_sigma.topRows(N) * phi
            VectorXd v_foil = ue_sig.topRows(N) * phi;  // [N]
            VectorXd v_ext = VectorXd::Zero(AIC_dim);
            v_ext.head(N) = v_foil;

            // B1a: mu_f^T * dA * v_ext  via AICMatVecVJP
            vjp += AICMatVecVJP(state, mu_f, v_ext);

            // B1b: mu_f^T * dB * phi
            //   B columns 0..N-2: foil constant-source panels (foil.x[k:k+2])
            //   B columns N-1..Npanel-1: wake linear-source panels (wake.x)
            //   For each eval i and panel k, accumulate:
            //     sens += mu_f[i] * phi[k] * dB[i, k]/dfoil_x
            //
            // Foil source panels first (no wake.x dependence — pure foil chain):
            for (int i = 0; i < N; ++i) {
                double mu_i = mu_f(i);
                if (std::abs(mu_i) < 1e-20) continue;
                Vector2d xi = state.foil.x.col(i);
                for (int k = 0; k < N - 1; ++k) {
                    double phi_k = phi(k);
                    if (std::abs(phi_k) < 1e-20) continue;
                    Matrix2d Xj;
                    Xj.col(0) = state.foil.x.col(k);
                    Xj.col(1) = state.foil.x.col(k + 1);
                    ConstSourceStreamDeriv d =
                        PanelConstSourceStreamDeriv(Xj, xi);
                    double w = mu_i * phi_k;
                    // xi derivatives (foil.x[i])
                    vjp(2 * i)     += w * d.val_xi(0);
                    vjp(2 * i + 1) += w * d.val_xi(1);
                    // Xj derivatives (foil.x[k], foil.x[k+1])
                    vjp(2 * k)         += w * d.val_Xj(0);
                    vjp(2 * k + 1)     += w * d.val_Xj(1);
                    vjp(2 * (k + 1))   += w * d.val_Xj(2);
                    vjp(2 * (k + 1) + 1) += w * d.val_Xj(3);
                }
            }
            // Wake source panels in B (cols N-1..Npanel-1):
            // CalcUeM builds these with half-panels at midpoints. Per panel j
            // (wake panel index 0..Nw-2), it computes left half [wake.x[j], Xm]
            // and right half [Xm, X_right] where Xm = 0.5*(wake.x[j]+wake.x[j+1])
            // and X_right is wake.x[j+1] (or ghost extension for j=Nw-2).
            //
            // Each call to PanelLinSourceStream returns (a, b) with derivatives
            // w.r.t. Xj endpoints (4-vector) and xi eval point (2-vector).
            // Chain through wake.x: dXj = (partial dwake.x for left endpoint,
            // 0.5*(dwake.x[j]+dwake.x[j+1]) for Xm endpoint, similar for right).
            MatrixXd dwake_x_dx = DWakeXDFoilX(state);  // [2*Nw x 2N]
            auto get_dwx = [&](int j) {
                MatrixXd dwj(2, twoN);
                dwj.row(0) = dwake_x_dx.row(0 * Nw + j);
                dwj.row(1) = dwake_x_dx.row(1 * Nw + j);
                return dwj;
            };
            for (int i = 0; i < N; ++i) {
                double mu_i = mu_f(i);
                if (std::abs(mu_i) < 1e-20) continue;
                Vector2d xi = state.foil.x.col(i);
                for (int j = 0; j < Nw - 1; ++j) {
                    Vector2d Xm = 0.5 * (state.wake.x.col(j) +
                                          state.wake.x.col(j + 1));
                    MatrixXd dXm = 0.5 * (get_dwx(j) + get_dwx(j + 1));

                    Vector2d X_right;
                    MatrixXd dXr(2, twoN);
                    if (j == Nw - 2) {
                        X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
                        dXr = 2.0 * get_dwx(j + 1) - dXm;
                    } else {
                        X_right = state.wake.x.col(j + 1);
                        dXr = get_dwx(j + 1);
                    }

                    // Left half panel
                    Matrix2d Xj_L;
                    Xj_L.col(0) = state.wake.x.col(j);
                    Xj_L.col(1) = Xm;
                    MatrixXd dXjL_0 = get_dwx(j);
                    MatrixXd dXjL_1 = dXm;
                    LinSourceStreamDeriv dL =
                        PanelLinSourceStreamDeriv(Xj_L, xi);
                    // Mapping a_left/b_left to B columns matches the +=
                    // structure in CalcUeM (lines 795-801).
                    double a_w_panel_curr_l, a_w_panel_prev_l, b_w_panel_curr_l;
                    if (j > 0) {
                        // B[i, N-1+j] += 0.5*a_left + b_left
                        // B[i, N-1+j-1] += 0.5*a_left
                        a_w_panel_curr_l = 0.5;  // coefficient of a in panel (N-1+j) col
                        a_w_panel_prev_l = 0.5;  // coefficient of a in panel (N-1+j-1) col
                        b_w_panel_curr_l = 1.0;  // coefficient of b in panel (N-1+j) col
                    } else {
                        // j=0: only b_left to panel (N-1+0)=N-1
                        a_w_panel_curr_l = 0.0;
                        a_w_panel_prev_l = 0.0;
                        b_w_panel_curr_l = 1.0;
                    }
                    int col_curr = N - 1 + j;
                    int col_prev = N - 1 + j - 1;
                    double phi_curr = phi(col_curr);
                    double phi_prev = (j > 0) ? phi(col_prev) : 0.0;
                    // Effective weight on the panel's (a, b) values:
                    double w_a = mu_i * (a_w_panel_curr_l * phi_curr +
                                          a_w_panel_prev_l * phi_prev);
                    double w_b = mu_i * (b_w_panel_curr_l * phi_curr);
                    // xi chain
                    vjp(2 * i)     += w_a * dL.a_xi(0) + w_b * dL.b_xi(0);
                    vjp(2 * i + 1) += w_a * dL.a_xi(1) + w_b * dL.b_xi(1);
                    // Xj chain — Xj cols 0,1 = first endpoint (wake.x[j]),
                    // cols 2,3 = second endpoint (Xm).
                    Vector4d dab = w_a * dL.a_Xj + w_b * dL.b_Xj;
                    // First endpoint = wake.x[j]: dXj rows = dXjL_0
                    vjp += dab(0) * dXjL_0.row(0).transpose() +
                           dab(1) * dXjL_0.row(1).transpose();
                    // Second endpoint = Xm: dXj rows = dXjL_1
                    vjp += dab(2) * dXjL_1.row(0).transpose() +
                           dab(3) * dXjL_1.row(1).transpose();

                    // Right half panel
                    Matrix2d Xj_R;
                    Xj_R.col(0) = Xm;
                    Xj_R.col(1) = X_right;
                    MatrixXd dXjR_0 = dXm;
                    MatrixXd dXjR_1 = dXr;
                    LinSourceStreamDeriv dR =
                        PanelLinSourceStreamDeriv(Xj_R, xi);
                    // B[i, N-1+j] += a_right + 0.5*b_right
                    // if (j<Nw-2): B[i, N-1+j+1] += 0.5*b_right
                    // else:       B[i, N-1+j] += 0.5*b_right
                    double a_w_panel_curr_r = 1.0;
                    double b_w_panel_curr_r;
                    double b_w_panel_next_r;
                    int col_next = N - 1 + j + 1;
                    double phi_next;
                    if (j < Nw - 2) {
                        b_w_panel_curr_r = 0.5;
                        b_w_panel_next_r = 0.5;
                        phi_next = phi(col_next);
                    } else {
                        b_w_panel_curr_r = 1.0;  // 0.5 + 0.5 = 1
                        b_w_panel_next_r = 0.0;
                        phi_next = 0.0;
                    }
                    w_a = mu_i * (a_w_panel_curr_r * phi_curr);
                    w_b = mu_i * (b_w_panel_curr_r * phi_curr +
                                  b_w_panel_next_r * phi_next);
                    vjp(2 * i)     += w_a * dR.a_xi(0) + w_b * dR.b_xi(0);
                    vjp(2 * i + 1) += w_a * dR.a_xi(1) + w_b * dR.b_xi(1);
                    dab = w_a * dR.a_Xj + w_b * dR.b_Xj;
                    // First endpoint = Xm
                    vjp += dab(0) * dXjR_0.row(0).transpose() +
                           dab(1) * dXjR_0.row(1).transpose();
                    // Second endpoint = X_right
                    vjp += dab(2) * dXjR_1.row(0).transpose() +
                           dab(3) * dXjR_1.row(1).transpose();
                }
            }
        }
    }
    // ---- B2a: dCgam direct chain for j > 0 wake rows -----------------------
    //
    // For each wake station j > 0:
    //   contribution = -lambda_w[j] * d(Cgam[j, :] * v)/dfoil_x
    //                = -lambda_w[j] * d(V_for_v . wake.t[j])/dfoil_x
    //                = -lambda_w[j] * (dV_for_v . wake.t[j] + V_for_v . dwake.t[j])
    // where V_for_v is the velocity at wake.x[j] computed with vortex
    // strengths = v (held fixed for this chain piece). dwake.t comes from the
    // REAL velocity chain (same as A2).
    //
    // For j = 0, Cgam[0, :] is replaced by Bp[N-1, :] in Dw[0]; handled in B1.
    {
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            VectorXd v_foil = ue_sig.topRows(N) * phi;  // [N]
            VectorXd v_padded(AIC_dim);
            v_padded.head(N) = v_foil;
            v_padded(AIC_dim - 1) = 0;
            MatrixXd dgam_zero = MatrixXd::Zero(AIC_dim, twoN);
            MatrixXd dgam_real = DGammaDFoilX(state, alpha);
            MatrixXd dwake_x_dx = DWakeXDFoilX(state);

            for (int j = 1; j < Nw; ++j) {
                double lw = lambda_ue(N + j);
                if (std::abs(lw) < 1e-20) continue;
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d t_wake = state.wake.t.col(j);

                // V_for_v at wake.x[j] (vortex strengths = v, dvortex = 0)
                VelDeriv vd_v = InviscidVelocityWithDerivExt(
                    state, xi_w, v_padded, dgam_zero, alpha);
                MatrixXd dwxj(2, twoN);
                dwxj.row(0) = dwake_x_dx.row(0 * Nw + j);
                dwxj.row(1) = dwake_x_dx.row(1 * Nw + j);
                MatrixXd dV_for_v = vd_v.V_xfoil + vd_v.V_xi * dwxj;

                vjp -= lw * (t_wake.transpose() * dV_for_v).transpose();

                // V_real chain for dwake.t (j < Nw-1: same as A2 simplification)
                VelDeriv vd_real = InviscidVelocityWithDeriv(
                    state, xi_w, dgam_real, alpha);
                MatrixXd dV_real = vd_real.V_xfoil + vd_real.V_xi * dwxj;
                double Vnorm = vd_real.V.norm();
                if (Vnorm > 1e-12) {
                    Matrix2d norm_jac = (Matrix2d::Identity() -
                                          t_wake * t_wake.transpose()) /
                                         Vnorm;
                    MatrixXd dwt = norm_jac * dV_real;
                    vjp -= lw * (vd_v.V.transpose() * dwt).transpose();
                }
            }
        }
    }

    // TODO: wake-side Csig chain (panel source velocity at wake stations,
    // foil-source panels + half-panel wake-source contributions). This is
    // the last missing chain piece.

    return vjp;
}

// Build DRueDFoilX [Nsys x 2N] by applying ComputeRueVJP with each standard
// basis vector as lambda. Slow (Nsys VJP calls) — for validation only.
MatrixXd DRueDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    MatrixXd dR_ue = MatrixXd::Zero(Nsys, 2 * N);
    VectorXd e = VectorXd::Zero(Nsys);
    for (int i = 0; i < Nsys; ++i) {
        e.setZero();
        e(i) = 1.0;
        dR_ue.row(i) = ComputeRueVJP(state, e, dfoil_s_dx).transpose();
    }
    return dR_ue;
}

}  // namespace cppfoil
