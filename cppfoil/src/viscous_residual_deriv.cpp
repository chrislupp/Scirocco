/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#include "cppfoil/core/viscous_residual_deriv.h"

#include <vector>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid_deriv.h"
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

    // ---- A2. Wake-side dueinv chain (uewi) — TODO ---------------------------
    // The wake inviscid edge velocity uewi[j] depends on foil.x (panel
    // influences from airfoil), on wake.x[:, j] (eval point), and on
    // wake.t[:, j] (tangent). Reverse-mode chain through these requires
    // panel-velocity derivatives at every wake station + DWakeXDFoilX.
    // For now this contribution is missing.

    // ---- B. dueisigma chain (panel velocity derivative chain) — TODO --------
    // ue_sigma's columns are source-panel-velocity influences at every
    // station, depending on panel positions and wake.x. Reverse-mode chain
    // requires propagating beta * phi back through the panel-influence
    // velocity derivatives. Largest missing piece for full accuracy.

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
