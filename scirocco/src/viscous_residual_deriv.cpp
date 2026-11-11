/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#include "scirocco/core/viscous_residual_deriv.h"

#include <vector>

#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/inviscid_deriv.h"
#include "scirocco/core/panel_influence_deriv.h"
#include "scirocco/core/viscous_geom_deriv.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

/**
 * @brief dR_BL/d(foil.x)
 *
 * dR_BL/d(foil.x)
 *   = R_x · DXiDFoilX                       (closed-form xi chain)
 *   + (wgap chain at non-initial wake stations from ResidualStation)
 *   + (hTE chain at the first wake station from WakeSys)
 *
 * R_x is the [3*Nsys x Nsys] sparse linearisation already assembled by
 * BuildGlobSys.
 *
 * wgap chain: ResidualStation modifies U via "U.row(1) -= Aux", so
 *   ∂R_station/∂Aux(0) = -∂R_station/∂U[1] at the first column (lin_U(:, 1)),
 *   ∂R_station/∂Aux(1) = -∂R_station/∂U[1] at the second column (lin_U(:, 5))
 * These derivatives live in R_U at rows 3*Is[2][ip]..3*Is[2][ip]+2,
 * columns 4*Is[2][ip-1]+1 and 4*Is[2][ip]+1 respectively. Aux at wake-station
 * local index `i` equals wgap(i), so the chain through DWgapDFoilX is:
 *   contribution to dR_BL row r = -R_U(r, 4*Is[2][ip-1]+1) * dwgap(ip-1)
 *                                 -R_U(r, 4*Is[2][ip]+1)   * dwgap(ip)
 * for each wake-station index ip = 1..Nw-1 (the first wake station is
 * handled by WakeSys, not ResidualStation).
 *
 * hTE chain at WakeSys: residual is R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE),
 * so ∂R(1)/∂hTE = -1, and ∂hTE/∂foil.x[k]_dim comes from GetTEInfoDeriv.
 * Only the second component of the first-wake-station residual triple is
 * affected.
 *
 * @param state        Fully assembled foil state (R_x and R_U populated).
 * @param dfoil_s_dx   Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 * @return dR_BL/d(foil.x) [3*Nsys x 2N].
 */
MatrixXd DRblDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    if (state.glob.R_x.rows() != 3 * Nsys || state.glob.R_x.cols() != Nsys) {
        throw std::runtime_error(
            "DRblDFoilX: state.glob.R_x not assembled — call BuildGlobSys "
            "first");
    }

    MatrixXd dxi_dx = DXiDFoilX(state, dfoil_s_dx);     /* [Nsys x 2N] */
    MatrixXd dwgap_dx = DWgapDFoilX(state, dfoil_s_dx); /* [Nw   x 2N] */

    /* ---- R_x · DXi ---- */
    MatrixXd dRbl = state.glob.R_x * dxi_dx; /* sparse * dense → dense */

    /* ---- wgap chain at wake stations ip = 1..Nw-1 ----
     *
     * ResidualStation modifies U via "U.row(1) -= Aux" AND uses Aux directly
     * in GetHw (Hw = Aux/U[0]). So Aux has TWO chain paths into R:
     *   1. U_mod[1] = U[1] - Aux → dU_mod[1]/dAux = -1, contributes -R_U[:,
     *   U1_col]
     *   2. Hw = Aux/U[0] → dHw/dAux = 1/U[0], contributes ∂R/∂Hw * 1/U[0]
     *
     * Hw appears in (Rmom, Rshape, Rlag) with coefficients (uelog, -uelog, 0).
     * Total ∂R/∂Aux = -R_U[:, U1_col] + [uelog, -uelog, 0] / U[0].
     *
     * For wake stations with cf=0, these two paths CANCEL for Rmom and Rshape
     * (the U_mod chain through H = (U[1]-Aux)/U[0] gives ±uelog/U[0] that
     * exactly cancels the Hw-direct path). Missing the Hw path produces a
     * spurious 0.5-0.9-magnitude chain term — caught at the Phase 2 test. */
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

        /* Compute uelog and U[0] at the two endpoints (Karman-Tsien-corrected
         * edge velocity uk via GetUk to match residual_station internals). */
        auto [uk_prev, unused_uk_U_prev] =
            GetUk(state.glob.U(3, node_prev), state.param);
        auto [uk_curr, unused_uk_U_curr] =
            GetUk(state.glob.U(3, node_curr), state.param);
        double uelog = std::log(uk_curr / uk_prev);
        double U0_prev = state.glob.U(0, node_prev);
        double U0_curr = state.glob.U(0, node_curr);

        /* Hw in residual is the AVERAGE: Hw_avg = 0.5*(Hw1 + Hw2). So the
         * per-endpoint Hw chain contributes 0.5*∂R/∂Hw_avg, not the full thing.
         * R coefficients for Hw_avg: +uelog (Rmom), -uelog (Rshape), 0 (Rlag).
         * Per-endpoint Hw_i: factor 0.5 from the average. */
        Vector3d Hw_coef(0.5 * uelog, -0.5 * uelog, 0.0);

        for (int r = 0; r < 3; ++r) {
            double cprev = RU.coeff(Ig + r, col_prev);
            double ccurr = RU.coeff(Ig + r, col_curr);
            /* Correct ∂R/∂Aux = (-R_U via U_mod chain) + (Hw chain via 1/U[0])
             */
            double dRdAux_prev = -cprev + Hw_coef(r) / U0_prev;
            double dRdAux_curr = -ccurr + Hw_coef(r) / U0_curr;
            /* Chain to foil.x: contribute (∂R/∂Aux) * dwgap_dx */
            if (dRdAux_prev != 0.0) {
                dRbl.row(Ig + r) += dRdAux_prev * dwgap_dx.row(ip - 1);
            }
            if (dRdAux_curr != 0.0) {
                dRbl.row(Ig + r) += dRdAux_curr * dwgap_dx.row(ip);
            }
        }
    }

    /* ---- hTE chain at the first wake station (WakeSys) ----
     * R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE), so dR(1)/dhTE = -1.
     * hTE depends on foil.x[k]_dim only for k ∈ {0, 1, N-2, N-1}. */
    int Ig_wake_init =
        3 * Iw[0] + 1; /* second residual of first wake station */
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

/**
 * @brief Reverse-mode VJP through R_ue chain.
 *
 * Returns lambda_ue^T * dR_ue/dfoil_x as a [2N] vector.
 *
 * R_ue[i] = U[3, i] - ueinv[i] - sgue[i] * (ue_sigma * phi)[i]
 *   phi   = sigma_m * mass,   mass[k] = U[1, k] * U[3, k]  (held fixed at U)
 *   sgue  = [sgnue (foil); ones (wake)]
 *   ueinv = uea for foil; for wake[0] continuity uea[N-1]; for wake[j>0] uewi.
 *
 * Chain pieces:
 *   A1. Foil-side dueinv via dgam (also the wake[0] continuity term).
 *   A2. Wake-side dueinv via uewi = V(wake.x[j]) . wake.t[j].
 *   B1. ue_sigma foil block (Bp = -A^{-1} B) for foil rows AND for the first
 *       wake row (Dw[0, :] = Bp[N-1, :]). Folds Cgam[j>0, :]*Bp into the foil
 *       chain via an effective foil-lambda.
 *   B2a. dCgam direct chain for wake rows j>0 (V_G*v . wake.t).
 *   B2b. dCsig foil-source chain.
 *   B2c. dCsig wake-source chain (non-self + interior self-influence log).
 *   C. sigma_m chain via foil.s/wake.s panel-length sensitivities.
 *
 * @param state        Fully assembled foil state.
 * @param lambda_ue    Adjoint vector [Nsys].
 * @param ctx          BuildRueVJPContext(state, dfoil_s_dx).
 * @return VJP result [2N].
 */
VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const RueVJPContext& ctx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    int twoN = 2 * N;
    double alpha = state.oper.alpha;

    if (lambda_ue.size() != Nsys) {
        throw std::runtime_error("ComputeRueVJP: lambda_ue size != Nsys");
    }
    if (ctx.dgam_dx.cols() != twoN || ctx.dwake_x_dx.rows() != 2 * Nw) {
        throw std::invalid_argument(
            "ComputeRueVJP: context was built for a different state");
    }
    const MatrixXd& dfoil_s_dx = ctx.dfoil_s_dx;

    VectorXd vjp = VectorXd::Zero(twoN);

    /* ---- A1. Foil-side dueinv chain via dgam
     * --------------------------------- lambda_g (N-vector) = sgnue[0..N) .*
     * lambda_ue[0..N) PLUS continuity contribution from the first wake row:
     *   get_ueinv sets ueinv[N+0] = uea[N-1] = sgnue[N-1] * gam[N-1]  (NOT
     *   uewi[0]). So d(ueinv[N+0]) chains through gam[N-1], same as if we added
     *   sgnue[N-1] * lambda_ue[N] to lambda_g[N-1].
     * Extended to N+1 with 0 in the streamfunction slot. */
    int AIC_dim = static_cast<int>(state.isol.AIC.rows()); /* N + 1 */
    VectorXd lambda_g_ext = VectorXd::Zero(AIC_dim);
    for (int i = 0; i < N; ++i) {
        lambda_g_ext(i) = state.isol.sgnue(i) * lambda_ue(i);
    }
    lambda_g_ext(N - 1) += state.isol.sgnue(N - 1) * lambda_ue(N);
    /* mu = AIC^{-T} * lambda_g_ext */
    const Eigen::PartialPivLU<MatrixXd>& lu = ctx.aic_t_lu;
    VectorXd mu = lu.solve(lambda_g_ext);
    /* Contribution: +AICResidualSensVJP(state, mu, alpha)
     *   (derivation: -lambda^T * dueinv_foil = -lambda_g^T * dgam,
     *    dgam = -(AIC)^{-1} * dR_AIC/dfoil_x, so -lambda_g^T * dgam =
     *    +(AIC^{-T} lambda_g)^T * dR_AIC = +AICResidualSensVJP.) */
    vjp += AICResidualSensVJP(state, mu, alpha);

    /* ---- C. dphi chain via dsigma_m ----------------------------------------
     * beta = ue_sigma^T * (sgue .* lambda_ue),  shape [Npanel = N+Nw-2] */
    VectorXd sgue_lambda(Nsys);
    for (int i = 0; i < N; ++i)
        sgue_lambda(i) = state.isol.sgnue(i) * lambda_ue(i);
    for (int i = N; i < Nsys; ++i) sgue_lambda(i) = lambda_ue(i);
    VectorXd beta =
        state.vsol.ue_sigma.transpose() * sgue_lambda; /* [Npanel] */

    /* phi = sigma_m * mass,  mass[k] = U[1, k] * U[3, k] */
    VectorXd mass(Nsys);
    for (int i = 0; i < Nsys; ++i)
        mass(i) = state.glob.U(1, i) * state.glob.U(3, i);
    VectorXd phi = state.vsol.sigma_m * mass; /* [Npanel] */

    /* Foil panels k = 0 .. N-2: ds_k = foil.s[k+1] - foil.s[k]
     * dsigma_m[k, i]/d(foil.x) = +sgnue[i]/ds_k^2 * d(ds_k)/d(foil.x)  for i =
     * k or k+1 sign-flipped dphi[k] = sum_i dsigma_m[k, i] * mass[i] =
     * -phi[k]/ds_k * d(ds_k) (Mass contribution multiplied through cancels via
     * -1/ds^2 * (sgnue[k+1]*mass[k+1] - sgnue[k]*mass[k]) = -phi[k]/ds.)
     *
     * Total contribution to VJP for panel k:
     *   -lambda_sg^T * ue_sigma * dphi[k]
     *     = -beta[k] * dphi[k]
     *     = -beta[k] * (-phi[k]/ds_k * d(ds_k)/dfoil_x)
     *     = (beta[k] * phi[k] / ds_k) * d(ds_k)/dfoil_x */
    for (int k = 0; k < N - 1; ++k) {
        double ds_k = state.foil.s(k + 1) - state.foil.s(k);
        if (ds_k < 1e-14) continue;
        double coef = beta(k) * phi(k) / ds_k;
        /* d(ds_k)/dfoil_x = dfoil_s_dx.row(k+1) - dfoil_s_dx.row(k) */
        vjp += coef * (dfoil_s_dx.row(k + 1).transpose() -
                       dfoil_s_dx.row(k).transpose());
    }
    /* Wake panels k = N-1 .. Npanel-1: ds_w = wake.s[w+1] - wake.s[w]
     * where w = k - (N-1) is the wake panel index.
     * Same structure: need d(ds_w)/d(foil.x) = dwake_s_dx.row(w+1) -
     * dwake_s_dx.row(w) */
    const MatrixXd& dwake_s_dx = ctx.dwake_s_dx;
    for (int w = 0; w < Nw - 1; ++w) {
        int k = (N - 1) + w;
        double ds_w = state.wake.s(w + 1) - state.wake.s(w);
        if (ds_w < 1e-14) continue;
        double coef = beta(k) * phi(k) / ds_w;
        vjp += coef * (dwake_s_dx.row(w + 1).transpose() -
                       dwake_s_dx.row(w).transpose());
    }

    /* ---- A2. Wake-side dueinv chain (uewi) ---------------------------------
     * uewi[j] = V(foil.x, gam, Vinf, alpha, wake.x[j]) . wake.t[j]
     *
     *   d(uewi)/dfoil_x = t_wake^T * dV + V^T * dwake.t
     *
     * For j < Nw-1, V is aligned with t_wake (BuildWake sets t = V/|V| at the
     * same point), so V^T * (I - t t^T) = 0 and the V . dwake.t term vanishes.
     * For j == Nw-1, wake.t is set from the predictor V/|V| at x_pred, NOT at
     * wake.x[Nw-1] — so V and t are NOT aligned and the V . dwake.t term is
     * nonzero. DWakeTDFoilX gives the correct dwake.t through the predictor
     * chain.
     *
     * SKIP j=0: ueinv[N+0] is overwritten by uea[N-1] in get_ueinv (continuity
     * BC), so its chain goes through gam[N-1] in A1, NOT through V at
     * wake.x[0]. */
    const MatrixXd& dgam_dx_full = ctx.dgam_dx;
    const MatrixXd& dwake_x_dx = ctx.dwake_x_dx; /* [2*Nw x 2N] */
    const MatrixXd& dwake_t_dx = ctx.dwake_t_dx; /* [2*Nw x 2N] */
    {
        for (int j = 1; j < Nw; ++j) {
            VelDeriv vd = InviscidVelocityWithDeriv(state, state.wake.x.col(j),
                                                    dgam_dx_full, alpha);
            MatrixXd dwxj(2, twoN);
            dwxj.row(0) = dwake_x_dx.row(0 * Nw + j);
            dwxj.row(1) = dwake_x_dx.row(1 * Nw + j);
            MatrixXd dV = vd.V_xfoil + vd.V_xi * dwxj; /* [2 x 2N] */
            Vector2d t_wake = state.wake.t.col(j);
            /* First term: t_wake^T * dV */
            vjp -= lambda_ue(N + j) * (t_wake.transpose() * dV).transpose();
            /* Second term: V^T * dwake.t (zero for j < Nw-1 since V || t). */
            MatrixXd dwt(2, twoN);
            dwt.row(0) = dwake_t_dx.row(0 * Nw + j);
            dwt.row(1) = dwake_t_dx.row(1 * Nw + j);
            vjp -= lambda_ue(N + j) * (vd.V.transpose() * dwt).transpose();
        }
    }

    /* ---- B. dueisigma chain (panel velocity derivative chain) -------------
     *
     * ue_sigma = [Bp; Dw] where Bp = -A^{-1} B (airfoil source-to-ue) and
     * Dw = Cgam * Bp + Csig (wake source-to-ue) with row 0 overwritten by
     * Bp[N-1, :] for continuity.
     *
     * B1: foil-side dBp = -A^{-1}(dA Bp + dB) chain via AICMatVecVJP + dB
     *     panel-streamfunction chain. Effective foil-lambda absorbs:
     *       (i) sgnue * lambda_ue[0:N] (foil rows),
     *       (ii) lambda_ue[N] at index N-1 (first wake row = Bp[N-1, :]),
     *       (iii) sum_{j>0} lambda_ue[N+j] * Cgam[j, k] (Cgam*Bp folded in).
     * B2a: dCgam direct chain for wake rows j>0 — d(V_G*v . wake.t).
     * B2b: dCsig foil-source-panel chain (V at wake.x[j_eval] from foil panel).
     * B2c: dCsig wake-source-panel chain (non-self + interior self log term).
     */
    {
        /* Extract Bp from state.vsol.ue_sigma (top N rows). */
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            /* Build effective foil lambda including first wake row.
             * lambda_f_eff[i] = sgnue[i] * lambda_ue[i] (foil-side as before)
             * Plus lambda_w[0] * Bp[N-1, :] contribution: same as adding
             * lambda_w[0] * 1 to row N-1 of the foil lambda (since
             * Dw[0]=Bp[N-1]). Wait — Dw[0] uses lambda_w[0] directly (no sgnue
             * scaling since sgue[N+0]=1). So we add lambda_w[0]*Bp[N-1, :]
             * which means we include +lambda_w[0] at index N-1 of the effective
             * lambda (with sgue=1, no sgnue scaling here). */
            VectorXd lambda_f_eff(N);
            for (int i = 0; i < N; ++i) {
                lambda_f_eff(i) = state.isol.sgnue(i) * lambda_ue(i);
            }
            /* Add Dw[0] = Bp[N-1] contribution: lambda_w[0] * 1 * Bp[N-1] */
            lambda_f_eff(N - 1) += lambda_ue(N); /* sgue at wake = 1 */

            /* For j > 0 wake rows, Dw[j, :] = Cgam[j, :] * Bp + Csig[j, :].
             * The Cgam*Bp part folds into the foil-side dBp chain via
             * effective lambda for foil index k:
             *   lambda_f_eff[k] += sum_{j>0} lambda_ue[N+j] * Cgam[j, k]
             * (We need to rebuild Cgam here since CalcUeM doesn't save it.) */
            MatrixXd Cgam(Nw, N);
            for (int j = 0; j < Nw; ++j) {
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d ti_w = state.wake.t.col(j);
                MatrixXd V_G;
                InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi_w,
                                 true, &V_G);
                for (int k = 0; k < N; ++k) {
                    Cgam(j, k) = V_G(0, k) * ti_w(0) + V_G(1, k) * ti_w(1);
                }
            }
            /* Cgam[0, :] contribution is replaced by Bp[N-1, :] in Dw[0], so
             * only j > 0 rows of Cgam enter the dBp chain. */
            for (int j = 1; j < Nw; ++j) {
                double lw = lambda_ue(N + j);
                if (std::abs(lw) < 1e-20) continue;
                for (int k = 0; k < N; ++k) {
                    lambda_f_eff(k) += lw * Cgam(j, k);
                }
            }

            /* mu_f = A^{-T} * [lambda_f_eff; 0] */
            VectorXd lambda_f_ext = VectorXd::Zero(AIC_dim);
            lambda_f_ext.head(N) = lambda_f_eff;
            VectorXd mu_f = lu.solve(lambda_f_ext);

            /* v = Bp * phi (foil columns of ue_sigma * phi gives v)
             * ue_sigma's top N rows = Bp, so v = ue_sigma.topRows(N) * phi */
            VectorXd v_foil = ue_sig.topRows(N) * phi; /* [N] */
            VectorXd v_ext = VectorXd::Zero(AIC_dim);
            v_ext.head(N) = v_foil;

            /* B1a: mu_f^T * dA * v_ext  via AICMatVecVJP */
            vjp += AICMatVecVJP(state, mu_f, v_ext);

            /* B1b: mu_f^T * dB * phi
             *   B columns 0..N-2: foil constant-source panels (foil.x[k:k+2])
             *   B columns N-1..Npanel-1: wake linear-source panels (wake.x)
             *   For each eval i and panel k, accumulate:
             *     sens += mu_f[i] * phi[k] * dB[i, k]/dfoil_x
             *
             * Foil source panels first (no wake.x dependence — pure foil
             * chain): */
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
                    /* xi derivatives (foil.x[i]) */
                    vjp(2 * i) += w * d.val_xi(0);
                    vjp(2 * i + 1) += w * d.val_xi(1);
                    /* Xj derivatives (foil.x[k], foil.x[k+1]) */
                    vjp(2 * k) += w * d.val_Xj(0);
                    vjp(2 * k + 1) += w * d.val_Xj(1);
                    vjp(2 * (k + 1)) += w * d.val_Xj(2);
                    vjp(2 * (k + 1) + 1) += w * d.val_Xj(3);
                }
            }
            /* Wake source panels in B (cols N-1..Npanel-1):
             * CalcUeM builds these with half-panels at midpoints. Per panel j
             * (wake panel index 0..Nw-2), it computes left half [wake.x[j], Xm]
             * and right half [Xm, X_right] where Xm =
             * 0.5*(wake.x[j]+wake.x[j+1]) and X_right is wake.x[j+1] (or ghost
             * extension for j=Nw-2).
             *
             * Each call to PanelLinSourceStream returns (a, b) with derivatives
             * w.r.t. Xj endpoints (4-vector) and xi eval point (2-vector).
             * Chain through wake.x: dXj = (partial dwake.x for left endpoint,
             * 0.5*(dwake.x[j]+dwake.x[j+1]) for Xm endpoint, similar for
             * right). Uses dwake_x_dx from the context (see A2). */
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
                    Vector2d Xm =
                        0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
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

                    /* Left half panel */
                    Matrix2d Xj_L;
                    Xj_L.col(0) = state.wake.x.col(j);
                    Xj_L.col(1) = Xm;
                    MatrixXd dXjL_0 = get_dwx(j);
                    MatrixXd dXjL_1 = dXm;
                    LinSourceStreamDeriv dL =
                        PanelLinSourceStreamDeriv(Xj_L, xi);
                    /* Mapping a_left/b_left to B columns matches the +=
                     * structure in CalcUeM (lines 795-801). */
                    double a_w_panel_curr_l, a_w_panel_prev_l, b_w_panel_curr_l;
                    if (j > 0) {
                        /* B[i, N-1+j] += 0.5*a_left + b_left
                         * B[i, N-1+j-1] += 0.5*a_left */
                        a_w_panel_curr_l =
                            0.5; /* coefficient of a in panel (N-1+j) col */
                        a_w_panel_prev_l =
                            0.5; /* coefficient of a in panel (N-1+j-1) col */
                        b_w_panel_curr_l =
                            1.0; /* coefficient of b in panel (N-1+j) col */
                    } else {
                        /* j=0: only b_left to panel (N-1+0)=N-1 */
                        a_w_panel_curr_l = 0.0;
                        a_w_panel_prev_l = 0.0;
                        b_w_panel_curr_l = 1.0;
                    }
                    int col_curr = N - 1 + j;
                    int col_prev = N - 1 + j - 1;
                    double phi_curr = phi(col_curr);
                    double phi_prev = (j > 0) ? phi(col_prev) : 0.0;
                    /* Effective weight on the panel's (a, b) values: */
                    double w_a = mu_i * (a_w_panel_curr_l * phi_curr +
                                         a_w_panel_prev_l * phi_prev);
                    double w_b = mu_i * (b_w_panel_curr_l * phi_curr);
                    /* xi chain */
                    vjp(2 * i) += w_a * dL.a_xi(0) + w_b * dL.b_xi(0);
                    vjp(2 * i + 1) += w_a * dL.a_xi(1) + w_b * dL.b_xi(1);
                    /* Xj chain — Xj cols 0,1 = first endpoint (wake.x[j]),
                     * cols 2,3 = second endpoint (Xm). */
                    Vector4d dab = w_a * dL.a_Xj + w_b * dL.b_Xj;
                    /* First endpoint = wake.x[j]: dXj rows = dXjL_0 */
                    vjp += dab(0) * dXjL_0.row(0).transpose() +
                           dab(1) * dXjL_0.row(1).transpose();
                    /* Second endpoint = Xm: dXj rows = dXjL_1 */
                    vjp += dab(2) * dXjL_1.row(0).transpose() +
                           dab(3) * dXjL_1.row(1).transpose();

                    /* Right half panel */
                    Matrix2d Xj_R;
                    Xj_R.col(0) = Xm;
                    Xj_R.col(1) = X_right;
                    MatrixXd dXjR_0 = dXm;
                    MatrixXd dXjR_1 = dXr;
                    LinSourceStreamDeriv dR =
                        PanelLinSourceStreamDeriv(Xj_R, xi);
                    /* B[i, N-1+j] += a_right + 0.5*b_right
                     * if (j<Nw-2): B[i, N-1+j+1] += 0.5*b_right
                     * else:       B[i, N-1+j] += 0.5*b_right */
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
                        b_w_panel_curr_r = 1.0; /* 0.5 + 0.5 = 1 */
                        b_w_panel_next_r = 0.0;
                        phi_next = 0.0;
                    }
                    w_a = mu_i * (a_w_panel_curr_r * phi_curr);
                    w_b = mu_i * (b_w_panel_curr_r * phi_curr +
                                  b_w_panel_next_r * phi_next);
                    vjp(2 * i) += w_a * dR.a_xi(0) + w_b * dR.b_xi(0);
                    vjp(2 * i + 1) += w_a * dR.a_xi(1) + w_b * dR.b_xi(1);
                    dab = w_a * dR.a_Xj + w_b * dR.b_Xj;
                    /* First endpoint = Xm */
                    vjp += dab(0) * dXjR_0.row(0).transpose() +
                           dab(1) * dXjR_0.row(1).transpose();
                    /* Second endpoint = X_right */
                    vjp += dab(2) * dXjR_1.row(0).transpose() +
                           dab(3) * dXjR_1.row(1).transpose();
                }
            }
        }
    }
    /* ---- B2a: dCgam direct chain for j > 0 wake rows -----------------------
     *
     * For each wake station j > 0:
     *   contribution = -lambda_w[j] * d(Cgam[j, :] * v)/dfoil_x
     *                = -lambda_w[j] * d(V_for_v . wake.t[j])/dfoil_x
     *                = -lambda_w[j] * (dV_for_v . wake.t[j] + V_for_v .
     *                dwake.t[j])
     * where V_for_v is the velocity at wake.x[j] computed with vortex
     * strengths = v (held fixed for this chain piece). dwake.t comes from the
     * REAL velocity chain (same as A2).
     *
     * For j = 0, Cgam[0, :] is replaced by Bp[N-1, :] in Dw[0]; handled in B1.
     */
    {
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            VectorXd v_foil = ue_sig.topRows(N) * phi; /* [N] */
            VectorXd v_padded(AIC_dim);
            v_padded.head(N) = v_foil;
            v_padded(AIC_dim - 1) = 0;
            MatrixXd dgam_zero = MatrixXd::Zero(AIC_dim, twoN);

            for (int j = 1; j < Nw; ++j) {
                double lw = lambda_ue(N + j);
                if (std::abs(lw) < 1e-20) continue;
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d t_wake = state.wake.t.col(j);

                /* V_for_v at wake.x[j] (vortex strengths = v, dvortex = 0).
                 * NOTE: vd_v.V from InviscidVelocityWithDerivExt INCLUDES the
                 * freestream Vinf*(cos α, sin α). Cgam[j, :] * v =
                 * (V_G * v).t_w which is gamma-induced only — no freestream.
                 * Subtract Vinf before using vd_v.V in the dwake.t chain. */
                VelDeriv vd_v = InviscidVelocityWithDerivExt(
                    state, xi_w, v_padded, dgam_zero, alpha);
                Vector2d Vinf_vec(state.oper.Vinf * Cosd(alpha),
                                  state.oper.Vinf * Sind(alpha));
                Vector2d V_G_v = vd_v.V - Vinf_vec;
                MatrixXd dwxj(2, twoN);
                dwxj.row(0) = dwake_x_dx.row(0 * Nw + j);
                dwxj.row(1) = dwake_x_dx.row(1 * Nw + j);
                MatrixXd dV_for_v = vd_v.V_xfoil + vd_v.V_xi * dwxj;

                vjp -= lw * (t_wake.transpose() * dV_for_v).transpose();

                /* dwake.t[j] chain via DWakeTDFoilX (exact predictor chain). */
                MatrixXd dwt(2, twoN);
                dwt.row(0) = dwake_t_dx.row(0 * Nw + j);
                dwt.row(1) = dwake_t_dx.row(1 * Nw + j);
                vjp -= lw * (V_G_v.transpose() * dwt).transpose();
            }
        }
    }

    /* ---- B2b: dCsig direct chain for j > 0 wake rows (foil-source part) ---
     *
     * Csig[j, k] = V_at_wake_j(from panel k) . wake.t[j]
     * For k = 0..N-2 (foil source panels), PanelConstSourceVelocity gives V
     * from a foil panel evaluated at wake.x[j]. Chain through Xj=foil.x[k:k+2]
     * (direct), xi=wake.x[j] (via dwake.x[j]), and wake.t[j] (via dV_real
     * chain).
     *
     * Skip cases from CalcUeM: (j == 0) skips foil panels k=0 and k=N-2.
     * For j > 0, all foil source panels contribute. */
    {
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            for (int j = 1; j < Nw; ++j) {
                double lw = lambda_ue(N + j);
                if (std::abs(lw) < 1e-20) continue;
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d t_w = state.wake.t.col(j);

                MatrixXd dwxj(2, twoN);
                dwxj.row(0) = dwake_x_dx.row(0 * Nw + j);
                dwxj.row(1) = dwake_x_dx.row(1 * Nw + j);

                /* dwake.t[j] via DWakeTDFoilX (predictor-aware). */
                MatrixXd dwt(2, twoN);
                dwt.row(0) = dwake_t_dx.row(0 * Nw + j);
                dwt.row(1) = dwake_t_dx.row(1 * Nw + j);

                /* Iterate over foil source panels (k = 0..N-2) */
                for (int k = 0; k < N - 1; ++k) {
                    double phi_k = phi(k);
                    if (std::abs(phi_k) < 1e-20) continue;
                    double w = -lw * phi_k; /* weight for d(Csig[j,k]) */

                    Matrix2d Xj_foil;
                    Xj_foil.col(0) = state.foil.x.col(k);
                    Xj_foil.col(1) = state.foil.x.col(k + 1);
                    ConstSourceVelocityDeriv csv =
                        PanelConstSourceVelocityDeriv(Xj_foil, xi_w);
                    Vector2d V_panel = csv.val;

                    /* Chain 1: Xj derivatives — w * (val_Xj.col(c) . t)
                     * contribute to (foil.x[k], foil.x[k+1]) positions */
                    for (int c = 0; c < 2; ++c) {
                        double d_proj = csv.val_Xj(0, c) * t_w(0) +
                                        csv.val_Xj(1, c) * t_w(1);
                        vjp(2 * k + c) += w * d_proj;
                    }
                    for (int c = 2; c < 4; ++c) {
                        double d_proj = csv.val_Xj(0, c) * t_w(0) +
                                        csv.val_Xj(1, c) * t_w(1);
                        vjp(2 * (k + 1) + (c - 2)) += w * d_proj;
                    }

                    /* Chain 2: xi derivatives chained through dwake.x[j]
                     * d(V.t)/dxi [1x2] = t^T * val_xi */
                    Vector2d dprojdxi(
                        csv.val_xi(0, 0) * t_w(0) + csv.val_xi(1, 0) * t_w(1),
                        csv.val_xi(0, 1) * t_w(0) + csv.val_xi(1, 1) * t_w(1));
                    vjp += w * (dprojdxi.transpose() * dwxj).transpose();

                    /* Chain 3: V_panel . dwake.t[j] = V_panel^T * dwt */
                    vjp += w * (V_panel.transpose() * dwt).transpose();
                }
            }
        }
    }

    /* ---- B2c: dCsig wake-source-panel chain --------------------------------
     *
     * For each wake eval j_eval > 0 and each wake source panel index
     * wp (0..Nw-1, corresponds to Csig column N-1+wp-1 / N-1+wp / N-1+wp+1
     * depending on the half-panel decomposition; matches CalcUeM
     * mfoil.py:1533-1569).
     *
     * Cases:
     *   self-influence interior (1 ≤ j_eval = wp ≤ Nw-2): log-formula
     *     contribution implemented below (the standalone block right after
     *     this comment, guarded by `j_eval >= 1 && j_eval <= Nw - 2`).
     *   self-influence wp == 0 (TE three-panel system): unreachable — the
     *     outer loop starts at j_eval = 1, and Dw[0,:] is overridden by
     *     Bp[N-1,:] for continuity (handled in B1; see lambda_f_eff bump at
     *     line ~316), so Csig[0,:] never enters the VJP.
     *   self-influence wp == Nw-1 (ghost-extended last panel): identically
     *     zero (mfoil.py:1550-1551, "no self effect of last pan").
     *   non-self wp == 0 (first wake panel): half-panel decomp on right side
     *     only.
     *   non-self wp == Nw-1 (last wake panel): const-source full panel.
     *   non-self interior: left + right half panels.
     *
     * Chain through wake.x (panel coords + midpoint coords) and wake.t (dwt).
     */
    {
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            auto get_dwx2 = [&](int j) -> MatrixXd {
                MatrixXd dwj(2, twoN);
                dwj.row(0) = dwake_x_dx.row(0 * Nw + j);
                dwj.row(1) = dwake_x_dx.row(1 * Nw + j);
                return dwj;
            };

            for (int j_eval = 1; j_eval < Nw; ++j_eval) {
                double lw = lambda_ue(N + j_eval);
                if (std::abs(lw) < 1e-20) continue;
                Vector2d xi_w = state.wake.x.col(j_eval);
                Vector2d t_w = state.wake.t.col(j_eval);

                MatrixXd dwxj = get_dwx2(j_eval);
                /* dwake.t[j_eval] via DWakeTDFoilX (predictor-aware) */
                MatrixXd dwt(2, twoN);
                dwt.row(0) = dwake_t_dx.row(0 * Nw + j_eval);
                dwt.row(1) = dwake_t_dx.row(1 * Nw + j_eval);

                /* Self-influence: interior wp == j_eval case (skipped for
                 * endpoints). */
                if (j_eval >= 1 && j_eval <= Nw - 2) {
                    int wp = j_eval;
                    Vector2d Xl =
                        0.5 * (state.wake.x.col(wp - 1) + state.wake.x.col(wp));
                    Vector2d Xm = state.wake.x.col(wp);
                    Vector2d Xr =
                        0.5 * (state.wake.x.col(wp) + state.wake.x.col(wp + 1));
                    MatrixXd dXl = 0.5 * (get_dwx2(wp - 1) + get_dwx2(wp));
                    MatrixXd dXm = get_dwx2(wp);
                    MatrixXd dXr = 0.5 * (get_dwx2(wp) + get_dwx2(wp + 1));
                    Vector2d e1 = Xm - Xl;
                    Vector2d e2 = Xr - Xm;
                    double d1 = e1.norm();
                    double d2 = e2.norm();
                    if (d1 > 1e-14 && d2 > 1e-14) {
                        Vector2d e1u = e1 / d1;
                        Vector2d e2u = e2 / d2;
                        MatrixXd de1 = dXm - dXl;
                        MatrixXd de2 = dXr - dXm;
                        VectorXd dd1_dx = (e1u.transpose() * de1).transpose();
                        VectorXd dd2_dx = (e2u.transpose() * de2).transpose();
                        VectorXd daa_dx =
                            (0.25 / kPi) * (dd1_dx / d1 - dd2_dx / d2);
                        int col_prev = N - 1 + wp - 1;
                        int col_curr = N - 1 + wp;
                        double w = -lw * (phi(col_prev) + phi(col_curr));
                        vjp += w * daa_dx;
                    }
                }

                /* Iterate over each wake source panel wp (0..Nw-1 wake nodes)
                 * Skip self-influence (handled above for interior). */
                for (int wp = 0; wp < Nw; ++wp) {
                    if (wp == j_eval) continue;

                    /* CalcUeM's wake source assembly: see lines 877-908. */
                    if (wp == 0) {
                        /* First wake "panel": Xj = [wake.x[0],
                         * 0.5*(wake.x[0]+wake.x[1])] Csig contributions:
                         * col(N-1) += b.dot(ti); col(0) += a.dot(ti); col(N-2)
                         * += a.dot(ti) */
                        Vector2d Xm0 =
                            0.5 * (state.wake.x.col(0) + state.wake.x.col(1));
                        MatrixXd dXm0 = 0.5 * (get_dwx2(0) + get_dwx2(1));
                        Matrix2d Xj;
                        Xj.col(0) = state.wake.x.col(0);
                        Xj.col(1) = Xm0;
                        MatrixXd dXj0 = get_dwx2(0);
                        const MatrixXd& dXj1 = dXm0;
                        LinSourceVelocityDeriv lvd =
                            PanelLinSourceVelocityDeriv(Xj, xi_w);
                        /* Csig coefficients on (a, b) for the 3 affected cols:
                         * col(N-1): w_b * (b.dot(ti)) → weight phi[N-1]
                         * col(0): w_a * (a.dot(ti)) → weight phi[0]
                         * col(N-2): w_a * (a.dot(ti)) → weight phi[N-2] */
                        double phi_N_1 = phi(N - 1);
                        double phi_0 = phi(0);
                        double phi_N_2 = phi(N - 2);
                        double w_a =
                            -lw *
                            (phi_0 +
                             phi_N_2); /* a contributes to cols 0 and N-2 */
                        double w_b = -lw * phi_N_1;
                        /* Project a, b onto t for VJP weights */
                        auto add_velocity_chain =
                            [&](double wa, double wb,
                                const LinSourceVelocityDeriv& d,
                                const MatrixXd& dXj0_, const MatrixXd& dXj1_) {
                                /* Chain 1: Xj derivatives */
                                for (int c = 0; c < 2; ++c) {
                                    double pa = wa * (d.a_Xj(0, c) * t_w(0) +
                                                      d.a_Xj(1, c) * t_w(1));
                                    double pb = wb * (d.b_Xj(0, c) * t_w(0) +
                                                      d.b_Xj(1, c) * t_w(1));
                                    /* Map to first endpoint (dXj0) */
                                    vjp += pa * dXj0_.row(c).transpose();
                                    vjp += pb * dXj0_.row(c).transpose();
                                }
                                for (int c = 2; c < 4; ++c) {
                                    double pa = wa * (d.a_Xj(0, c) * t_w(0) +
                                                      d.a_Xj(1, c) * t_w(1));
                                    double pb = wb * (d.b_Xj(0, c) * t_w(0) +
                                                      d.b_Xj(1, c) * t_w(1));
                                    vjp += pa * dXj1_.row(c - 2).transpose();
                                    vjp += pb * dXj1_.row(c - 2).transpose();
                                }
                                /* Chain 2: xi derivatives via dwake.x[j_eval]
                                 */
                                Vector2d dpa_dxi(wa * (d.a_xi(0, 0) * t_w(0) +
                                                       d.a_xi(1, 0) * t_w(1)),
                                                 wa * (d.a_xi(0, 1) * t_w(0) +
                                                       d.a_xi(1, 1) * t_w(1)));
                                Vector2d dpb_dxi(wb * (d.b_xi(0, 0) * t_w(0) +
                                                       d.b_xi(1, 0) * t_w(1)),
                                                 wb * (d.b_xi(0, 1) * t_w(0) +
                                                       d.b_xi(1, 1) * t_w(1)));
                                vjp += (dpa_dxi.transpose() * dwxj).transpose();
                                vjp += (dpb_dxi.transpose() * dwxj).transpose();
                            };
                        add_velocity_chain(w_a, w_b, lvd, dXj0, dXj1);
                        /* V_panel . dwake.t[j_eval] chain.
                         * V contributions to Csig cols:
                         *   col(N-1): b   col(0): a   col(N-2): a
                         * V_total = phi[N-1]*b + (phi[0]+phi[N-2])*a */
                        {
                            Vector2d V_total = phi(N - 1) * lvd.b +
                                               (phi(0) + phi(N - 2)) * lvd.a;
                            vjp +=
                                -lw * (V_total.transpose() * dwt).transpose();
                        }
                    } else if (wp == Nw - 1) {
                        /* Last wake panel: full constant-source from [Xl, Xr]
                         * Csig col (N + Nw - 3 = Npanel - 1) gets v.dot(ti) */
                        Vector2d Xl = 0.5 * (state.wake.x.col(wp - 1) +
                                             state.wake.x.col(wp));
                        Vector2d Xm = state.wake.x.col(wp);
                        Vector2d Xr = 2.0 * Xm - Xl;
                        MatrixXd dXl = 0.5 * (get_dwx2(wp - 1) + get_dwx2(wp));
                        MatrixXd dXm = get_dwx2(wp);
                        MatrixXd dXr = 2.0 * dXm - dXl;
                        Matrix2d Xj;
                        Xj.col(0) = Xl;
                        Xj.col(1) = Xr;
                        ConstSourceVelocityDeriv csv =
                            PanelConstSourceVelocityDeriv(Xj, xi_w);
                        int col_last = Npanel - 1;
                        double w = -lw * phi(col_last);
                        /* Chain 1: Xj derivatives (Xl and Xr; map their
                         * dXj0/dXj1 via dXl/dXr) */
                        for (int c = 0; c < 2; ++c) {
                            double d_proj = csv.val_Xj(0, c) * t_w(0) +
                                            csv.val_Xj(1, c) * t_w(1);
                            vjp += w * d_proj * dXl.row(c).transpose();
                        }
                        for (int c = 2; c < 4; ++c) {
                            double d_proj = csv.val_Xj(0, c) * t_w(0) +
                                            csv.val_Xj(1, c) * t_w(1);
                            vjp += w * d_proj * dXr.row(c - 2).transpose();
                        }
                        /* Chain 2: xi via dwxj */
                        Vector2d d_proj_xi(csv.val_xi(0, 0) * t_w(0) +
                                               csv.val_xi(1, 0) * t_w(1),
                                           csv.val_xi(0, 1) * t_w(0) +
                                               csv.val_xi(1, 1) * t_w(1));
                        vjp += w * (d_proj_xi.transpose() * dwxj).transpose();
                        /* Chain 3: V_panel . dwake.t[j_eval] */
                        vjp += w * (csv.val.transpose() * dwt).transpose();
                    } else {
                        /* Interior wake panel wp: left half [Xl, Xm] + right
                         * half [Xm, Xr] */
                        Vector2d Xl = 0.5 * (state.wake.x.col(wp - 1) +
                                             state.wake.x.col(wp));
                        Vector2d Xm = state.wake.x.col(wp);
                        Vector2d Xr = 0.5 * (state.wake.x.col(wp) +
                                             state.wake.x.col(wp + 1));
                        MatrixXd dXl = 0.5 * (get_dwx2(wp - 1) + get_dwx2(wp));
                        MatrixXd dXm = get_dwx2(wp);
                        MatrixXd dXr = 0.5 * (get_dwx2(wp) + get_dwx2(wp + 1));

                        /* Csig accumulation per CalcUeM (lines 893-908):
                         *   Csig(j_eval, N-1+wp-1) += a1.dot(ti) + 0.5 *
                         *   b1.dot(ti) Csig(j_eval, N-1+wp)   += 0.5 *
                         *   a2.dot(ti) + b2.dot(ti)
                         * So phi-weighted: w_a1 = phi[N-1+wp-1], w_b1 =
                         * 0.5*phi[N-1+wp-1]
                         *                  w_a2 = 0.5*phi[N-1+wp],  w_b2 =
                         *                  phi[N-1+wp] */
                        double phi_prev = phi(N - 1 + wp - 1);
                        double phi_curr = phi(N - 1 + wp);
                        double w_a1 = -lw * phi_prev;
                        double w_b1 = -lw * 0.5 * phi_prev;
                        double w_a2 = -lw * 0.5 * phi_curr;
                        double w_b2 = -lw * phi_curr;

                        /* Left half [Xl, Xm] */
                        Matrix2d Xj_L;
                        Xj_L.col(0) = Xl;
                        Xj_L.col(1) = Xm;
                        LinSourceVelocityDeriv dL =
                            PanelLinSourceVelocityDeriv(Xj_L, xi_w);
                        /* Xj endpoints: 0,1 = Xl (via dXl), 2,3 = Xm (via dXm)
                         */
                        for (int c = 0; c < 2; ++c) {
                            double pa = w_a1 * (dL.a_Xj(0, c) * t_w(0) +
                                                dL.a_Xj(1, c) * t_w(1));
                            double pb = w_b1 * (dL.b_Xj(0, c) * t_w(0) +
                                                dL.b_Xj(1, c) * t_w(1));
                            vjp += (pa + pb) * dXl.row(c).transpose();
                        }
                        for (int c = 2; c < 4; ++c) {
                            double pa = w_a1 * (dL.a_Xj(0, c) * t_w(0) +
                                                dL.a_Xj(1, c) * t_w(1));
                            double pb = w_b1 * (dL.b_Xj(0, c) * t_w(0) +
                                                dL.b_Xj(1, c) * t_w(1));
                            vjp += (pa + pb) * dXm.row(c - 2).transpose();
                        }
                        Vector2d dxi_L(w_a1 * (dL.a_xi(0, 0) * t_w(0) +
                                               dL.a_xi(1, 0) * t_w(1)) +
                                           w_b1 * (dL.b_xi(0, 0) * t_w(0) +
                                                   dL.b_xi(1, 0) * t_w(1)),
                                       w_a1 * (dL.a_xi(0, 1) * t_w(0) +
                                               dL.a_xi(1, 1) * t_w(1)) +
                                           w_b1 * (dL.b_xi(0, 1) * t_w(0) +
                                                   dL.b_xi(1, 1) * t_w(1)));
                        vjp += (dxi_L.transpose() * dwxj).transpose();

                        /* Right half [Xm, Xr] */
                        Matrix2d Xj_R;
                        Xj_R.col(0) = Xm;
                        Xj_R.col(1) = Xr;
                        LinSourceVelocityDeriv dR =
                            PanelLinSourceVelocityDeriv(Xj_R, xi_w);
                        for (int c = 0; c < 2; ++c) {
                            double pa = w_a2 * (dR.a_Xj(0, c) * t_w(0) +
                                                dR.a_Xj(1, c) * t_w(1));
                            double pb = w_b2 * (dR.b_Xj(0, c) * t_w(0) +
                                                dR.b_Xj(1, c) * t_w(1));
                            vjp += (pa + pb) * dXm.row(c).transpose();
                        }
                        for (int c = 2; c < 4; ++c) {
                            double pa = w_a2 * (dR.a_Xj(0, c) * t_w(0) +
                                                dR.a_Xj(1, c) * t_w(1));
                            double pb = w_b2 * (dR.b_Xj(0, c) * t_w(0) +
                                                dR.b_Xj(1, c) * t_w(1));
                            vjp += (pa + pb) * dXr.row(c - 2).transpose();
                        }
                        Vector2d dxi_R(w_a2 * (dR.a_xi(0, 0) * t_w(0) +
                                               dR.a_xi(1, 0) * t_w(1)) +
                                           w_b2 * (dR.b_xi(0, 0) * t_w(0) +
                                                   dR.b_xi(1, 0) * t_w(1)),
                                       w_a2 * (dR.a_xi(0, 1) * t_w(0) +
                                               dR.a_xi(1, 1) * t_w(1)) +
                                           w_b2 * (dR.b_xi(0, 1) * t_w(0) +
                                                   dR.b_xi(1, 1) * t_w(1)));
                        vjp += (dxi_R.transpose() * dwxj).transpose();
                        /* V_panel . dwake.t[j_eval] chain (left + right half).
                         */
                        {
                            Vector2d V_left =
                                phi_prev * dL.a + 0.5 * phi_prev * dL.b;
                            Vector2d V_right =
                                0.5 * phi_curr * dR.a + phi_curr * dR.b;
                            Vector2d V_total = V_left + V_right;
                            vjp +=
                                -lw * (V_total.transpose() * dwt).transpose();
                        }
                    }
                }
            }
        }
    }
    return vjp;
}

RueVJPContext BuildRueVJPContext(const FoilState& state,
                                 const MatrixXd& dfoil_s_dx) {
    RueVJPContext ctx;
    ctx.dfoil_s_dx = dfoil_s_dx;
    ctx.dgam_dx = DGammaDFoilX(state, state.oper.alpha);
    /* The wake Jacobians need d(gamma)/d(foil.x) at the tracing alpha. That
     * is the same matrix unless the wake was traced at another alpha. */
    const MatrixXd dgam_dx_wake = WakeTracedAtOtherAlpha(state)
                                      ? DGammaDFoilXForWake(state)
                                      : ctx.dgam_dx;
    ctx.dwake_x_dx = DWakeXDFoilX(state, dfoil_s_dx, dgam_dx_wake);
    ctx.dwake_t_dx = DWakeTDFoilX(state, dfoil_s_dx, dgam_dx_wake);
    ctx.dwake_s_dx = DWakeSDFoilX(state, dfoil_s_dx);
    ctx.aic_t_lu.compute(MatrixXd(state.isol.AIC.transpose()));
    return ctx;
}

VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const MatrixXd& dfoil_s_dx) {
    return ComputeRueVJP(state, lambda_ue,
                         BuildRueVJPContext(state, dfoil_s_dx));
}

/**
 * @brief Build DRueDFoilX [Nsys x 2N] by applying ComputeRueVJP with each
 * standard basis vector as lambda. Slow (Nsys VJP calls) — for validation only.
 *
 * @param state        Fully assembled foil state.
 * @param dfoil_s_dx   Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 * @return dR_ue/d(foil.x) [Nsys x 2N].
 */
MatrixXd DRueDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    const RueVJPContext ctx = BuildRueVJPContext(state, dfoil_s_dx);
    MatrixXd dR_ue = MatrixXd::Zero(Nsys, 2 * N);
    VectorXd e = VectorXd::Zero(Nsys);
    for (int i = 0; i < Nsys; ++i) {
        e.setZero();
        e(i) = 1.0;
        dR_ue.row(i) = ComputeRueVJP(state, e, ctx).transpose();
    }
    return dR_ue;
}

} /* namespace scirocco */
