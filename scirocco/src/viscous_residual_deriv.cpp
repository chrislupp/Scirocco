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

    if (N < 2 || Nw < 1) {
        throw std::invalid_argument("ComputeRueVJP: state has no viscous wake");
    }
    if (lambda_ue.size() != Nsys) {
        throw std::runtime_error("ComputeRueVJP: lambda_ue size != Nsys");
    }
    if (ctx.dgam_dx.cols() != twoN || ctx.dgam_dx.rows() != N + 1 ||
        ctx.dwake_x_dx.rows() != 2 * Nw) {
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
    /* Contribution: -lambda^T * dueinv_foil = -lambda_g^T * dgam, with
     * dgam = d(gamma)/d(foil.x) from the context. (Equivalently
     * +AICResidualSensVJP(state, AIC^{-T} lambda_g, alpha).) */
    vjp -= (lambda_g_ext.transpose() * ctx.dgam_dx).transpose();
    const Eigen::PartialPivLU<MatrixXd>& lu = ctx.aic_t_lu;

    /* ---- C. dphi chain via dsigma_m ----------------------------------------
     * beta = ue_sigma^T * (sgue .* lambda_ue),  shape [Npanel = N+Nw-2] */
    VectorXd sgue_lambda(Nsys);
    for (int i = 0; i < N; ++i)
        sgue_lambda(i) = state.isol.sgnue(i) * lambda_ue(i);
    for (int i = N; i < Nsys; ++i) sgue_lambda(i) = lambda_ue(i);
    VectorXd beta =
        state.vsol.ue_sigma.transpose() * sgue_lambda; /* [Npanel] */

    /* phi = sigma_m * mass,  mass[k] = U[1, k] * U[3, k] */
    const VectorXd& phi = ctx.phi; /* [Npanel] */

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
     * wake.x[0].
     *
     * Row j of d(uewi)/d(foil.x) does not depend on lambda, so the context
     * holds it (together with B2a below) in ctx.wake_rows. */
    vjp += (lambda_ue.tail(Nw).transpose() * ctx.wake_rows).transpose();

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
             * Cgam[0, :] is replaced by Bp[N-1, :] in Dw[0], so only j > 0
             * rows enter. The context rebuilds Cgam, since CalcUeM doesn't
             * save it. */
            lambda_f_eff += (lambda_ue.tail(Nw - 1).transpose() *
                             ctx.cgam.bottomRows(Nw - 1))
                                .transpose();

            /* mu_f = A^{-T} * [lambda_f_eff; 0] */
            VectorXd lambda_f_ext = VectorXd::Zero(AIC_dim);
            lambda_f_ext.head(N) = lambda_f_eff;
            VectorXd mu_f = lu.solve(lambda_f_ext);

            /* B1a + B1b: mu_f^T * d(A*v + B*phi)/d(foil.x), v = [Bp*phi; 0]. */
            vjp += (mu_f.head(N).transpose() * ctx.dsrc_dx).transpose();

            /* B1b: mu_f^T * dB * phi, foil and wake source panels. Its rows
             * do not depend on lambda; the context assembles them together
             * with B1a into dsrc_dx = d(A*v + B*phi)/d(foil.x). */
        }
    }
    /* ---- B2a: dCgam direct chain for j > 0 wake rows -----------------------
     *
     * For each wake station j > 0:
     *   contribution = -lambda_w[j] * d(V_for_v . wake.t[j])/dfoil_x
     * where V_for_v is the gamma-induced velocity at wake.x[j] with vortex
     * strengths v = Bp * phi held fixed. Row j does not depend on lambda, so
     * it is folded into ctx.wake_rows (applied with A2 above).
     *
     * For j = 0, Cgam[0, :] is replaced by Bp[N-1, :] in Dw[0]; handled in B1.
     */
    /* ---- B2b, B2c: dCsig chains for j > 0 wake rows ----------------------
     * Linear in lambda_w with lambda-independent rows; folded into
     * ctx.wake_rows by BuildRueVJPContext. */
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

    const int N = state.foil.N;
    const int Nw = state.wake.N;
    const int Nsys = N + Nw;
    const int twoN = 2 * N;
    const int AIC_dim = static_cast<int>(state.isol.AIC.rows());
    const double alpha = state.oper.alpha;

    /* phi = sigma_m * mass,  mass[k] = U[1, k] * U[3, k] */
    VectorXd mass(Nsys);
    for (int i = 0; i < Nsys; ++i)
        mass(i) = state.glob.U(1, i) * state.glob.U(3, i);
    ctx.phi = state.vsol.sigma_m * mass;

    /* wake_rows[j] = d(lambda_w[j] terms)/d(foil.x) for j > 0, applied as
     * vjp += wake_rows^T lambda_w. A2 and B2a enter with a minus sign.
     *
     * A2: row j of d(uewi[j])/d(foil.x) = t^T dV + V^T dwake.t, j > 0. */
    ctx.wake_rows = MatrixXd::Zero(Nw, twoN);
    for (int j = 1; j < Nw; ++j) {
        VelDeriv vd = InviscidVelocityWithDeriv(state, state.wake.x.col(j),
                                                ctx.dgam_dx, alpha);
        MatrixXd dwxj(2, twoN);
        dwxj.row(0) = ctx.dwake_x_dx.row(0 * Nw + j);
        dwxj.row(1) = ctx.dwake_x_dx.row(1 * Nw + j);
        MatrixXd dwt(2, twoN);
        dwt.row(0) = ctx.dwake_t_dx.row(0 * Nw + j);
        dwt.row(1) = ctx.dwake_t_dx.row(1 * Nw + j);
        MatrixXd dV = vd.V_xfoil + vd.V_xi * dwxj; /* [2 x 2N] */
        Vector2d t_wake = state.wake.t.col(j);
        ctx.wake_rows.row(j) -= t_wake.transpose() * dV;
        /* V^T * dwake.t (zero for j < Nw-1 since V || t). */
        ctx.wake_rows.row(j) -= vd.V.transpose() * dwt;
    }

    /* The per-output loops (B1b, B2b, B2c) are linear in lambda with
     * lambda-independent coefficients, so their rows are assembled once here.
     * The loop bodies are those of the lambda form with the lambda weight set
     * to 1: scatters into foil.x go to R(row, :), and chains through wake.x
     * and wake.t collect 2-vector weights in GX and GT, which multiply the
     * wake Jacobians at the end. Column d*Nw + j matches their row layout. */
    const VectorXd& phi = ctx.phi;
    int row = 0;
    MatrixXd R;
    MatrixXd GX;
    MatrixXd GT;
    auto add_x = [&](int j, const Vector2d& g) {
        GX(row, j) += g(0);
        GX(row, Nw + j) += g(1);
    };
    auto add_t = [&](int j, const Vector2d& g) {
        GT(row, j) += g(0);
        GT(row, Nw + j) += g(1);
    };
    /* Weight on the midpoint 0.5*(wake.x[a] + wake.x[b]). */
    auto add_mid = [&](int a, int b, const Vector2d& g) {
        add_x(a, 0.5 * g);
        add_x(b, 0.5 * g);
    };
    /* Weight on the ghost point 2*wake.x[b] - 0.5*(wake.x[a] + wake.x[b]). */
    auto add_ghost = [&](int a, int b, const Vector2d& g) {
        add_x(b, 1.5 * g);
        add_x(a, -0.5 * g);
    };

    const MatrixXd& ue_sig = state.vsol.ue_sigma;
    const int Npanel = N + Nw - 2;
    if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
        /* v = Bp * phi; ue_sigma's top N rows are Bp. */
        VectorXd v_ext = VectorXd::Zero(AIC_dim);
        v_ext.head(N) = ue_sig.topRows(N) * ctx.phi;
        const MatrixXd daic_v_dx = AICMatVecJacobianX(state, v_ext);

        /* B1a + B1b rows: d(A*v + B*phi)/d(foil.x), foil rows only. */
        R = MatrixXd::Zero(N, twoN);
        GX = MatrixXd::Zero(N, 2 * Nw);
        /* B1b: mu_f^T * dB * phi
         *   B columns 0..N-2: foil constant-source panels (foil.x[k:k+2])
         *   B columns N-1..Npanel-1: wake linear-source panels (wake.x)
         *   For each eval i and panel k, accumulate:
         *     sens += mu_f[i] * phi[k] * dB[i, k]/dfoil_x
         *
         * Foil source panels first (no wake.x dependence — pure foil
         * chain): */
        for (int i = 0; i < N; ++i) {
            const double mu_i = 1.0; /* row i of the Jacobian */
            row = i;
            Vector2d xi = state.foil.x.col(i);
            for (int k = 0; k < N - 1; ++k) {
                double phi_k = phi(k);
                if (std::abs(phi_k) < 1e-20) continue;
                Matrix2d Xj;
                Xj.col(0) = state.foil.x.col(k);
                Xj.col(1) = state.foil.x.col(k + 1);
                ConstSourceStreamDeriv d = PanelConstSourceStreamDeriv(Xj, xi);
                double w = mu_i * phi_k;
                /* xi derivatives (foil.x[i]) */
                R(row, 2 * i) += w * d.val_xi(0);
                R(row, 2 * i + 1) += w * d.val_xi(1);
                /* Xj derivatives (foil.x[k], foil.x[k+1]) */
                R(row, 2 * k) += w * d.val_Xj(0);
                R(row, 2 * k + 1) += w * d.val_Xj(1);
                R(row, 2 * (k + 1)) += w * d.val_Xj(2);
                R(row, 2 * (k + 1) + 1) += w * d.val_Xj(3);
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
         * right). Chains through wake.x go to the accumulators. */
        for (int i = 0; i < N; ++i) {
            const double mu_i = 1.0; /* row i of the Jacobian */
            row = i;
            Vector2d xi = state.foil.x.col(i);
            for (int j = 0; j < Nw - 1; ++j) {
                Vector2d Xm =
                    0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));

                Vector2d X_right;
                const bool ghost_right = (j == Nw - 2);
                if (ghost_right) {
                    X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
                } else {
                    X_right = state.wake.x.col(j + 1);
                }

                /* Left half panel */
                Matrix2d Xj_L;
                Xj_L.col(0) = state.wake.x.col(j);
                Xj_L.col(1) = Xm;
                LinSourceStreamDeriv dL = PanelLinSourceStreamDeriv(Xj_L, xi);
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
                R(row, 2 * i) += w_a * dL.a_xi(0) + w_b * dL.b_xi(0);
                R(row, 2 * i + 1) += w_a * dL.a_xi(1) + w_b * dL.b_xi(1);
                /* Xj chain — Xj cols 0,1 = first endpoint (wake.x[j]),
                 * cols 2,3 = second endpoint (Xm). */
                Vector4d dab = w_a * dL.a_Xj + w_b * dL.b_Xj;
                /* First endpoint = wake.x[j], second = Xm. */
                add_x(j, dab.head<2>());
                add_mid(j, j + 1, dab.tail<2>());

                /* Right half panel */
                Matrix2d Xj_R;
                Xj_R.col(0) = Xm;
                Xj_R.col(1) = X_right;
                LinSourceStreamDeriv dR = PanelLinSourceStreamDeriv(Xj_R, xi);
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
                R(row, 2 * i) += w_a * dR.a_xi(0) + w_b * dR.b_xi(0);
                R(row, 2 * i + 1) += w_a * dR.a_xi(1) + w_b * dR.b_xi(1);
                dab = w_a * dR.a_Xj + w_b * dR.b_Xj;
                /* First endpoint = Xm, second = X_right. */
                add_mid(j, j + 1, dab.head<2>());
                if (ghost_right) {
                    add_ghost(j, j + 1, dab.tail<2>());
                } else {
                    add_x(j + 1, dab.tail<2>());
                }
            }
        }
        ctx.dsrc_dx = daic_v_dx.topRows(N) + R + GX * ctx.dwake_x_dx;

        /* Cgam[j, k]: gamma_k's induced velocity at wake.x[j] along t[j]. */
        ctx.cgam.resize(Nw, N);
        for (int j = 0; j < Nw; ++j) {
            Vector2d xi_w = state.wake.x.col(j);
            Vector2d ti_w = state.wake.t.col(j);
            MatrixXd V_G;
            InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi_w, true,
                             &V_G);
            for (int k = 0; k < N; ++k) {
                ctx.cgam(j, k) = V_G(0, k) * ti_w(0) + V_G(1, k) * ti_w(1);
            }
        }

        /* B2a: row j of d(V_for_v . wake.t[j])/d(foil.x), j > 0, where
         * V_for_v is the velocity induced by vortex strengths v (held fixed).
         * InviscidVelocityWithDerivExt includes the freestream in .V, so
         * subtract it: Cgam[j, :] * v is gamma-induced only. */
        const MatrixXd dgam_zero = MatrixXd::Zero(AIC_dim, twoN);
        const Vector2d Vinf_vec(state.oper.Vinf * Cosd(alpha),
                                state.oper.Vinf * Sind(alpha));
        for (int j = 1; j < Nw; ++j) {
            Vector2d xi_w = state.wake.x.col(j);
            Vector2d t_wake = state.wake.t.col(j);
            VelDeriv vd_v = InviscidVelocityWithDerivExt(state, xi_w, v_ext,
                                                         dgam_zero, alpha);
            Vector2d V_G_v = vd_v.V - Vinf_vec;
            MatrixXd dwxj(2, twoN);
            dwxj.row(0) = ctx.dwake_x_dx.row(0 * Nw + j);
            dwxj.row(1) = ctx.dwake_x_dx.row(1 * Nw + j);
            MatrixXd dwt(2, twoN);
            dwt.row(0) = ctx.dwake_t_dx.row(0 * Nw + j);
            dwt.row(1) = ctx.dwake_t_dx.row(1 * Nw + j);
            MatrixXd dV_for_v = vd_v.V_xfoil + vd_v.V_xi * dwxj;
            ctx.wake_rows.row(j) -= t_wake.transpose() * dV_for_v;
            ctx.wake_rows.row(j) -= V_G_v.transpose() * dwt;
        }
    }

    /* B2b and B2c rows. */
    R = MatrixXd::Zero(Nw, twoN);
    GX = MatrixXd::Zero(Nw, 2 * Nw);
    GT = MatrixXd::Zero(Nw, 2 * Nw);
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
                const double lw = 1.0; /* row j of the wake block */
                row = j;
                Vector2d xi_w = state.wake.x.col(j);
                Vector2d t_w = state.wake.t.col(j);

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
                        R(row, 2 * k + c) += w * d_proj;
                    }
                    for (int c = 2; c < 4; ++c) {
                        double d_proj = csv.val_Xj(0, c) * t_w(0) +
                                        csv.val_Xj(1, c) * t_w(1);
                        R(row, 2 * (k + 1) + (c - 2)) += w * d_proj;
                    }

                    /* Chain 2: xi derivatives chained through dwake.x[j]
                     * d(V.t)/dxi [1x2] = t^T * val_xi */
                    Vector2d dprojdxi(
                        csv.val_xi(0, 0) * t_w(0) + csv.val_xi(1, 0) * t_w(1),
                        csv.val_xi(0, 1) * t_w(0) + csv.val_xi(1, 1) * t_w(1));
                    add_x(j, w * dprojdxi);

                    /* Chain 3: V_panel . dwake.t[j] */
                    add_t(j, w * V_panel);
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
     * Chain through wake.x (panel coords + midpoint coords) and wake.t.
     */
    {
        const MatrixXd& ue_sig = state.vsol.ue_sigma;
        int Npanel = N + Nw - 2;
        if (ue_sig.cols() == Npanel && ue_sig.rows() == Nsys) {
            for (int j_eval = 1; j_eval < Nw; ++j_eval) {
                const double lw = 1.0; /* row j_eval of the wake block */
                row = j_eval;
                Vector2d xi_w = state.wake.x.col(j_eval);
                Vector2d t_w = state.wake.t.col(j_eval);

                /* Self-influence: interior wp == j_eval case (skipped for
                 * endpoints). */
                if (j_eval >= 1 && j_eval <= Nw - 2) {
                    int wp = j_eval;
                    Vector2d Xl =
                        0.5 * (state.wake.x.col(wp - 1) + state.wake.x.col(wp));
                    Vector2d Xm = state.wake.x.col(wp);
                    Vector2d Xr =
                        0.5 * (state.wake.x.col(wp) + state.wake.x.col(wp + 1));
                    Vector2d e1 = Xm - Xl;
                    Vector2d e2 = Xr - Xm;
                    double d1 = e1.norm();
                    double d2 = e2.norm();
                    if (d1 > 1e-14 && d2 > 1e-14) {
                        Vector2d e1u = e1 / d1;
                        Vector2d e2u = e2 / d2;
                        /* aa = (0.25/pi) * (log-term), with
                         *   d(aa) = (0.25/pi) * (e1u.de1/d1 - e2u.de2/d2),
                         *   de1 = dXm - dXl, de2 = dXr - dXm. */
                        int col_prev = N - 1 + wp - 1;
                        int col_curr = N - 1 + wp;
                        double w = -lw * (phi(col_prev) + phi(col_curr));
                        const Vector2d g1 = (w * 0.25 / kPi / d1) * e1u;
                        const Vector2d g2 = (-w * 0.25 / kPi / d2) * e2u;
                        add_x(wp, g1 - g2);       /* Xm */
                        add_mid(wp - 1, wp, -g1); /* Xl */
                        add_mid(wp, wp + 1, g2);  /* Xr */
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
                        Matrix2d Xj;
                        Xj.col(0) = state.wake.x.col(0);
                        Xj.col(1) = Xm0;
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
                                const LinSourceVelocityDeriv& d) {
                                /* Chain 1: Xj derivatives. First endpoint is
                                 * wake.x[0], second is Xm0. */
                                Vector2d g0, g1;
                                for (int c = 0; c < 2; ++c) {
                                    double pa = wa * (d.a_Xj(0, c) * t_w(0) +
                                                      d.a_Xj(1, c) * t_w(1));
                                    double pb = wb * (d.b_Xj(0, c) * t_w(0) +
                                                      d.b_Xj(1, c) * t_w(1));
                                    g0(c) = pa + pb;
                                }
                                for (int c = 2; c < 4; ++c) {
                                    double pa = wa * (d.a_Xj(0, c) * t_w(0) +
                                                      d.a_Xj(1, c) * t_w(1));
                                    double pb = wb * (d.b_Xj(0, c) * t_w(0) +
                                                      d.b_Xj(1, c) * t_w(1));
                                    g1(c - 2) = pa + pb;
                                }
                                add_x(0, g0);
                                add_mid(0, 1, g1);
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
                                add_x(j_eval, dpa_dxi + dpb_dxi);
                            };
                        add_velocity_chain(w_a, w_b, lvd);
                        /* V_panel . dwake.t[j_eval] chain.
                         * V contributions to Csig cols:
                         *   col(N-1): b   col(0): a   col(N-2): a
                         * V_total = phi[N-1]*b + (phi[0]+phi[N-2])*a */
                        {
                            Vector2d V_total = phi(N - 1) * lvd.b +
                                               (phi(0) + phi(N - 2)) * lvd.a;
                            add_t(j_eval, -lw * V_total);
                        }
                    } else if (wp == Nw - 1) {
                        /* Last wake panel: full constant-source from [Xl, Xr]
                         * Csig col (N + Nw - 3 = Npanel - 1) gets v.dot(ti) */
                        Vector2d Xl = 0.5 * (state.wake.x.col(wp - 1) +
                                             state.wake.x.col(wp));
                        Vector2d Xm = state.wake.x.col(wp);
                        Vector2d Xr = 2.0 * Xm - Xl;
                        Matrix2d Xj;
                        Xj.col(0) = Xl;
                        Xj.col(1) = Xr;
                        ConstSourceVelocityDeriv csv =
                            PanelConstSourceVelocityDeriv(Xj, xi_w);
                        int col_last = Npanel - 1;
                        double w = -lw * phi(col_last);
                        /* Chain 1: Xj derivatives. Xl is a midpoint, Xr
                         * the ghost point 2*Xm - Xl. */
                        Vector2d gl, gr;
                        for (int c = 0; c < 2; ++c) {
                            gl(c) = w * (csv.val_Xj(0, c) * t_w(0) +
                                         csv.val_Xj(1, c) * t_w(1));
                        }
                        for (int c = 2; c < 4; ++c) {
                            gr(c - 2) = w * (csv.val_Xj(0, c) * t_w(0) +
                                             csv.val_Xj(1, c) * t_w(1));
                        }
                        add_mid(wp - 1, wp, gl);
                        add_ghost(wp - 1, wp, gr);
                        /* Chain 2: xi via wake.x[j_eval] */
                        Vector2d d_proj_xi(csv.val_xi(0, 0) * t_w(0) +
                                               csv.val_xi(1, 0) * t_w(1),
                                           csv.val_xi(0, 1) * t_w(0) +
                                               csv.val_xi(1, 1) * t_w(1));
                        add_x(j_eval, w * d_proj_xi);
                        /* Chain 3: V_panel . dwake.t[j_eval] */
                        add_t(j_eval, w * csv.val);
                    } else {
                        /* Interior wake panel wp: left half [Xl, Xm] + right
                         * half [Xm, Xr] */
                        Vector2d Xl = 0.5 * (state.wake.x.col(wp - 1) +
                                             state.wake.x.col(wp));
                        Vector2d Xm = state.wake.x.col(wp);
                        Vector2d Xr = 0.5 * (state.wake.x.col(wp) +
                                             state.wake.x.col(wp + 1));

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
                        Vector2d gl, gm;
                        for (int c = 0; c < 2; ++c) {
                            double pa = w_a1 * (dL.a_Xj(0, c) * t_w(0) +
                                                dL.a_Xj(1, c) * t_w(1));
                            double pb = w_b1 * (dL.b_Xj(0, c) * t_w(0) +
                                                dL.b_Xj(1, c) * t_w(1));
                            gl(c) = pa + pb;
                        }
                        for (int c = 2; c < 4; ++c) {
                            double pa = w_a1 * (dL.a_Xj(0, c) * t_w(0) +
                                                dL.a_Xj(1, c) * t_w(1));
                            double pb = w_b1 * (dL.b_Xj(0, c) * t_w(0) +
                                                dL.b_Xj(1, c) * t_w(1));
                            gm(c - 2) = pa + pb;
                        }
                        add_mid(wp - 1, wp, gl);
                        add_x(wp, gm);
                        Vector2d dxi_L(w_a1 * (dL.a_xi(0, 0) * t_w(0) +
                                               dL.a_xi(1, 0) * t_w(1)) +
                                           w_b1 * (dL.b_xi(0, 0) * t_w(0) +
                                                   dL.b_xi(1, 0) * t_w(1)),
                                       w_a1 * (dL.a_xi(0, 1) * t_w(0) +
                                               dL.a_xi(1, 1) * t_w(1)) +
                                           w_b1 * (dL.b_xi(0, 1) * t_w(0) +
                                                   dL.b_xi(1, 1) * t_w(1)));
                        add_x(j_eval, dxi_L);

                        /* Right half [Xm, Xr] */
                        Matrix2d Xj_R;
                        Xj_R.col(0) = Xm;
                        Xj_R.col(1) = Xr;
                        LinSourceVelocityDeriv dR =
                            PanelLinSourceVelocityDeriv(Xj_R, xi_w);
                        Vector2d gm2, gr;
                        for (int c = 0; c < 2; ++c) {
                            double pa = w_a2 * (dR.a_Xj(0, c) * t_w(0) +
                                                dR.a_Xj(1, c) * t_w(1));
                            double pb = w_b2 * (dR.b_Xj(0, c) * t_w(0) +
                                                dR.b_Xj(1, c) * t_w(1));
                            gm2(c) = pa + pb;
                        }
                        for (int c = 2; c < 4; ++c) {
                            double pa = w_a2 * (dR.a_Xj(0, c) * t_w(0) +
                                                dR.a_Xj(1, c) * t_w(1));
                            double pb = w_b2 * (dR.b_Xj(0, c) * t_w(0) +
                                                dR.b_Xj(1, c) * t_w(1));
                            gr(c - 2) = pa + pb;
                        }
                        add_x(wp, gm2);
                        add_mid(wp, wp + 1, gr);
                        Vector2d dxi_R(w_a2 * (dR.a_xi(0, 0) * t_w(0) +
                                               dR.a_xi(1, 0) * t_w(1)) +
                                           w_b2 * (dR.b_xi(0, 0) * t_w(0) +
                                                   dR.b_xi(1, 0) * t_w(1)),
                                       w_a2 * (dR.a_xi(0, 1) * t_w(0) +
                                               dR.a_xi(1, 1) * t_w(1)) +
                                           w_b2 * (dR.b_xi(0, 1) * t_w(0) +
                                                   dR.b_xi(1, 1) * t_w(1)));
                        add_x(j_eval, dxi_R);
                        /* V_panel . dwake.t[j_eval] chain (left + right half).
                         */
                        {
                            Vector2d V_left =
                                phi_prev * dL.a + 0.5 * phi_prev * dL.b;
                            Vector2d V_right =
                                0.5 * phi_curr * dR.a + phi_curr * dR.b;
                            Vector2d V_total = V_left + V_right;
                            add_t(j_eval, -lw * V_total);
                        }
                    }
                }
            }
        }
    }

    ctx.wake_rows += R + GX * ctx.dwake_x_dx + GT * ctx.dwake_t_dx;
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
