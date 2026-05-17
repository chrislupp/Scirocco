/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "cppfoil/core/inviscid_deriv.h"

#include <cmath>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/panel_influence_deriv.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

VectorXd AICResidualSensVJP(const FoilState& state, const VectorXd& lambda,
                             double alpha) {
    // Computes lambda^T * d(A*gamma - rhs)/dx_foil  [2*N]
    //
    // The inviscid residual is: R = A*gamma - rhs = 0
    // where A is the AIC matrix, gamma are the vortex strengths, and
    // rhs depends on the evaluation point coordinates.
    //
    // Perturbing foil node k affects:
    //   (a) Row k: xi changes for all panels (evaluation point)
    //   (b) Panels k-1->k and k->k+1: Xj endpoint changes
    //   (c) TE source/vortex panels (nodes 0 and N-1)
    //   (d) rhs: rhs(i,0) = -xi(1), rhs(i,1) = xi(0)
    //   (e) Kutta condition row N: no coordinate dependence
    //   (f) No-gap row N-1: depends on gamma extrapolation (no coord dependence)

    int N = state.foil.N;
    VectorXd sens = VectorXd::Zero(2 * N);

    const VectorXd& gam = state.isol.gam;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    TEInfo te = GetTEInfo(state.foil.x);
    bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);

    // For each node k, accumulate:
    // sens(2*k + 0) = d/d(x_k_x), sens(2*k + 1) = d/d(x_k_z)
    //
    // The VJP is: sum_i lambda(i) * dR(i)/dx_k for each k

    // --- (1) Linear vortex panel contributions ---
    // Row i, panel j (endpoints j, j+1):
    //   A(i,j) += a, A(i,j+1) += b from PanelLinVortexStream(Xj, xi)
    //   R(i) += a*gam(j) + b*gam(j+1)
    //
    // When we perturb node k:
    //   - As xi (row i=k): dR(k)/dx_k via dxi for ALL panels j
    //   - As Xj endpoint: for panel j=k-1 (endpoint x2=k) and j=k (endpoint x1=k)

    for (int i = 0; i < N; ++i) {
        if (nogap && i == N - 1) continue;  // Skip no-gap row
        double lam_i = lambda(i);
        if (std::abs(lam_i) < 1e-20) continue;

        Vector2d xi = state.foil.x.col(i);

        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);

            auto deriv = PanelLinVortexStreamDeriv(Xj, xi);

            // Contribution to gamma: dR(i)/d_something += da*gam(j) + db*gam(j+1)
            double gj = gam(j);
            double gjp = gam(j + 1);

            // xi derivatives (affect node i's coordinates)
            // sens(2*i+0) += lam_i * (da/dxi_x * gj + db/dxi_x * gjp)
            // sens(2*i+1) += lam_i * (da/dxi_z * gj + db/dxi_z * gjp)
            sens(2 * i) += lam_i * (deriv.a_xi(0) * gj + deriv.b_xi(0) * gjp);
            sens(2 * i + 1) += lam_i * (deriv.a_xi(1) * gj + deriv.b_xi(1) * gjp);

            // Xj derivatives (affect panel endpoint coordinates)
            // Xj = [x1_x, x1_z, x2_x, x2_z] = node j and j+1
            Vector4d dR_dXj = deriv.a_Xj * gj + deriv.b_Xj * gjp;
            sens(2 * j) += lam_i * dR_dXj(0);       // x1_x = node j x
            sens(2 * j + 1) += lam_i * dR_dXj(1);   // x1_z = node j z
            sens(2 * (j + 1)) += lam_i * dR_dXj(2);     // x2_x = node j+1 x
            sens(2 * (j + 1) + 1) += lam_i * dR_dXj(3); // x2_z = node j+1 z
        }

        // --- (2) TE source panel contribution ---
        // TE panel: Xj = [node N-1, node 0]
        // A(i,0) += -a_src * (0.5*tcp)
        // A(i,N-1) += a_src * (0.5*tcp)
        {
            Matrix2d XjTE;
            XjTE.col(0) = state.foil.x.col(N - 1);
            XjTE.col(1) = state.foil.x.col(0);

            auto deriv = PanelConstSourceStreamDeriv(XjTE, xi);

            // R contribution: a_src * 0.5*tcp * (gam(N-1) - gam(0))
            double sigma_factor = 0.5 * te.tcp * (gam(N - 1) - gam(0));

            // xi derivatives
            sens(2 * i) += lam_i * deriv.val_xi(0) * sigma_factor;
            sens(2 * i + 1) += lam_i * deriv.val_xi(1) * sigma_factor;

            // Xj derivatives for TE panel: Xj = [N-1, 0]
            // (dtcp/dx and dtdp/dx are handled separately in section (5) below)
            sens(2 * (N - 1)) += lam_i * deriv.val_Xj(0) * sigma_factor;
            sens(2 * (N - 1) + 1) += lam_i * deriv.val_Xj(1) * sigma_factor;
            sens(0) += lam_i * deriv.val_Xj(2) * sigma_factor;
            sens(1) += lam_i * deriv.val_Xj(3) * sigma_factor;
        }

        // --- (3) TE vortex panel contribution ---
        {
            Matrix2d XjTE;
            XjTE.col(0) = state.foil.x.col(N - 1);
            XjTE.col(1) = state.foil.x.col(0);

            auto deriv = PanelLinVortexStreamDeriv(XjTE, xi);

            // The TE vortex panel uses (a+b) as a uniform vortex with
            // strength factor: 0.5*tdp for node 0 and -0.5*tdp for node N-1
            // R contribution from TE vortex:
            // A(i,0) += -(a+b)*(-0.5*tdp) = (a+b)*0.5*tdp
            // A(i,N-1) += (a+b)*(-0.5*tdp)
            double gam_factor = 0.5 * te.tdp * (gam(0) - gam(N - 1));
            // Total R contribution = (a+b) * gam_factor

            // xi derivatives (a_xi and b_xi are Vector2d: [d/dxi_x, d/dxi_z])
            double dab_dxi_x = deriv.a_xi(0) + deriv.b_xi(0);
            double dab_dxi_z = deriv.a_xi(1) + deriv.b_xi(1);
            sens(2 * i) += lam_i * dab_dxi_x * gam_factor;
            sens(2 * i + 1) += lam_i * dab_dxi_z * gam_factor;

            // Xj derivatives: Xj = [N-1, 0]
            Vector4d dab_dXj = deriv.a_Xj + deriv.b_Xj;
            sens(2 * (N - 1)) += lam_i * dab_dXj(0) * gam_factor;
            sens(2 * (N - 1) + 1) += lam_i * dab_dXj(1) * gam_factor;
            sens(0) += lam_i * dab_dXj(2) * gam_factor;
            sens(1) += lam_i * dab_dXj(3) * gam_factor;
        }

        // --- (4) RHS derivatives ---
        // rhs(i,0) = -xi(1), rhs(i,1) = xi(0)
        // rhs(i) = rhs(i,0)*cos(alpha) + rhs(i,1)*sin(alpha)
        //        = -xi(1)*cos(alpha) + xi(0)*sin(alpha)
        // R = A*gam - rhs, so dR/dx_k_x = -drhs/dx_k_x at row i=k
        // drhs(k)/dx_k_x = sin(alpha), drhs(k)/dx_k_z = -cos(alpha)
        // => -drhs/dx_k_x = -sin(alpha), -drhs/dx_k_z = cos(alpha)
        sens(2 * i) += lam_i * (-sina);
        sens(2 * i + 1) += lam_i * cosa;
    }

    // --- (5) TE factor derivatives: d(tcp)/dx and d(tdp)/dx ---
    // The TE source/vortex panel strengths are proportional to tcp/tdp,
    // which depend on nodes 0, 1, N-2, N-1.
    {
        Vector2d s = state.foil.x.col(N - 1) - state.foil.x.col(0);
        double sn = s.norm();
        if (sn > 1e-10) {
            Matrix2d XjTE;
            XjTE.col(0) = state.foil.x.col(N - 1);
            XjTE.col(1) = state.foil.x.col(0);

            // Precompute TE influence at each node, weighted by lambda
            // TE source: sum_i lambda(i) * a_src(i) * 0.5 * (gam(N-1)-gam(0))
            // TE vortex: sum_i lambda(i) * (a_vor+b_vor)(i) * 0.5 * (gam(0)-gam(N-1))
            double tcp_sum = 0.0, tdp_sum = 0.0;
            for (int i = 0; i < N; ++i) {
                if (std::abs(lambda(i)) < 1e-20) continue;
                Vector2d xi_i = state.foil.x.col(i);
                double a_src = PanelConstSourceStream(XjTE, xi_i);
                tcp_sum += lambda(i) * a_src;
                double a_v, b_v;
                PanelLinVortexStream(XjTE, xi_i, a_v, b_v);
                tdp_sum += lambda(i) * (a_v + b_v);
            }
            tcp_sum *= 0.5 * (gam(N - 1) - gam(0));
            tdp_sum *= 0.5 * (gam(0) - gam(N - 1));

            // For each TE-affecting node, add dtcp*tcp_sum + dtdp*tdp_sum
            int te_nodes[] = {0, 1, N - 2, N - 1};
            for (int kk : te_nodes) {
                if (kk < 0 || kk >= N) continue;
                for (int dd = 0; dd < 2; ++dd) {
                    // Compute d(tcp)/d(x[kk]_dd) and d(tdp)/d(x[kk]_dd)
                    // Inline the derivative computation
                    Vector2d v1 = state.foil.x.col(0) - state.foil.x.col(1);
                    Vector2d v2 = state.foil.x.col(N-1) - state.foil.x.col(N-2);
                    double v1n = v1.norm(), v2n = v2.norm();
                    if (v1n < 1e-10 || v2n < 1e-10) continue;
                    Vector2d t1 = v1/v1n, t2 = v2/v2n;
                    Vector2d traw = 0.5*(t1+t2);
                    double tnorm = traw.norm();
                    if (tnorm < 1e-10) continue;
                    Vector2d tvec = traw/tnorm;
                    Vector2d pvec = s/sn;

                    // dt1/d(x[kk]_dd)
                    Matrix2d dt1_dv1 = (Matrix2d::Identity() - t1*t1.transpose())/v1n;
                    Vector2d dv1 = Vector2d::Zero();
                    if (kk == 0) dv1(dd) = 1.0;
                    else if (kk == 1) dv1(dd) = -1.0;
                    Vector2d dt1_d = dt1_dv1 * dv1;

                    Matrix2d dt2_dv2 = (Matrix2d::Identity() - t2*t2.transpose())/v2n;
                    Vector2d dv2 = Vector2d::Zero();
                    if (kk == N-1) dv2(dd) = 1.0;
                    else if (kk == N-2) dv2(dd) = -1.0;
                    Vector2d dt2_d = dt2_dv2 * dv2;

                    Vector2d dtraw_d = 0.5*(dt1_d + dt2_d);
                    Matrix2d dt_dtraw = (Matrix2d::Identity() - tvec*tvec.transpose())/tnorm;
                    Vector2d dt_d = dt_dtraw * dtraw_d;

                    Vector2d ds_d = Vector2d::Zero();
                    if (kk == N-1) ds_d(dd) = 1.0;
                    else if (kk == 0) ds_d(dd) = -1.0;
                    Matrix2d dp_ds = (Matrix2d::Identity() - pvec*pvec.transpose())/sn;
                    Vector2d dp_d = dp_ds * ds_d;

                    double cross_val = tvec(0)*pvec(1) - tvec(1)*pvec(0);
                    double sign_c = (cross_val >= 0) ? 1.0 : -1.0;
                    double dcross = dt_d(0)*pvec(1) + tvec(0)*dp_d(1)
                                  - dt_d(1)*pvec(0) - tvec(1)*dp_d(0);
                    double dtcp_d = sign_c * dcross;
                    double dtdp_d = dt_d.dot(pvec) + tvec.dot(dp_d);

                    sens(2*kk + dd) += dtcp_d * tcp_sum + dtdp_d * tdp_sum;
                }
            }
        }
    }

    return sens;
}

// ============================================================================
// AICMatVecVJP — lambda^T * d(A*v)/dx_foil
//
// Same panel-chain logic as AICResidualSensVJP but:
//   - gam → v (caller supplies any vector of length N+1)
//   - no rhs term (we differentiate only A*v, not A*v - rhs)
//
// Used by Phase 2b to chain through Bp = -A^{-1} * B:
//   dBp = -A^{-1} (dA * Bp + dB)
//   lambda_f^T * dBp * phi reduces to
//     mu^T (dA * v + dB * phi)  where v = Bp * phi, mu = A^{-T} lambda_f.
// This function computes the mu^T dA v half.
// ============================================================================
VectorXd AICMatVecVJP(const FoilState& state, const VectorXd& lambda,
                      const VectorXd& v) {
    int N = state.foil.N;
    VectorXd sens = VectorXd::Zero(2 * N);

    TEInfo te = GetTEInfo(state.foil.x);
    bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);

    for (int i = 0; i < N; ++i) {
        if (nogap && i == N - 1) continue;
        double lam_i = lambda(i);
        if (std::abs(lam_i) < 1e-20) continue;
        Vector2d xi = state.foil.x.col(i);

        // (1) Linear vortex panels
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);
            auto deriv = PanelLinVortexStreamDeriv(Xj, xi);
            double vj = v(j), vjp = v(j + 1);
            sens(2 * i)     += lam_i * (deriv.a_xi(0) * vj + deriv.b_xi(0) * vjp);
            sens(2 * i + 1) += lam_i * (deriv.a_xi(1) * vj + deriv.b_xi(1) * vjp);
            Vector4d dR_dXj = deriv.a_Xj * vj + deriv.b_Xj * vjp;
            sens(2 * j)         += lam_i * dR_dXj(0);
            sens(2 * j + 1)     += lam_i * dR_dXj(1);
            sens(2 * (j + 1))   += lam_i * dR_dXj(2);
            sens(2 * (j + 1) + 1) += lam_i * dR_dXj(3);
        }

        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);

        // (2) TE source panel
        {
            auto deriv = PanelConstSourceStreamDeriv(XjTE, xi);
            double sigma_factor = 0.5 * te.tcp * (v(N - 1) - v(0));
            sens(2 * i)         += lam_i * deriv.val_xi(0) * sigma_factor;
            sens(2 * i + 1)     += lam_i * deriv.val_xi(1) * sigma_factor;
            sens(2 * (N - 1))   += lam_i * deriv.val_Xj(0) * sigma_factor;
            sens(2 * (N - 1) + 1) += lam_i * deriv.val_Xj(1) * sigma_factor;
            sens(0)             += lam_i * deriv.val_Xj(2) * sigma_factor;
            sens(1)             += lam_i * deriv.val_Xj(3) * sigma_factor;
        }

        // (3) TE vortex panel
        {
            auto deriv = PanelLinVortexStreamDeriv(XjTE, xi);
            double v_factor = 0.5 * te.tdp * (v(0) - v(N - 1));
            double dab_dxi_x = deriv.a_xi(0) + deriv.b_xi(0);
            double dab_dxi_z = deriv.a_xi(1) + deriv.b_xi(1);
            sens(2 * i)         += lam_i * dab_dxi_x * v_factor;
            sens(2 * i + 1)     += lam_i * dab_dxi_z * v_factor;
            Vector4d dab_dXj = deriv.a_Xj + deriv.b_Xj;
            sens(2 * (N - 1))   += lam_i * dab_dXj(0) * v_factor;
            sens(2 * (N - 1) + 1) += lam_i * dab_dXj(1) * v_factor;
            sens(0)             += lam_i * dab_dXj(2) * v_factor;
            sens(1)             += lam_i * dab_dXj(3) * v_factor;
        }

        // (4) NO rhs term (intentionally skipped — we want d(A*v) only).
    }

    // (5) TE factor derivatives (tcp, tdp) — same as AICResidualSensVJP
    // but with v in place of gam.
    {
        Vector2d s = state.foil.x.col(N - 1) - state.foil.x.col(0);
        double sn = s.norm();
        if (sn > 1e-10) {
            Matrix2d XjTE;
            XjTE.col(0) = state.foil.x.col(N - 1);
            XjTE.col(1) = state.foil.x.col(0);
            double tcp_sum = 0.0, tdp_sum = 0.0;
            for (int i = 0; i < N; ++i) {
                if (std::abs(lambda(i)) < 1e-20) continue;
                Vector2d xi_i = state.foil.x.col(i);
                double a_src = PanelConstSourceStream(XjTE, xi_i);
                tcp_sum += lambda(i) * a_src;
                double a_v_, b_v_;
                PanelLinVortexStream(XjTE, xi_i, a_v_, b_v_);
                tdp_sum += lambda(i) * (a_v_ + b_v_);
            }
            tcp_sum *= 0.5 * (v(N - 1) - v(0));
            tdp_sum *= 0.5 * (v(0) - v(N - 1));

            int te_nodes[] = {0, 1, N - 2, N - 1};
            for (int kk : te_nodes) {
                if (kk < 0 || kk >= N) continue;
                for (int dd = 0; dd < 2; ++dd) {
                    Vector2d v1 = state.foil.x.col(0) - state.foil.x.col(1);
                    Vector2d v2 = state.foil.x.col(N - 1) - state.foil.x.col(N - 2);
                    double v1n = v1.norm(), v2n = v2.norm();
                    if (v1n < 1e-10 || v2n < 1e-10) continue;
                    Vector2d t1 = v1 / v1n, t2 = v2 / v2n;
                    Vector2d traw = 0.5 * (t1 + t2);
                    double tnorm = traw.norm();
                    if (tnorm < 1e-10) continue;
                    Vector2d tvec = traw / tnorm;
                    Vector2d pvec = s / sn;
                    Matrix2d dt1_dv1 = (Matrix2d::Identity() - t1 * t1.transpose()) / v1n;
                    Vector2d dv1 = Vector2d::Zero();
                    if (kk == 0) dv1(dd) = 1.0;
                    else if (kk == 1) dv1(dd) = -1.0;
                    Vector2d dt1_d = dt1_dv1 * dv1;
                    Matrix2d dt2_dv2 = (Matrix2d::Identity() - t2 * t2.transpose()) / v2n;
                    Vector2d dv2 = Vector2d::Zero();
                    if (kk == N - 1) dv2(dd) = 1.0;
                    else if (kk == N - 2) dv2(dd) = -1.0;
                    Vector2d dt2_d = dt2_dv2 * dv2;
                    Vector2d dtraw_d = 0.5 * (dt1_d + dt2_d);
                    Matrix2d dt_dtraw = (Matrix2d::Identity() - tvec * tvec.transpose()) / tnorm;
                    Vector2d dt_d = dt_dtraw * dtraw_d;
                    Vector2d ds_d = Vector2d::Zero();
                    if (kk == N - 1) ds_d(dd) = 1.0;
                    else if (kk == 0) ds_d(dd) = -1.0;
                    Matrix2d dp_ds = (Matrix2d::Identity() - pvec * pvec.transpose()) / sn;
                    Vector2d dp_d = dp_ds * ds_d;
                    double cross_val = tvec(0) * pvec(1) - tvec(1) * pvec(0);
                    double sign_c = (cross_val >= 0) ? 1.0 : -1.0;
                    double dcross = dt_d(0) * pvec(1) + tvec(0) * dp_d(1) -
                                    dt_d(1) * pvec(0) - tvec(1) * dp_d(0);
                    double dtcp_d = sign_c * dcross;
                    double dtdp_d = dt_d.dot(pvec) + tvec.dot(dp_d);
                    sens(2 * kk + dd) += dtcp_d * tcp_sum + dtdp_d * tdp_sum;
                }
            }
        }
    }

    return sens;
}

}  // namespace cppfoil
