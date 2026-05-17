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
#include "cppfoil/core/adjoint.h"

#include <cmath>

#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/viscous_geom_deriv.h"
#include "cppfoil/core/viscous_residual_deriv.h"
#include "cppfoil/core/inviscid_deriv.h"
#include "cppfoil/core/panel_influence.h"
#include "cppfoil/core/panel_influence_deriv.h"
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

// ============================================================================
// TE geometry derivatives: d(tcp)/d(x[k]) and d(tdp)/d(x[k])
// tcp and tdp depend on nodes 0, 1, N-2, N-1 through:
//   t1 = normalize(x[0] - x[1]), t2 = normalize(x[N-1] - x[N-2])
//   t = normalize(0.5*(t1+t2))
//   s = x[N-1] - x[0], p = s/|s|
//   tcp = |t × p|, tdp = t · p
//
// Returns d(tcp)/d(x[k]_dim) and d(tdp)/d(x[k]_dim) for a given node k and dim.
// Only nonzero for k ∈ {0, 1, N-2, N-1}.
// ============================================================================
static void GetTEInfoDeriv(const MatrixXd& x, int k, int dim,
                            double& dtcp, double& dtdp) {
    int N = static_cast<int>(x.cols());
    dtcp = 0.0;
    dtdp = 0.0;

    // Only nodes 0, 1, N-2, N-1 affect tcp/tdp
    if (k != 0 && k != 1 && k != N - 2 && k != N - 1) return;

    // --- Forward computation (mirrors GetTEInfo) ---
    Vector2d v1 = x.col(0) - x.col(1);
    Vector2d v2 = x.col(N - 1) - x.col(N - 2);
    double v1n = v1.norm();
    double v2n = v2.norm();
    if (v1n < 1e-10 || v2n < 1e-10) return;

    Vector2d t1 = v1 / v1n;
    Vector2d t2 = v2 / v2n;

    Vector2d t_raw = 0.5 * (t1 + t2);
    double tn = t_raw.norm();
    if (tn < 1e-10) return;
    Vector2d t = t_raw / tn;

    Vector2d s = x.col(N - 1) - x.col(0);
    double sn = s.norm();
    if (sn < 1e-10) return;
    Vector2d p = s / sn;

    double cross = t(0) * p(1) - t(1) * p(0);
    double tcp_val = std::abs(cross);
    double sign_cross = (cross >= 0) ? 1.0 : -1.0;
    // double tdp_val = t.dot(p);

    // --- d(t1)/d(x[k]_dim) ---
    // t1 = v1/|v1|, v1 = x[0]-x[1]
    // d(v1)/d(x[0]_dim) = e_dim, d(v1)/d(x[1]_dim) = -e_dim
    // d(t1)/d(v1) = (I - t1*t1^T)/|v1|
    Matrix2d dt1_dv1 = (Matrix2d::Identity() - t1 * t1.transpose()) / v1n;
    Vector2d dv1(0, 0);
    if (k == 0) dv1(dim) = 1.0;
    else if (k == 1) dv1(dim) = -1.0;
    Vector2d dt1 = dt1_dv1 * dv1;  // d(t1)/d(x[k]_dim)

    // --- d(t2)/d(x[k]_dim) ---
    Matrix2d dt2_dv2 = (Matrix2d::Identity() - t2 * t2.transpose()) / v2n;
    Vector2d dv2(0, 0);
    if (k == N - 1) dv2(dim) = 1.0;
    else if (k == N - 2) dv2(dim) = -1.0;
    Vector2d dt2 = dt2_dv2 * dv2;

    // --- d(t)/d(x[k]_dim) ---
    // t = normalize(0.5*(t1+t2))
    Vector2d dt_raw = 0.5 * (dt1 + dt2);
    Matrix2d dt_dtraw = (Matrix2d::Identity() - t * t.transpose()) / tn;
    Vector2d dt = dt_dtraw * dt_raw;

    // --- d(p)/d(x[k]_dim) ---
    // p = s/|s|, s = x[N-1]-x[0]
    Vector2d ds(0, 0);
    if (k == N - 1) ds(dim) = 1.0;
    else if (k == 0) ds(dim) = -1.0;
    Matrix2d dp_ds = (Matrix2d::Identity() - p * p.transpose()) / sn;
    Vector2d dp = dp_ds * ds;

    // --- d(tcp)/d(x[k]_dim) ---
    // tcp = |cross| where cross = t(0)*p(1) - t(1)*p(0)
    double dcross = dt(0) * p(1) + t(0) * dp(1) - dt(1) * p(0) - t(1) * dp(0);
    dtcp = sign_cross * dcross;

    // --- d(tdp)/d(x[k]_dim) ---
    // tdp = t·p
    dtdp = dt.dot(p) + t.dot(dp);
}








// ============================================================================
// CalcForceLin: compute partial derivatives of cl, cm, cdpi w.r.t. foil.x
// holding ue (gamma) fixed. Also computes cm_ue and cdpi_ue linearizations.
// ============================================================================
void CalcForceLin(FoilState& state) {
    double chord = state.geom.chord;
    Vector2d xref = state.geom.xref;
    double Vinf = state.param.Vinf > 0 ? state.param.Vinf : 1.0;
    double alpha = state.oper.alpha;
    int N = state.foil.N;

    // Get edge velocity (use inviscid ue = gamma for inviscid case)
    VectorXd ue;
    if (state.oper.viscous && state.glob.U.cols() > 0) {
        ue = state.glob.U.row(3).transpose();
    } else {
        ue = GetUeInv(state);
    }

    // Compute cp at each node (same as CalcForce)
    VectorXd cp(N);
    for (int i = 0; i < N; ++i) {
        double cp_inc = 1.0 - (ue(i) / Vinf) * (ue(i) / Vinf);
        if (state.param.KTl > 0) {
            double denom = state.param.KTb + state.param.KTl * (1.0 + state.param.KTb) * 0.5 * cp_inc;
            if (std::abs(denom) > 1e-10) {
                cp(i) = cp_inc / denom;
            } else {
                cp(i) = cp_inc;
            }
        } else {
            cp(i) = cp_inc;
        }
    }

    // Initialize partial derivative matrices
    state.post.dcl_dx_partial = MatrixXd::Zero(2, N);
    state.post.dcm_dx_partial = MatrixXd::Zero(2, N);
    state.post.dcdpi_dx_partial = MatrixXd::Zero(2, N);
    state.post.cm_ue = VectorXd::Zero(N);
    state.post.cdpi_ue = VectorXd::Zero(N);

    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    // The force integration loops over panel pairs (ip, i) matching CalcForce
    // For each panel:
    //   dx = -(x2(0)-x1(0))*cosa - (x2(1)-x1(1))*sina   (lift direction)
    //   dz = (x2(1)-x1(1))*cosa - (x2(0)-x1(0))*sina    (drag direction)
    //   cl += dx * cpbar / chord
    //   cm += moment terms / chord^2
    //   cdpi += dz * cpbar / chord
    //
    // Partial w.r.t. coordinates (holding cp fixed):
    //   dcl/dx1_x = cpbar * cosa / chord  (from d(-dxv(0))/dx1_x = 1 => no, dxv=x2-x1)
    //   dxv = x2 - x1 => d(dxv)/d(x1) = [-1, 0; 0, -1], d(dxv)/d(x2) = [1, 0; 0, 1]
    //   dx = -dxv(0)*cosa - dxv(1)*sina
    //   d(dx)/d(x1_x) = cosa, d(dx)/d(x1_z) = sina
    //   d(dx)/d(x2_x) = -cosa, d(dx)/d(x2_z) = -sina

    for (int i0 = 1; i0 <= N; ++i0) {
        int i = (i0 == N) ? 0 : i0;
        int ip = (i0 == N) ? N - 1 : i0 - 1;

        Vector2d x1 = state.foil.x.col(ip);
        Vector2d x2 = state.foil.x.col(i);
        Vector2d dxv = x2 - x1;
        Vector2d dx1 = x1 - xref;
        Vector2d dx2 = x2 - xref;

        double dz = dxv(1) * cosa - dxv(0) * sina;

        double cp1 = cp(ip);
        double cp2 = cp(i);
        double cpbar = 0.5 * (cp1 + cp2);

        // --- cl partial w.r.t. coordinates ---
        // dcl/d(x1_x) = cpbar * cosa / chord  (from ddx/dx1_x = cosa)
        // dcl/d(x1_z) = cpbar * sina / chord
        // dcl/d(x2_x) = cpbar * (-cosa) / chord
        // dcl/d(x2_z) = cpbar * (-sina) / chord
        state.post.dcl_dx_partial(0, ip) += cpbar * cosa / chord;
        state.post.dcl_dx_partial(1, ip) += cpbar * sina / chord;
        state.post.dcl_dx_partial(0, i) += cpbar * (-cosa) / chord;
        state.post.dcl_dx_partial(1, i) += cpbar * (-sina) / chord;

        // --- cdpi partial w.r.t. coordinates ---
        // dz = dxv(1)*cosa - dxv(0)*sina
        // d(dz)/d(x1_x) = sina, d(dz)/d(x1_z) = -cosa
        // d(dz)/d(x2_x) = -sina, d(dz)/d(x2_z) = cosa
        state.post.dcdpi_dx_partial(0, ip) += cpbar * sina / chord;
        state.post.dcdpi_dx_partial(1, ip) += cpbar * (-cosa) / chord;
        state.post.dcdpi_dx_partial(0, i) += cpbar * (-sina) / chord;
        state.post.dcdpi_dx_partial(1, i) += cpbar * cosa / chord;

        // --- cm partial w.r.t. coordinates ---
        // cm += (cp1 * dx1nds/3 + cp1 * dx2nds/6 + cp2 * dx1nds/6 + cp2 * dx2nds/3) / chord^2
        // where dx1nds = dxv(0)*dx1(0) + dxv(1)*dx1(1)
        //       dx2nds = dxv(0)*dx2(0) + dxv(1)*dx2(1)
        //
        // This depends on both dxv and dx1, dx2. The derivatives are:
        // d(dx1nds)/d(x1_x) = -dx1(0) + dxv(0)*(-1) = -dx1(0) - dxv(0)  [d(dxv)/dx1_x = -1, d(dx1)/dx1_x = 1... wait]
        // dxv = x2 - x1, dx1 = x1 - xref, dx2 = x2 - xref
        // d(dxv(0))/d(x1_x) = -1, d(dx1(0))/d(x1_x) = 1
        // dx1nds = dxv(0)*dx1(0) + dxv(1)*dx1(1)
        // d(dx1nds)/d(x1_x) = (-1)*dx1(0) + dxv(0)*(1) = dxv(0) - dx1(0)
        double inv_chord2 = 1.0 / (chord * chord);

        // d(dx1nds)/d(x1_x) = dxv(0) - dx1(0)  ... actually let me be more careful
        // dx1nds = dxv . dx1
        // d(dx1nds)/d(x1_x) = d(dxv)/d(x1_x) . dx1 + dxv . d(dx1)/d(x1_x)
        //                   = (-1,0) . dx1 + dxv . (1,0) = -dx1(0) + dxv(0)
        double ddx1nds_dx1x = -dx1(0) + dxv(0);
        double ddx1nds_dx1z = -dx1(1) + dxv(1);
        double ddx1nds_dx2x = dx1(0);  // d(dxv)/d(x2_x) = (1,0), d(dx1)/d(x2_x) = 0
        double ddx1nds_dx2z = dx1(1);

        // d(dx2nds)/d(x1_x) = (-1,0).dx2 + dxv.(0,0) = -dx2(0)
        double ddx2nds_dx1x = -dx2(0);
        double ddx2nds_dx1z = -dx2(1);
        // d(dx2nds)/d(x2_x) = (1,0).dx2 + dxv.(1,0) = dx2(0) + dxv(0)
        double ddx2nds_dx2x = dx2(0) + dxv(0);
        double ddx2nds_dx2z = dx2(1) + dxv(1);

        double w1 = cp1 / 3.0 + cp2 / 6.0;
        double w2 = cp1 / 6.0 + cp2 / 3.0;

        state.post.dcm_dx_partial(0, ip) += (w1 * ddx1nds_dx1x + w2 * ddx2nds_dx1x) * inv_chord2;
        state.post.dcm_dx_partial(1, ip) += (w1 * ddx1nds_dx1z + w2 * ddx2nds_dx1z) * inv_chord2;
        state.post.dcm_dx_partial(0, i) += (w1 * ddx1nds_dx2x + w2 * ddx2nds_dx2x) * inv_chord2;
        state.post.dcm_dx_partial(1, i) += (w1 * ddx1nds_dx2z + w2 * ddx2nds_dx2z) * inv_chord2;

        // --- cm_ue linearization ---
        // cm depends on cp which depends on ue: dcp/due = -2*ue/(Vinf^2)
        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);

        auto [cp_val_ip, cp_U_ip] = GetCpFromUe(ue(ip), state.param);
        auto [cp_val_i, cp_U_i] = GetCpFromUe(ue(i), state.param);
        double dcp_due_ip = cp_U_ip(3);
        double dcp_due_i = cp_U_i(3);

        state.post.cm_ue(ip) += (dcp_due_ip / 3.0 * dx1nds + dcp_due_ip / 6.0 * dx2nds) * inv_chord2;
        state.post.cm_ue(i) += (dcp_due_i / 6.0 * dx1nds + dcp_due_i / 3.0 * dx2nds) * inv_chord2;

        // --- cdpi_ue linearization ---
        state.post.cdpi_ue(ip) += dz * 0.5 * dcp_due_ip / chord;
        state.post.cdpi_ue(i) += dz * 0.5 * dcp_due_i / chord;
    }

    // ========================================================================
    // Viscous output partials: dJ/dV at fixed foil.x, V = U flattened [4*Nsys].
    // - cl, cm, cdpi: only depend on U[3, i] for i ∈ [0, N) (foil edge velocity)
    //   so cl_V[4*i + 3] = cl_ue[i] (similar cm, cdpi). All other entries zero.
    // - cd: Squire-Young from end-of-wake state Uw only; cd_V is sparse at iw.
    // - dcd/dx_partial = 0 because Squire-Young has no direct foil.x dependence.
    // ========================================================================
    if (state.oper.viscous && state.glob.U.cols() > 0) {
        int Nsys = static_cast<int>(state.glob.U.cols());
        state.post.cl_V = VectorXd::Zero(4 * Nsys);
        state.post.cm_V = VectorXd::Zero(4 * Nsys);
        state.post.cd_V = VectorXd::Zero(4 * Nsys);
        for (int i = 0; i < N; ++i) {
            state.post.cl_V(4 * i + 3) = state.post.cl_ue(i);
            state.post.cm_V(4 * i + 3) = state.post.cm_ue(i);
        }

        // Squire-Young cd from end-of-wake state.
        state.post.cd_U.setZero();
        state.post.dcd_dx_partial = MatrixXd::Zero(2, N);
        if (state.vsol.Is.size() > 2 && !state.vsol.Is[2].empty()) {
            int iw = state.vsol.Is[2].back();
            if (iw < state.glob.U.cols()) {
                Vector4d Uw = state.glob.U.col(iw);
                auto [H, H_U] = GetH(Uw);
                auto [uk, uk_U] = GetUk(Uw(3), state.param);
                double q = uk / Vinf;
                double p = 0.5 * (5.0 + H);
                double qp = std::pow(q, p);
                double lnq = std::log(std::max(q, 1e-300));

                // cd = 2 * Uw[0] * q^p
                // dcd/dUw[k] = 2*δ_k0 * q^p + 2*Uw[0] * d(q^p)/dUw[k]
                // d(q^p)/dUw[k] = q^p * (dp/dUw[k] * lnq + (p/q) * dq/dUw[k])
                //               = q^p * (0.5*H_U[k]*lnq + (p/q)*(1/Vinf)*uk_U[k])
                for (int k = 0; k < 4; ++k) {
                    double dqp = qp * (0.5 * H_U(k) * lnq +
                                       (p / q) * (1.0 / Vinf) * uk_U(k));
                    state.post.cd_U(k) = 2.0 * Uw(0) * dqp;
                    if (k == 0) state.post.cd_U(k) += 2.0 * qp;
                    state.post.cd_V(4 * iw + k) = state.post.cd_U(k);
                }
            }
        }
    }
}

// ============================================================================
// InviscidSensitivities: compute dJ/dx using the adjoint method
// ============================================================================
Sensitivity InviscidSensitivities(FoilState& state) {
    int N = state.foil.N;
    if (N <= 0 || state.isol.AIC.rows() == 0) {
        throw std::runtime_error("InviscidSensitivities: no inviscid solution");
    }

    // Step 1: Compute force partial derivatives
    CalcForceLin(state);

    double alpha = state.oper.alpha;

    // Step 2: Form dJ/dgamma for each objective
    // For inviscid flow, ue(i) = sgnue(i) * gamma(i)
    // cl_ue is dcl/due [N], and due/dgamma(i) = sgnue(i)
    // So dcl/dgamma = cl_ue .* sgnue (elementwise)
    // Similarly for cm, cdpi

    // The AIC system is [N+1 x N+1] (includes streamfunction unknown and Kutta)
    // gamma is [N], extended to [N+1] with streamfunction as last entry
    // dJ/dgamma_extended = [dJ/dgamma; 0] (streamfunction doesn't affect forces directly)

    VectorXd dcl_dgam = VectorXd::Zero(N + 1);
    VectorXd dcm_dgam = VectorXd::Zero(N + 1);
    VectorXd dcdpi_dgam = VectorXd::Zero(N + 1);

    for (int i = 0; i < N; ++i) {
        double sg = state.isol.sgnue(i);
        dcl_dgam(i) = state.post.cl_ue(i) * sg;
        dcm_dgam(i) = state.post.cm_ue(i) * sg;
        dcdpi_dgam(i) = state.post.cdpi_ue(i) * sg;
    }

    // Step 3: Solve adjoint equations: AIC^T * lambda = dJ/dgamma
    // Note: the residual is R = A*gamma - rhs = 0
    // Adjoint: A^T * lambda = (dJ/dgamma)^T
    // But since R = A*gam - rhs, dR/dgamma = A, so adjoint is A^T * lambda = dJ/dgamma.
    // The sign convention: dJ/dx_total = dJ/dx_partial - lambda^T * dR/dx
    // Since dR = d(A*gam - rhs), we need to subtract.
    // Actually, for R = A*gam - rhs = 0, implicit function theorem gives:
    //   dgam/dx = -A^{-1} * dR/dx  (where dR/dx = dA/dx*gam - drhs/dx)
    //   dJ/dx = dJ/dx_partial + dJ/dgam * dgam/dx
    //         = dJ/dx_partial - dJ/dgam * A^{-1} * dR/dx
    //         = dJ/dx_partial - lambda^T * dR/dx
    // where A^T * lambda = dJ/dgam
    // So: dJ/dx_total = dJ/dx_partial - lambda^T * dR/dx_foil

    const MatrixXd& AIC = state.isol.AIC;
    MatrixXd AIC_T = AIC.transpose();

    // Solve for adjoint vectors
    auto lu = AIC_T.partialPivLu();

    VectorXd lambda_cl = lu.solve(dcl_dgam);
    VectorXd lambda_cm = lu.solve(dcm_dgam);
    VectorXd lambda_cdpi = lu.solve(dcdpi_dgam);

    // Step 4: Compute VJP: lambda^T * dR/dx_foil
    VectorXd vjp_cl = AICResidualSensVJP(state, lambda_cl, alpha);
    VectorXd vjp_cm = AICResidualSensVJP(state, lambda_cm, alpha);
    VectorXd vjp_cdpi = AICResidualSensVJP(state, lambda_cdpi, alpha);

    // Step 5: Total sensitivity = partial - VJP
    Sensitivity sens;
    sens.dcl_dx = MatrixXd::Zero(2, N);
    sens.dcm_dx = MatrixXd::Zero(2, N);
    sens.dcdpi_dx = MatrixXd::Zero(2, N);
    sens.dcd_dx = MatrixXd::Zero(2, N);  // Zero for inviscid (no viscous drag)

    for (int k = 0; k < N; ++k) {
        sens.dcl_dx(0, k) = state.post.dcl_dx_partial(0, k) - vjp_cl(2 * k);
        sens.dcl_dx(1, k) = state.post.dcl_dx_partial(1, k) - vjp_cl(2 * k + 1);
        sens.dcm_dx(0, k) = state.post.dcm_dx_partial(0, k) - vjp_cm(2 * k);
        sens.dcm_dx(1, k) = state.post.dcm_dx_partial(1, k) - vjp_cm(2 * k + 1);
        sens.dcdpi_dx(0, k) = state.post.dcdpi_dx_partial(0, k) - vjp_cdpi(2 * k);
        sens.dcdpi_dx(1, k) = state.post.dcdpi_dx_partial(1, k) - vjp_cdpi(2 * k + 1);
    }

    return sens;
}


void AlphaSensitivities(FoilState& state, Sensitivity& sens) {
    double alpha = state.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    // dgam/dalpha = gamref[:,0]*(-sin a) + gamref[:,1]*(cos a), times deg2rad
    VectorXd dgam_dalpha = (-sina * state.isol.gamref.col(0) +
                             cosa * state.isol.gamref.col(1)) * kDegToRad;

    // dueinv/dalpha = sgnue .* dgam/dalpha
    VectorXd dueinv_dalpha = state.isol.sgnue.cwiseProduct(dgam_dalpha);

    sens.dcl_dalpha = state.post.cl_alpha +
                      state.post.cl_ue.dot(dueinv_dalpha);
    sens.dcm_dalpha = 0.0 + state.post.cm_ue.dot(dueinv_dalpha);
    sens.dcd_dalpha = state.post.cdpi_alpha +
                      state.post.cdpi_ue.dot(dueinv_dalpha);
}

static void ComputeCpDerivsMa(double ue, const Param& param, double Ma,
                               double& dcp_dMa, double& d2cp_duedMa) {
    double Vinf = std::max(param.Vinf, 1e-10);
    double cp_inc = 1.0 - (ue / Vinf) * (ue / Vinf);
    double dcp_due_inc = -2.0 * ue / (Vinf * Vinf);

    if (Ma <= 0 || param.KTl <= 0) {
        dcp_dMa = 0.0;
        d2cp_duedMa = 0.0;
        return;
    }

    double b = param.KTb;  // sqrt(1-Ma^2)
    double l = param.KTl;  // Ma^2/(1+b)^2

    // db/dMa = -Ma/b
    double db = -Ma / b;
    // dl/dMa = 2*Ma/(1+b)^2 - 2*Ma^2/(1+b)^3 * db
    //        = 2*Ma/(1+b)^2 + 2*Ma^3/((1+b)^3 * b)
    double opb = 1.0 + b;
    double dl = 2.0 * Ma / (opb * opb) + 2.0 * Ma * Ma * Ma / (opb * opb * opb * b);

    double den = b + 0.5 * l * (1.0 + b) * cp_inc;
    double den_cp = 0.5 * l * (1.0 + b);

    // cp = cp_inc / den
    double cp = cp_inc / den;

    // dden/dMa = db + 0.5*(dl*(1+b) + l*db)*cp_inc
    double dden_dMa = db + 0.5 * (dl * opb + l * db) * cp_inc;
    // dden_cp/dMa = 0.5*(dl*(1+b) + l*db)
    double dden_cp_dMa = 0.5 * (dl * opb + l * db);

    // dcp/dMa = -cp_inc / den^2 * dden_dMa  (cp_inc doesn't depend on Ma)
    dcp_dMa = -cp / den * dden_dMa;

    // d(dcp/due)/dMa: dcp/due = dcp_due_inc * (1 - cp*den_cp) / den
    // d(dcp/due)/dMa = dcp_due_inc * d/dMa[(1 - cp*den_cp)/den]
    double dfac_dMa = (-(dcp_dMa * den_cp + cp * dden_cp_dMa) * den -
                       (1.0 - cp * den_cp) * dden_dMa) / (den * den);
    d2cp_duedMa = dcp_due_inc * dfac_dMa;
}

// ============================================================================
// MachSensitivities: compute dJ/dMa for cl, cd, cm

void MachSensitivities(FoilState& state, Sensitivity& sens) {
    int N = state.foil.N;
    double Ma = state.oper.Ma;
    double chord = state.geom.chord;
    double alpha = state.oper.alpha;

    if (Ma <= 0) {
        sens.dcl_dMa = 0.0;
        sens.dcd_dMa = 0.0;
        sens.dcm_dMa = 0.0;
        return;
    }

    VectorXd ue = GetUeInv(state);

    double dcl = 0.0, dcm = 0.0, dcdpi = 0.0;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);
    Vector2d xref = state.geom.xref;

    for (int i0 = 1; i0 <= N; ++i0) {
        int i = (i0 == N) ? 0 : i0;
        int ip = (i0 == N) ? N - 1 : i0 - 1;

        Vector2d x1 = state.foil.x.col(ip);
        Vector2d x2 = state.foil.x.col(i);
        Vector2d dxv = x2 - x1;
        Vector2d dx1 = x1 - xref;
        Vector2d dx2 = x2 - xref;

        double dx_lift = -dxv(0) * cosa - dxv(1) * sina;
        double dz = dxv(1) * cosa - dxv(0) * sina;

        double dcp_ip, d2cp_ip;
        double dcp_i, d2cp_i;
        ComputeCpDerivsMa(ue(ip), state.param, Ma, dcp_ip, d2cp_ip);
        ComputeCpDerivsMa(ue(i), state.param, Ma, dcp_i, d2cp_i);

        double dcpbar = 0.5 * (dcp_ip + dcp_i);
        dcl += dx_lift * dcpbar;
        dcdpi += dz * dcpbar;

        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
        dcm += dcp_ip / 3.0 * dx1nds + dcp_ip / 6.0 * dx2nds +
               dcp_i / 6.0 * dx1nds + dcp_i / 3.0 * dx2nds;
    }

    sens.dcl_dMa = dcl / chord;
    sens.dcm_dMa = dcm / (chord * chord);
    sens.dcd_dMa = dcdpi / chord;
}

// ============================================================================
// ViscousSensitivities — adjoint-based dJ/dx_foil for J ∈ {cl, cd, cm}.
//
// Steps:
//   1. CalcForceLin → populates cl_V, cm_V, cd_V (∂J/∂V) and ∂J/∂x_partial.
//   2. StoreAdjointFactorization → R_V LU at the converged state.
//   3. For each J: λ_J = R_V^{-T} ∂J/∂V via SolveViscousAdjoint.
//   4. VJP λ_J^T · ∂R_V/∂x_foil
//        = λ_bl^T · DRblDFoilX(state, dfoil_s_dx)             [3*Nsys x 2N]
//        + ComputeRueVJP(state, λ_ue, dfoil_s_dx)             [2N]
//   5. dJ/dx_foil = ∂J/∂x_partial - vjp_J.
// ============================================================================
Sensitivity ViscousSensitivities(FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;

    if (!state.oper.viscous || state.glob.U.cols() == 0) {
        throw std::runtime_error(
            "ViscousSensitivities: viscous primal not converged");
    }

    // 1. Force partials
    CalcForceLin(state);

    // 2. Adjoint factorization
    StoreAdjointFactorization(state);
    if (!state.glob.lu_factorization) {
        throw std::runtime_error(
            "ViscousSensitivities: adjoint LU factorization failed");
    }

    // 3. Solve λ_J = R_V^{-T} ∂J/∂V for J ∈ {cl, cm, cd}
    VectorXd lambda_cl = SolveViscousAdjoint(state, state.post.cl_V);
    VectorXd lambda_cm = SolveViscousAdjoint(state, state.post.cm_V);
    VectorXd lambda_cd = SolveViscousAdjoint(state, state.post.cd_V);

    // 4. VJPs — DRbl as dense matrix * lambda_bl (lambda_bl size 3*Nsys),
    //    plus ComputeRueVJP for the wake-side chain.
    MatrixXd dfoil_s_dx = DFoilSDFoilXCumulative(state.foil);
    MatrixXd dRbl_dx = DRblDFoilX(state, dfoil_s_dx);  // [3*Nsys x 2N]

    auto vjp = [&](const VectorXd& lambda) -> VectorXd {
        VectorXd lambda_bl = lambda.head(3 * Nsys);
        VectorXd lambda_ue = lambda.tail(Nsys);
        VectorXd out = (lambda_bl.transpose() * dRbl_dx).transpose();  // [2N]
        out += ComputeRueVJP(state, lambda_ue, dfoil_s_dx);
        return out;
    };
    VectorXd vjp_cl = vjp(lambda_cl);
    VectorXd vjp_cm = vjp(lambda_cm);
    VectorXd vjp_cd = vjp(lambda_cd);

    // 5. Total dJ/dx = partial - vjp
    Sensitivity sens;
    sens.dcl_dx = MatrixXd::Zero(2, N);
    sens.dcm_dx = MatrixXd::Zero(2, N);
    sens.dcd_dx = MatrixXd::Zero(2, N);
    sens.dcdpi_dx = MatrixXd::Zero(2, N);
    for (int k = 0; k < N; ++k) {
        sens.dcl_dx(0, k) = state.post.dcl_dx_partial(0, k) - vjp_cl(2 * k);
        sens.dcl_dx(1, k) = state.post.dcl_dx_partial(1, k) - vjp_cl(2 * k + 1);
        sens.dcm_dx(0, k) = state.post.dcm_dx_partial(0, k) - vjp_cm(2 * k);
        sens.dcm_dx(1, k) = state.post.dcm_dx_partial(1, k) - vjp_cm(2 * k + 1);
        // dcd/dx_partial is zero for Squire-Young, so this is just -vjp_cd.
        sens.dcd_dx(0, k) = state.post.dcd_dx_partial(0, k) - vjp_cd(2 * k);
        sens.dcd_dx(1, k) = state.post.dcd_dx_partial(1, k) - vjp_cd(2 * k + 1);
    }

    // 6. Scalar parameter chains: dJ/d{alpha, Ma, Re} via adjoint VJP.
    //    dJ/dp = ∂J/∂p_partial - λ^T · ∂R_V/∂p.
    //
    // alpha: ∂cl/∂alpha = cl_alpha; ∂cm/∂alpha currently 0; ∂cd/∂alpha = 0
    //        (Squire-Young uses Uw only, no direct alpha).
    {
        VectorXd dRV_dalpha = DRVDAlpha(state);
        sens.dcl_dalpha = state.post.cl_alpha - lambda_cl.dot(dRV_dalpha);
        sens.dcm_dalpha = 0.0 - lambda_cm.dot(dRV_dalpha);
        sens.dcd_dalpha = 0.0 - lambda_cd.dot(dRV_dalpha);
    }

    // Ma: ∂J/∂Ma_partial via dcp/dMa chain at fixed U[3] (= ue). Reuse the
    // inviscid MachSensitivities formula: it sums cp_Ma over panels for
    // each force coefficient. For viscous cd (Squire-Young), ∂cd/∂Ma chains
    // through uk(Ma) and H_U·(uk depends on Ma): captured by AlphaCalcMa()
    // helper below. For now use a quick partial chain via FD on Ma at fixed
    // U for ∂J/∂Ma_partial.
    {
        VectorXd dRV_dMa = DRVDMa(state);
        // ∂J/∂Ma_partial via FD on Ma at fixed U: re-evaluate cl/cd/cm with
        // perturbed Minf (and corresponding KTb/KTl), holding U/foil.x fixed.
        double Ma0 = state.oper.Ma;
        // FD on Ma with FULL InitThermo update (KTb, KTl, H0, cps, mu0, rho0).
        // Updating only KTb/KTl was wrong because Squire-Young cd uses uk
        // which depends on H0 (compressibility correction in GetMach2/GetUk).
        auto eval_J = [&](double Ma) -> Vector3d {
            FoilState sp = state;
            sp.oper.Ma = Ma;
            InitThermo(sp);
            CalcForce(sp);
            return Vector3d(sp.post.cl, sp.post.cd, sp.post.cm);
        };
        double h = std::max(std::abs(Ma0) * 1e-2, 1e-3);
        Vector3d Jp2 = eval_J(Ma0 + 2 * h);
        Vector3d Jp1 = eval_J(Ma0 + h);
        Vector3d Jm1 = eval_J(Ma0 - h);
        Vector3d Jm2 = eval_J(Ma0 - 2 * h);
        Vector3d dJ_dMa_partial = (-Jp2 + 8.0 * Jp1 - 8.0 * Jm1 + Jm2) / (12.0 * h);
        sens.dcl_dMa = dJ_dMa_partial(0) - lambda_cl.dot(dRV_dMa);
        sens.dcd_dMa = dJ_dMa_partial(1) - lambda_cd.dot(dRV_dMa);
        sens.dcm_dMa = dJ_dMa_partial(2) - lambda_cm.dot(dRV_dMa);
    }

    // Re: only chain through Ret in BL closures; no direct Re in cl/cd/cm.
    {
        VectorXd dRV_dRe = DRVDRe(state);
        sens.dcl_dRe = -lambda_cl.dot(dRV_dRe);
        sens.dcd_dRe = -lambda_cd.dot(dRV_dRe);
        sens.dcm_dRe = -lambda_cm.dot(dRV_dRe);
    }

    return sens;
}

// ============================================================================
// DRVDAlpha — analytical ∂R_V/∂alpha as a [4*Nsys] vector.
// Matches the same chain that SolveGlob uses in its constrained-cl Newton
// block (coupled_solver.cpp:529-543).
// ============================================================================
VectorXd DRVDAlpha(const FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    double alpha = state.oper.alpha;
    double sa = Sind(alpha);
    double ca = Cosd(alpha);

    VectorXd dR = VectorXd::Zero(4 * Nsys);
    for (int i = 0; i < Nsys; ++i) {
        double dueinv_dalpha = 0.0;
        if (i < N) {
            // Foil: ueinv[i] = sgnue[i] * gamma[i]
            double dgam = -state.isol.gamref(i, 0) * sa +
                          state.isol.gamref(i, 1) * ca;
            dueinv_dalpha = state.isol.sgnue(i) * dgam;
        } else {
            int iw = i - N;
            if (iw == 0) {
                // Wake station 0 uses continuity BC ueinv[N+0] = uea[N-1]
                // = sgnue[N-1] * gamma[N-1]. The uewiref[0,:] velocity in
                // state.isol.uewi(0) is OVERWRITTEN by this in get_ueinv,
                // so the residual chain goes through gamma, not uewiref.
                double dgam = -state.isol.gamref(N - 1, 0) * sa +
                              state.isol.gamref(N - 1, 1) * ca;
                dueinv_dalpha = state.isol.sgnue(N - 1) * dgam;
            } else if (iw < static_cast<int>(state.isol.uewiref.rows())) {
                // Wake j>0: uewi[j] = uewiref[j,0]*cos α + uewiref[j,1]*sin α
                dueinv_dalpha = -state.isol.uewiref(iw, 0) * sa +
                                state.isol.uewiref(iw, 1) * ca;
            }
        }
        dueinv_dalpha *= kDegToRad;
        // R_ue[i] = U[3,i] - ueinv[i] - ..., so ∂R_ue/∂alpha = -∂ueinv/∂alpha
        dR(3 * Nsys + i) = -dueinv_dalpha;
    }
    return dR;
}

// ============================================================================
// DRVDScalar — central FD of R_V w.r.t. a scalar parameter (Ma or mu0).
//
// Holds U, foil.x, wake.x, ue_sigma fixed at converged state and rebuilds
// R_BL via BuildGlobSys (which re-evaluates BL closures with perturbed
// param) plus the R_ue residual analytically. Returns the [4*Nsys] FD.
// ============================================================================
static VectorXd EvaluateRV(FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    BuildGlobSys(state);
    VectorXd R = VectorXd::Zero(4 * Nsys);
    R.head(3 * Nsys) = state.glob.R;
    // R_ue[i] = U[3,i] - ueinv[i] - ue_m[i,:] * (ds.*ue)
    VectorXd mass(Nsys);
    for (int j = 0; j < Nsys; ++j) {
        mass(j) = state.glob.U(1, j) * state.glob.U(3, j);
    }
    VectorXd src = state.vsol.ue_m * mass;
    for (int i = 0; i < Nsys; ++i) {
        double ueinv_i;
        if (i < N) {
            ueinv_i = state.isol.sgnue(i) * state.isol.gam(i);
        } else if (i - N == 0) {
            ueinv_i = state.isol.sgnue(N - 1) * state.isol.gam(N - 1);
        } else {
            ueinv_i = state.isol.uewi(i - N);
        }
        R(3 * Nsys + i) = state.glob.U(3, i) - ueinv_i - src(i);
    }
    return R;
}

// 4th-order central FD: df/dx ≈ (-f(+2h) + 8f(+h) - 8f(-h) + f(-2h)) / (12h)
// Truncation O(h^4), so the optimal h is much larger than 2nd-order central
// — pushing cancellation noise (eps/h) well below truncation while keeping h
// well clear of the parameter's natural scale. For mu0 ~ 1e-7 with Newton
// residual floor ~3e-12, optimal h ≈ mu0 * 1e-2 yields ~3e-10 rel precision.
template <typename EvalR>
static VectorXd CentralFD4(double x0, double h, EvalR&& eval) {
    VectorXd Rp2 = eval(x0 + 2 * h);
    VectorXd Rp1 = eval(x0 + h);
    VectorXd Rm1 = eval(x0 - h);
    VectorXd Rm2 = eval(x0 - 2 * h);
    return (-Rp2 + 8.0 * Rp1 - 8.0 * Rm1 + Rm2) / (12.0 * h);
}

VectorXd DRVDMa(FoilState& state, double h) {
    double Ma0 = state.oper.Ma;
    // Save ALL Ma-dependent derived params. InitThermo recomputes:
    //   Minf, KTb, KTl, H0, cps, mu0, rho0, muinf
    // when state.oper.Ma changes. Earlier only updating KTb/KTl missed H0,
    // cps, mu0, rho0 — which break Mach-dependent BL closures (GetMach2,
    // GetHss, GetRet under compressible path).
    Param param0 = state.param;
    if (h <= 0) h = std::max(std::abs(Ma0) * 1e-2, 1e-3);
    auto eval = [&](double Ma) {
        state.oper.Ma = Ma;
        InitThermo(state);
        state.glob.realloc = true;
        return EvaluateRV(state);
    };
    VectorXd dR = CentralFD4(Ma0, h, eval);
    // Restore baseline
    state.oper.Ma = Ma0;
    state.param = param0;
    state.glob.realloc = true;
    BuildGlobSys(state);
    return dR;
}

VectorXd DRVDRe(FoilState& state, double h) {
    double mu0_0 = state.param.mu0;
    double Re0 = state.oper.Re;
    // R depends on mu0 only through Ret = X/mu0. So R^(n)(mu0) ~ R/mu0^n,
    // and the 4th-order central FD optimal h is roughly mu0 * eps^(1/5).
    // For mu0 ~ 1e-7, Newton floor ~3e-12: h_opt ~ mu0 * (3e-12)^(1/5) ~
    // mu0 * 1e-2. Use that. At this h, truncation + cancellation both ~
    // 3e-10 rel on ∂R/∂mu0.
    if (h <= 0) h = std::abs(mu0_0) * 1e-2;
    double dmu0_dRe = -mu0_0 / Re0;
    auto eval = [&](double mu0) {
        state.param.mu0 = mu0;
        state.glob.realloc = true;
        return EvaluateRV(state);
    };
    VectorXd dR_dmu0 = CentralFD4(mu0_0, h, eval);
    state.param.mu0 = mu0_0;
    state.glob.realloc = true;
    BuildGlobSys(state);
    return dR_dmu0 * dmu0_dRe;
}

// ============================================================================
// SolveViscousAdjoint — transpose solve through the converged viscous R_V.
//
// R_V is assembled and LU-factored by StoreAdjointFactorization. Eigen's
// SparseLU supports transpose solves via the `.transpose()` view, which is
// what we use here.
// ============================================================================
VectorXd SolveViscousAdjoint(FoilState& state, const VectorXd& rhs) {
    int Nsys = state.glob.Nsys;
    if (rhs.size() != 4 * Nsys) {
        throw std::runtime_error(
            "SolveViscousAdjoint: rhs size must equal 4*Nsys");
    }
    if (!state.glob.lu_factorization) {
        StoreAdjointFactorization(state);
    }
    if (!state.glob.lu_factorization) {
        throw std::runtime_error(
            "SolveViscousAdjoint: R_V factorization unavailable (singular?)");
    }
    return state.glob.lu_factorization->transpose().solve(rhs);
}

}  // namespace cppfoil
