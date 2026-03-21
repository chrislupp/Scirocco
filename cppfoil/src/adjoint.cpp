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
#include <cstring>
#include <stdexcept>

#include "cppfoil/core/boundary_layer.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/inviscid_deriv.h"
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

void GetTEInfoDerivPublic(const MatrixXd& x, int k, int dim,
                           double& dtcp, double& dtdp) {
    GetTEInfoDeriv(x, k, dim, dtcp, dtdp);
}

// ============================================================================
// GetTEGapDeriv: d(hTE)/d(x[k]_dim) and d(dtdx)/d(x[k]_dim)
// Same intermediate quantities as GetTEInfoDeriv but different outputs.
// Only nonzero for k ∈ {0, 1, N-2, N-1}.
// ============================================================================
static void GetTEGapDeriv(const MatrixXd& x, int k, int dim,
                           double& dhTE, double& ddtdx) {
    int N = static_cast<int>(x.cols());
    dhTE = 0.0;
    ddtdx = 0.0;

    if (k != 0 && k != 1 && k != N - 2 && k != N - 1) return;

    // Forward computation (same intermediates as GetTEInfoDeriv)
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

    // d(t1)/d(x[k]_dim)
    Matrix2d dt1_dv1 = (Matrix2d::Identity() - t1 * t1.transpose()) / v1n;
    Vector2d dv1(0, 0);
    if (k == 0) dv1(dim) = 1.0;
    else if (k == 1) dv1(dim) = -1.0;
    Vector2d dt1 = dt1_dv1 * dv1;

    // d(t2)/d(x[k]_dim)
    Matrix2d dt2_dv2 = (Matrix2d::Identity() - t2 * t2.transpose()) / v2n;
    Vector2d dv2(0, 0);
    if (k == N - 1) dv2(dim) = 1.0;
    else if (k == N - 2) dv2(dim) = -1.0;
    Vector2d dt2 = dt2_dv2 * dv2;

    // d(t)/d(x[k]_dim)
    Vector2d dt_raw = 0.5 * (dt1 + dt2);
    Matrix2d dt_dtraw = (Matrix2d::Identity() - t * t.transpose()) / tn;
    Vector2d dt = dt_dtraw * dt_raw;

    // d(s)/d(x[k]_dim)
    Vector2d ds(0, 0);
    if (k == N - 1) ds(dim) = 1.0;
    else if (k == 0) ds(dim) = -1.0;

    // d(hTE)/d(x[k]_dim)
    // hTE = -s(0)*t(1) + s(1)*t(0)
    dhTE = -ds(0) * t(1) - s(0) * dt(1) + ds(1) * t(0) + s(1) * dt(0);

    // d(dtdx)/d(x[k]_dim)
    // dtdx = t1(0)*t2(1) - t2(0)*t1(1)
    ddtdx = dt1(0) * t2(1) + t1(0) * dt2(1) - dt2(0) * t1(1) - t2(0) * dt1(1);
}

// ============================================================================
// PanelVelocityTangent: dV/d(foil.x[k,dim]) from panel Xj changes ONLY.
// Does NOT include the gamma chain (V_G * dgam). This is needed for the
// wake streamline tangent where the gamma effect on uewiref is already
// captured separately by the gamref chain at fixed wake positions.
// ============================================================================
static Vector2d PanelVelocityTangent(const FoilState& state, const Vector2d& xi,
                                      int k, int dim, const TEInfo& te) {
    int N = state.foil.N;
    Vector2d dV = Vector2d::Zero();

    for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
        int xj_idx = (k == j) ? dim : 2 + dim;
        Matrix2d Xj;
        Xj.col(0) = state.foil.x.col(j);
        Xj.col(1) = state.foil.x.col(j + 1);
        auto vd = PanelLinVortexVelocityDeriv(Xj, xi);
        dV += vd.a_Xj.col(xj_idx) * state.isol.gam(j) +
              vd.b_Xj.col(xj_idx) * state.isol.gam(j + 1);
    }

    if (k == 0 || k == N - 1) {
        int te_idx = (k == N - 1) ? dim : 2 + dim;
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        auto svd = PanelConstSourceVelocityDeriv(XjTE, xi);
        dV += svd.val_Xj.col(te_idx) * 0.5 * te.tcp *
              (state.isol.gam(N - 1) - state.isol.gam(0));
        auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi);
        dV += (vvd.a_Xj.col(te_idx) + vvd.b_Xj.col(te_idx)) * 0.5 * te.tdp *
              (state.isol.gam(0) - state.isol.gam(N - 1));
    }

    double dtcp_dk, dtdp_dk;
    GetTEInfoDeriv(state.foil.x, k, dim, dtcp_dk, dtdp_dk);
    if (std::abs(dtcp_dk) > 1e-20 || std::abs(dtdp_dk) > 1e-20) {
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        Vector2d v_src = PanelConstSourceVelocity(XjTE, xi, nullptr);
        dV += v_src * 0.5 * dtcp_dk * (state.isol.gam(N - 1) - state.isol.gam(0));
        Vector2d a_v, b_v;
        PanelLinVortexVelocity(XjTE, xi, nullptr, false, a_v, b_v);
        dV += (a_v + b_v) * 0.5 * dtdp_dk * (state.isol.gam(0) - state.isol.gam(N - 1));
    }

    return dV;
}

// PanelVelocityVJP: VJP of inviscid velocity w.r.t. foil.x (panel geometry
// + TE factor changes), at fixed gamma and fixed evaluation point xi.
// Given mu_V (adjoint of velocity at xi), accumulates into vjp [2*N].
static void PanelVelocityVJP(const FoilState& state, const Vector2d& xi,
                              const Vector2d& mu_V, const TEInfo& te,
                              VectorXd& vjp) {
    int N = state.foil.N;

    // Regular panels: dV = vd.a_Xj * gam(j) + vd.b_Xj * gam(j+1)
    for (int j = 0; j < N - 1; ++j) {
        Matrix2d Xj;
        Xj.col(0) = state.foil.x.col(j);
        Xj.col(1) = state.foil.x.col(j + 1);
        auto vd = PanelLinVortexVelocityDeriv(Xj, xi);

        double gj = state.isol.gam(j);
        double gjp = state.isol.gam(j + 1);
        for (int dim = 0; dim < 2; ++dim) {
            // Left endpoint (node j)
            Vector2d dV_l = vd.a_Xj.col(dim) * gj + vd.b_Xj.col(dim) * gjp;
            vjp(2 * j + dim) += mu_V.dot(dV_l);
            // Right endpoint (node j+1)
            Vector2d dV_r = vd.a_Xj.col(2 + dim) * gj +
                            vd.b_Xj.col(2 + dim) * gjp;
            vjp(2 * (j + 1) + dim) += mu_V.dot(dV_r);
        }
    }

    // TE source panel
    {
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        double sigma_te =
            0.5 * te.tcp * (state.isol.gam(N - 1) - state.isol.gam(0));

        auto svd = PanelConstSourceVelocityDeriv(XjTE, xi);
        for (int dim = 0; dim < 2; ++dim) {
            // Node N-1 (XjTE col 0)
            vjp(2 * (N - 1) + dim) +=
                mu_V.dot(svd.val_Xj.col(dim) * sigma_te);
            // Node 0 (XjTE col 1)
            vjp(dim) += mu_V.dot(svd.val_Xj.col(2 + dim) * sigma_te);
        }
    }

    // TE vortex panel
    {
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        double gam_te =
            0.5 * te.tdp * (state.isol.gam(0) - state.isol.gam(N - 1));

        auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi);
        for (int dim = 0; dim < 2; ++dim) {
            Vector2d dV = (vvd.a_Xj.col(dim) + vvd.b_Xj.col(dim)) * gam_te;
            vjp(2 * (N - 1) + dim) += mu_V.dot(dV);
            Vector2d dV_r =
                (vvd.a_Xj.col(2 + dim) + vvd.b_Xj.col(2 + dim)) * gam_te;
            vjp(dim) += mu_V.dot(dV_r);
        }
    }

    // TE factor derivatives: tcp and tdp depend on nodes {0, 1, N-2, N-1}
    {
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);

        // Precompute source and vortex velocities at xi from TE panel
        Vector2d v_src = PanelConstSourceVelocity(XjTE, xi, nullptr);
        Vector2d a_v, b_v;
        PanelLinVortexVelocity(XjTE, xi, nullptr, false, a_v, b_v);

        int te_nodes[] = {0, 1, N - 2, N - 1};
        for (int kk : te_nodes) {
            if (kk < 0 || kk >= N) continue;
            for (int dim = 0; dim < 2; ++dim) {
                double dtcp_dk, dtdp_dk;
                GetTEInfoDeriv(state.foil.x, kk, dim, dtcp_dk, dtdp_dk);
                if (std::abs(dtcp_dk) < 1e-20 && std::abs(dtdp_dk) < 1e-20)
                    continue;

                Vector2d dV = Vector2d::Zero();
                if (std::abs(dtcp_dk) > 1e-20) {
                    dV += v_src * 0.5 * dtcp_dk *
                          (state.isol.gam(N - 1) - state.isol.gam(0));
                }
                if (std::abs(dtdp_dk) > 1e-20) {
                    dV += (a_v + b_v) * 0.5 * dtdp_dk *
                          (state.isol.gam(0) - state.isol.gam(N - 1));
                }
                vjp(2 * kk + dim) += mu_V.dot(dV);
            }
        }
    }
}

// Velocity gradient dV/dxi at a point (from all panels + TE)
static Matrix2d VelocityGradient(const FoilState& state, const Vector2d& xi,
                                  const TEInfo& te) {
    int N = state.foil.N;
    Matrix2d G = Matrix2d::Zero();
    for (int j = 0; j < N - 1; ++j) {
        Matrix2d Xj;
        Xj.col(0) = state.foil.x.col(j);
        Xj.col(1) = state.foil.x.col(j + 1);
        auto vd = PanelLinVortexVelocityDeriv(Xj, xi);
        G += vd.a_xi * state.isol.gam(j) + vd.b_xi * state.isol.gam(j + 1);
    }
    {
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        auto svd = PanelConstSourceVelocityDeriv(XjTE, xi);
        G += svd.val_xi * 0.5 * te.tcp * (state.isol.gam(N-1) - state.isol.gam(0));
        auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi);
        G += (vvd.a_xi + vvd.b_xi) * 0.5 * te.tdp * (state.isol.gam(0) - state.isol.gam(N-1));
    }
    return G;
}

// Tangent of normalize(v): d(v/|v|) given dv
static Vector2d NormalizeTangent(const Vector2d& v, const Vector2d& dv) {
    double vn = v.norm();
    if (vn < 1e-10) return Vector2d::Zero();
    Vector2d n = v / vn;
    return (dv - n * n.dot(dv)) / vn;
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

// ============================================================================
// CalcCdLin: compute full dcd/dV [4*Nsys] and dcd/dx_partial [2 x N]
// ============================================================================
void CalcCdLin(FoilState& state) {
    int Nsys = state.glob.Nsys;
    int N = state.foil.N;
    double Vinf = state.param.Vinf > 0 ? state.param.Vinf : 1.0;

    state.post.cd_V = VectorXd::Zero(4 * Nsys);
    state.post.dcd_dx_partial = MatrixXd::Zero(2, N);

    if (!state.oper.viscous || state.glob.U.cols() == 0 ||
        state.vsol.Is.empty()) {
        return;
    }

    // --- Squire-Young cd: depends on last wake state ---
    if (state.vsol.Is.size() > 2 && !state.vsol.Is[2].empty()) {
        int iw = state.vsol.Is[2].back();
        if (iw < state.glob.U.cols()) {
            Vector4d Uw = state.glob.U.col(iw);
            auto [H, H_U] = GetH(Uw);
            auto [uk, uk_U] = GetUk(Uw(3), state.param);
            double uk_Vinf = uk / Vinf;
            double exponent = (5.0 + H) / 2.0;
            double base = std::pow(uk_Vinf, exponent);

            // cd = 2 * theta * (uk/Vinf)^((5+H)/2)
            // dcd/dtheta = 2 * base
            // dcd/dH = 2 * theta * base * 0.5 * log(uk/Vinf)
            // dcd/duk = 2 * theta * exponent * (uk/Vinf)^(exponent-1) / Vinf
            double log_uk_Vinf = SafeLog(uk_Vinf);

            // dcd/dU[iw] via chain rule
            Vector4d dcd_dUw = Vector4d::Zero();
            dcd_dUw(0) = 2.0 * base;  // dcd/dtheta
            dcd_dUw += 2.0 * Uw(0) * base * 0.5 * log_uk_Vinf * H_U;  // through H
            // through uk: dcd/duk * duk/due
            double dcd_duk = 2.0 * Uw(0) * exponent *
                             std::pow(uk_Vinf, exponent - 1.0) / Vinf;
            dcd_dUw(3) += dcd_duk * uk_U(3);  // uk_U is [4] but only ue component

            state.post.cd_V.segment<4>(4 * iw) = dcd_dUw;
        }
    }

    // Note: cd = cd_SY (Squire-Young) is the TOTAL drag from the wake.
    // It implicitly includes skin friction. The cdf/cdp breakdown is computed
    // separately in CalcForce but for the adjoint we only need dcd_SY/dV.
    // No additional skin friction linearization is needed.
}

// ============================================================================
// ComputeDxiDxFoil: sensitivity of arc-length coordinate xi to panel coords
// Returns sparse matrix [Nsys x 2*N]
//
// Chain: foil.x → foil.s → xi (through sstag interpolation)
//
// Arc-length: s(i) = sum_{j=1}^{i} L(j), where L(j) = |x[j] - x[j-1]|
// xi(i) = sgnue(i) * (s(i) - sstag) for foil, wake.s(w) - sstag for wake.
// sstag = w1*s(I0) + w2*s(I1), where I0,I1 are stagnation panel indices.
//
// NOTE: This function is no longer called in the main VJP path (we use the
// efficient ComputeArcLengthVJP instead), but is kept for testing/debugging.
// ============================================================================
SparseMatrixd ComputeDxiDxFoil(const FoilState& state) {
    int Nsys = state.glob.Nsys;
    int N = state.foil.N;
    SparseMatrixd dxi_dx(Nsys, 2 * N);
    // The efficient VJP path (ComputeArcLengthVJP) is used instead.
    return dxi_dx;
}

// ============================================================================
// ComputeArcLengthVJP: efficient VJP through all arc-length chains
//
// Two chains through foil.s:
//   Chain 1 (xi):       foil.x → foil.s → xi → R_BL (through R_x)
//   Chain 2 (sigma_m):  foil.x → foil.s → sigma_m → ue_m → R_ue
//
// Both chains share the same ds/dx_foil Jacobian, so we combine the
// per-node sensitivities q(k) and use a single suffix-sum pass.
//
// Uses suffix sums for O(N) complexity.
// ============================================================================
VectorXd ComputeArcLengthVJP(const FoilState& state,
                              const VectorXd& lambda_bl,
                              const VectorXd& lambda_ue) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    VectorXd vjp = VectorXd::Zero(2 * N);

    // =====================================================================
    // Chain 1: xi chain — lambda_bl^T * R_x * dxi/ds
    // =====================================================================
    VectorXd c = state.glob.R_x.transpose() * lambda_bl;

    int I0 = state.isol.Istag[0];
    int I1 = state.isol.Istag[1];

    // Stagnation interpolation weights: sstag = w1*s(I0) + w2*s(I1)
    // At convergence, StagPointMove uses ue-based weights (not gamma-based):
    //   w1 = ue(I1) / (ue(I0) + ue(I1)),  w2 = ue(I0) / (ue(I0) + ue(I1))
    double ue0 = state.glob.U(3, I0);
    double ue1 = state.glob.U(3, I1);
    double denom = ue0 + ue1;
    double w1 = ue1 / denom;  // weight on s(I0)
    double w2 = ue0 / denom;  // weight on s(I1)

    // dxi(i)/ds(j) = sgnue(i) * [delta_{ij} - w1*delta_{j,I0} - w2*delta_{j,I1}]
    // dxi(N+w)/ds(j) = delta_{j,N-1} - w1*delta_{j,I0} - w2*delta_{j,I1}
    double c_sgnue_sum = 0.0;
    for (int i = 0; i < N; ++i) {
        c_sgnue_sum += c(i) * state.isol.sgnue(i);
    }
    double c_wake_sum = 0.0;
    for (int w = 0; w < Nw; ++w) {
        c_wake_sum += c(N + w);
    }
    double sstag_coeff = c_sgnue_sum + c_wake_sum;

    // q_xi(j) = sensitivity of objective to foil.s(j) through the xi chain
    VectorXd q = VectorXd::Zero(N);
    for (int j = 0; j < N; ++j) {
        q(j) = c(j) * state.isol.sgnue(j);  // foil diagonal
    }
    q(I0) -= w1 * sstag_coeff;
    q(I1) -= w2 * sstag_coeff;
    q(N - 1) += c_wake_sum;  // wake.s(0) = foil.s(N-1)

    // =====================================================================
    // Chain 2: sigma_m chain — lambda_ue^T * dR_ue/ds
    //
    // R_ue(i) = ue(i) - ueinv(i) - sum_j ue_m(i,j) * mass(j)
    // where mass(j) = delta*(j) * ue(j)
    // ue_m = sgnue .* ue_sigma * sigma_m
    //
    // sigma_m(p,p) = -sgnue(p)/ds_p,  sigma_m(p,p+1) = sgnue(p+1)/ds_p
    // for foil panel p (0 <= p < N-1), where ds_p = s(p+1) - s(p).
    //
    // f(k) = lambda_ue^T * dR_ue/ds(k)
    //      = -a(k-1) * beta(k-1)/ds_{k-1}^2  [panel k-1, if k >= 1]
    //      +  a(k)   * beta(k)  /ds_k^2       [panel k,   if k < N-1]
    //
    // where a(p) = (lambda_ue .* sgnue_full)^T * ue_sigma(:,p)
    //       beta(p) = sgnue(p)*mass(p) - sgnue(p+1)*mass(p+1)
    // =====================================================================

    // Build mass vector: mass(j) = delta*(j) * ue(j)
    VectorXd mass(Nsys);
    for (int j = 0; j < Nsys; ++j) {
        mass(j) = state.glob.U(1, j) * state.glob.U(3, j);
    }

    // Build sgnue_full [Nsys]: sgnue for foil, +1 for wake
    VectorXd sgnue_full(Nsys);
    sgnue_full.head(N) = state.isol.sgnue;
    sgnue_full.tail(Nw).setOnes();

    // alpha = lambda_ue .* sgnue_full [Nsys]
    VectorXd alpha = lambda_ue.cwiseProduct(sgnue_full);

    // a = ue_sigma^T * alpha [Npanel]
    // Only need foil panel entries a(0)..a(N-2)
    VectorXd a = state.vsol.ue_sigma.transpose() * alpha;

    // Precompute foil panel spacings ds_p = s(p+1) - s(p) for p=0..N-2
    // and beta(p) = sgnue(p)*mass(p) - sgnue(p+1)*mass(p+1)
    for (int k = 0; k < N; ++k) {
        double fk = 0.0;

        // Panel p = k-1 contribution (if k >= 1 and k-1 < N-1, i.e., k < N)
        if (k >= 1) {
            int p = k - 1;
            double ds_p = state.foil.s(p + 1) - state.foil.s(p);
            if (std::abs(ds_p) > 1e-30) {
                double beta_p = state.isol.sgnue(p) * mass(p) -
                                state.isol.sgnue(p + 1) * mass(p + 1);
                fk += -a(p) * beta_p / (ds_p * ds_p);
            }
        }

        // Panel p = k contribution (if k < N-1)
        if (k < N - 1) {
            int p = k;
            double ds_p = state.foil.s(p + 1) - state.foil.s(p);
            if (std::abs(ds_p) > 1e-30) {
                double beta_p = state.isol.sgnue(p) * mass(p) -
                                state.isol.sgnue(p + 1) * mass(p + 1);
                fk += a(p) * beta_p / (ds_p * ds_p);
            }
        }

        q(k) += fk;
    }

    // Debug: print q(k) and its xi/sigma components at key nodes
    {
        int dbg_nodes[] = {20, 47, 48, 50, 80};
        for (int kk : dbg_nodes) {
            if (kk >= N) continue;
            // Recompute xi-only and sigma-only contributions
            double q_xi = c(kk) * state.isol.sgnue(kk);
            if (kk == I0) q_xi -= w1 * sstag_coeff;
            if (kk == I1) q_xi -= w2 * sstag_coeff;
            if (kk == N - 1) q_xi += c_wake_sum;
            double q_sigma = q(kk) - q_xi;
            std::fprintf(stderr, "  q(%d) = %.6e  (xi=%.6e  sigma=%.6e)\n",
                         kk, q(kk), q_xi, q_sigma);
        }
    }

    // =====================================================================
    // Chain both through ds/dx_foil using suffix sums
    // s(i) = sum_{j=1}^{i} L(j), where L(j) = |x[j] - x[j-1]|
    // sum_i q(i) * ds(i)/dL(j) = Q(j) = sum_{i>=j} q(i)  [suffix sum]
    // vjp(2*k+d) = Q(k)*dL(k)/dx(k,d) + Q(k+1)*dL(k+1)/dx(k,d)
    // =====================================================================

    VectorXd L(N);
    L(0) = 0.0;
    for (int j = 1; j < N; ++j) {
        L(j) = (state.foil.x.col(j) - state.foil.x.col(j - 1)).norm();
    }

    VectorXd Q(N);
    Q(N - 1) = q(N - 1);
    for (int j = N - 2; j >= 1; --j) {
        Q(j) = Q(j + 1) + q(j);
    }
    Q(0) = Q(1) + q(0);

    for (int d = 0; d < 2; ++d) {
        for (int k = 0; k < N; ++k) {
            double val = 0.0;

            if (k >= 1 && L(k) > 1e-30) {
                double dLk = (state.foil.x(d, k) - state.foil.x(d, k - 1)) / L(k);
                val += Q(k) * dLk;
            }

            if (k < N - 1 && L(k + 1) > 1e-30) {
                double dLk1 = -(state.foil.x(d, k + 1) - state.foil.x(d, k)) / L(k + 1);
                val += Q(k + 1) * dLk1;
            }

            vjp(2 * k + d) = val;
        }
    }

    return vjp;
}

// ============================================================================
// ComputeDRueDxFoil: sensitivity of 4th equation R_ue to panel coords
// R_ue(i) = ue(i) - ueinv(i) - sum_j ue_m(i,j) * ds(j) * ue(j)
// Returns [Nsys x 2*N]
// ============================================================================
MatrixXd ComputeDRueDxFoil(FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    double alpha = state.oper.alpha;

    MatrixXd dRue_dx = MatrixXd::Zero(Nsys, 2 * N);

    // --- (1) dueinv/dx_foil ---
    // For airfoil nodes: ueinv(i) = sgnue(i) * gamref(i,:) * [cos(a), sin(a)]
    // gamref depends on AIC which depends on x. Use implicit differentiation:
    // A * gamref = rhs => dgamref/dx = -A^{-1} * (dA/dx*gamref - drhs/dx)
    //
    // For efficiency, compute dueinv/dx_foil as a VJP through the AIC system.
    // For each airfoil node i:
    //   dueinv(i)/dx = sgnue(i) * [cos(a), sin(a)]^T * dgamref(i)/dx
    //
    // We need lambda such that A^T * lambda = e_i (unit vector at node i)
    // for each airfoil node. This is N solves — too expensive.
    //
    // Instead, compute the full dgamref/dx product:
    // dueinv/dx_foil = sgnue .* (cos(a) * dgamref0/dx + sin(a) * dgamref1/dx)
    //
    // dgamref_col/dx = -A^{-1} * d(A*gamref_col - rhs_col)/dx
    //
    // This requires N+1 x 2N matrix d(A*gamref_col)/dx for each ref column.
    // We can compute -A^{-1} * v for any v efficiently since A is already
    // factored. But d(A*gamref_col)/dx is itself an O(N^2) x 2N object.
    //
    // Efficient approach: for each adjoint vector lambda_ue (from the viscous
    // adjoint), the VJP lambda_ue^T * dueinv/dx_foil can be computed in O(N^2)
    // by first solving A^T * mu = diag(sgnue) * [cos(a)*lambda_ue, sin(a)*lambda_ue]
    // and then computing the AIC VJP.
    //
    // However, for the explicit matrix assembly approach, we compute column by
    // column using the AIC inverse which is already available.

    const MatrixXd& AIC = state.isol.AIC;
    auto lu = AIC.partialPivLu();

    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    // For each reference column (0° and 90°), compute d(A*gamref_col - rhs_col)/dx
    // and multiply by -A^{-1} to get dgamref_col/dx
    // This is expensive O(N^2 * 2N) but straightforward.
    // We'll compute the product dueinv/dx = sgnue .* (cosa*dgamref0/dx + sina*dgamref1/dx)
    // by computing the VJP for each design variable.

    // Actually, let's use a column-wise approach for the VJP.
    // For design variable x_k (k = 0..2N-1):
    //   d(A*gamref_col)/dx_k = dA/dx_k * gamref_col  (rhs deriv handled separately)
    //   dgamref_col/dx_k = -A^{-1} * (dA/dx_k * gamref_col - drhs_col/dx_k)
    //
    // This is still O(N^2) per design variable. For 2N variables, total O(N^3).
    // Since A^{-1} can be precomputed, each column is one mat-vec.

    // Precompute A^{-1} (dense, N+1 x N+1)
    MatrixXd Ainv = lu.solve(MatrixXd::Identity(N + 1, N + 1));

    // Precompute TE info
    TEInfo te = GetTEInfo(state.foil.x);
    bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);

    // For each foil node k and each coordinate dimension:
    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int xidx = 2 * k + dim;

            // Compute d(A*gamref[:,c] - rhs[:,c])/d(x_k_dim) for c=0,1
            VectorXd dAgam0 = VectorXd::Zero(N + 1);
            VectorXd dAgam1 = VectorXd::Zero(N + 1);

            // Perturbing node k affects:
            // (a) All rows i: via panels adjacent to k (Xj endpoint change)
            // (b) Row k only: via xi evaluation point change for ALL panels
            // (c) TE panels if k==0 or k==N-1

            // --- (a) Xj endpoint changes: panels k-1 and k ---
            for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
                int xj_idx = (k == j) ? dim : 2 + dim;

                for (int i = 0; i < N; ++i) {
                    if (nogap && i == N - 1) continue;
                    Vector2d xi = state.foil.x.col(i);
                    Matrix2d Xj;
                    Xj.col(0) = state.foil.x.col(j);
                    Xj.col(1) = state.foil.x.col(j + 1);

                    auto deriv = PanelLinVortexStreamDeriv(Xj, xi);

                    for (int c = 0; c < 2; ++c) {
                        double gj = state.isol.gamref(j, c);
                        double gjp = state.isol.gamref(j + 1, c);
                        double val = deriv.a_Xj(xj_idx) * gj +
                                    deriv.b_Xj(xj_idx) * gjp;
                        if (c == 0) dAgam0(i) += val;
                        else dAgam1(i) += val;
                    }
                }
            }

            // --- (b) xi evaluation point change: row k, ALL panels ---
            if (!nogap || k != N - 1) {
                Vector2d xi = state.foil.x.col(k);

                for (int j = 0; j < N - 1; ++j) {
                    Matrix2d Xj;
                    Xj.col(0) = state.foil.x.col(j);
                    Xj.col(1) = state.foil.x.col(j + 1);

                    auto deriv = PanelLinVortexStreamDeriv(Xj, xi);

                    for (int c = 0; c < 2; ++c) {
                        double gj = state.isol.gamref(j, c);
                        double gjp = state.isol.gamref(j + 1, c);
                        double val = deriv.a_xi(dim) * gj +
                                    deriv.b_xi(dim) * gjp;
                        if (c == 0) dAgam0(k) += val;
                        else dAgam1(k) += val;
                    }
                }

                // xi change also affects TE panels at row k
                {
                    Matrix2d XjTE;
                    XjTE.col(0) = state.foil.x.col(N - 1);
                    XjTE.col(1) = state.foil.x.col(0);

                    // TE source
                    auto src_d = PanelConstSourceStreamDeriv(XjTE, xi);
                    for (int c = 0; c < 2; ++c) {
                        double sigma_factor = 0.5 * te.tcp *
                            (state.isol.gamref(N - 1, c) - state.isol.gamref(0, c));
                        if (c == 0) dAgam0(k) += src_d.val_xi(dim) * sigma_factor;
                        else dAgam1(k) += src_d.val_xi(dim) * sigma_factor;
                    }

                    // TE vortex
                    auto vor_d = PanelLinVortexStreamDeriv(XjTE, xi);
                    for (int c = 0; c < 2; ++c) {
                        double gam_factor = 0.5 * te.tdp *
                            (state.isol.gamref(0, c) - state.isol.gamref(N - 1, c));
                        double dab_xi = vor_d.a_xi(dim) + vor_d.b_xi(dim);
                        if (c == 0) dAgam0(k) += dab_xi * gam_factor;
                        else dAgam1(k) += dab_xi * gam_factor;
                    }
                }
            }

            // --- (c) TE panel Xj endpoint changes (k==0 or k==N-1) ---
            if (k == 0 || k == N - 1) {
                Matrix2d XjTE;
                XjTE.col(0) = state.foil.x.col(N - 1);
                XjTE.col(1) = state.foil.x.col(0);

                int te_idx = (k == N - 1) ? dim : 2 + dim;

                for (int i = 0; i < N; ++i) {
                    if (nogap && i == N - 1) continue;
                    Vector2d xi = state.foil.x.col(i);

                    // TE source
                    auto src_d = PanelConstSourceStreamDeriv(XjTE, xi);
                    for (int c = 0; c < 2; ++c) {
                        double sigma_factor = 0.5 * te.tcp *
                            (state.isol.gamref(N - 1, c) - state.isol.gamref(0, c));
                        if (c == 0) dAgam0(i) += src_d.val_Xj(te_idx) * sigma_factor;
                        else dAgam1(i) += src_d.val_Xj(te_idx) * sigma_factor;
                    }

                    // TE vortex
                    auto vor_d = PanelLinVortexStreamDeriv(XjTE, xi);
                    for (int c = 0; c < 2; ++c) {
                        double gam_factor = 0.5 * te.tdp *
                            (state.isol.gamref(0, c) - state.isol.gamref(N - 1, c));
                        double dab = vor_d.a_Xj(te_idx) + vor_d.b_Xj(te_idx);
                        if (c == 0) dAgam0(i) += dab * gam_factor;
                        else dAgam1(i) += dab * gam_factor;
                    }
                }
            }

            // --- (d) TE factor derivatives: d(tcp)/dx and d(tdp)/dx ---
            // tcp and tdp enter the AIC through the TE source/vortex panels.
            // When tcp or tdp change, the AIC matrix changes at ALL rows.
            // Nodes 0, 1, N-2, N-1 can affect tcp/tdp.
            {
                double dtcp_dk, dtdp_dk;
                GetTEInfoDeriv(state.foil.x, k, dim, dtcp_dk, dtdp_dk);

                if (std::abs(dtcp_dk) > 1e-20 || std::abs(dtdp_dk) > 1e-20) {
                    Matrix2d XjTE;
                    XjTE.col(0) = state.foil.x.col(N - 1);
                    XjTE.col(1) = state.foil.x.col(0);

                    for (int i = 0; i < N; ++i) {
                        if (nogap && i == N - 1) continue;
                        Vector2d xi = state.foil.x.col(i);

                        // TE source: A contribution = a_src * 0.5*tcp*(gamref(N-1,c)-gamref(0,c))
                        // d/d(tcp) → a_src * 0.5*(gamref(N-1,c)-gamref(0,c)) * dtcp_dk
                        double a_src = PanelConstSourceStream(XjTE, xi);
                        for (int c = 0; c < 2; ++c) {
                            double dsigma = 0.5 * dtcp_dk *
                                (state.isol.gamref(N - 1, c) - state.isol.gamref(0, c));
                            if (c == 0) dAgam0(i) += a_src * dsigma;
                            else dAgam1(i) += a_src * dsigma;
                        }

                        // TE vortex: A contribution = (a_vor+b_vor) * 0.5*tdp*(gamref(0,c)-gamref(N-1,c))
                        double a_vor, b_vor;
                        PanelLinVortexStream(XjTE, xi, a_vor, b_vor);
                        double ab = a_vor + b_vor;
                        for (int c = 0; c < 2; ++c) {
                            double dgam = 0.5 * dtdp_dk *
                                (state.isol.gamref(0, c) - state.isol.gamref(N - 1, c));
                            if (c == 0) dAgam0(i) += ab * dgam;
                            else dAgam1(i) += ab * dgam;
                        }
                    }
                }
            }

            // --- RHS derivatives ---
            // rhs(k,0) = -xi_z, rhs(k,1) = xi_x (only row k affected)
            if (dim == 0) {
                dAgam1(k) -= 1.0;  // -drhs(k,1)/dx_k_x = -1
            } else {
                dAgam0(k) += 1.0;  // -drhs(k,0)/dx_k_z = +1
            }

            // dgamref_col/dx_k = -A^{-1} * dR_col
            VectorXd dgamref0 = -Ainv * dAgam0;
            VectorXd dgamref1 = -Ainv * dAgam1;




            // dRue(i)/dx_k = -dueinv(i)/dx_k
            // dueinv(i) = sgnue(i) * (cosa*gamref0(i) + sina*gamref1(i))
            for (int i = 0; i < N; ++i) {
                dRue_dx(i, xidx) -= state.isol.sgnue(i) *
                    (cosa * dgamref0(i) + sina * dgamref1(i));
            }

            // Wake ueinv derivatives (missing from original code)
            // Wake station 0: continuity ueinv(N) = sgnue(N-1) * gamref(N-1,:) * cs
            if (state.wake.N > 0) {
                dRue_dx(N, xidx) -= state.isol.sgnue(N - 1) *
                    (cosa * dgamref0(N - 1) + sina * dgamref1(N - 1));

                // Wake stations i > 0: uewiref(i,c) = V_c(wake.x[i]) · t[i]
                // Two contributions: (a) gamref chain, (b) panel Xj chain
                for (int i = 1; i < state.wake.N; ++i) {
                    Vector2d xi_w = state.wake.x.col(i);
                    Vector2d ti_w = state.wake.t.col(i);

                    // (a) Gamref chain: Cgam(i,:) · dgamref(:,c)
                    MatrixXd V_G;
                    InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0,
                                     xi_w, true, &V_G);
                    double duew0 = 0.0, duew1 = 0.0;
                    for (int j = 0; j < N; ++j) {
                        double Cgam_ij = V_G(0, j) * ti_w(0) + V_G(1, j) * ti_w(1);
                        duew0 += Cgam_ij * dgamref0(j);
                        duew1 += Cgam_ij * dgamref1(j);
                    }

                    // (b) Panel Xj chain: panels adjacent to k
                    for (int c = 0; c < 2; ++c) {
                        for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
                            int xj_idx = (k == j) ? dim : 2 + dim;
                            Matrix2d Xj;
                            Xj.col(0) = state.foil.x.col(j);
                            Xj.col(1) = state.foil.x.col(j + 1);
                            auto vd = PanelLinVortexVelocityDeriv(Xj, xi_w);
                            Vector2d dV = vd.a_Xj.col(xj_idx) * state.isol.gamref(j, c) +
                                          vd.b_Xj.col(xj_idx) * state.isol.gamref(j + 1, c);
                            if (c == 0) duew0 += dV.dot(ti_w);
                            else duew1 += dV.dot(ti_w);
                        }
                        // TE panels
                        if (k == 0 || k == N - 1) {
                            int te_idx = (k == N - 1) ? dim : 2 + dim;
                            Matrix2d XjTE;
                            XjTE.col(0) = state.foil.x.col(N - 1);
                            XjTE.col(1) = state.foil.x.col(0);
                            auto svd_te = PanelConstSourceVelocityDeriv(XjTE, xi_w);
                            double sigma_c = 0.5 * te.tcp *
                                (state.isol.gamref(N-1, c) - state.isol.gamref(0, c));
                            auto vvd_te = PanelLinVortexVelocityDeriv(XjTE, xi_w);
                            double gam_c = 0.5 * te.tdp *
                                (state.isol.gamref(0, c) - state.isol.gamref(N-1, c));
                            Vector2d dV_te = svd_te.val_Xj.col(te_idx) * sigma_c +
                                (vvd_te.a_Xj.col(te_idx) + vvd_te.b_Xj.col(te_idx)) * gam_c;
                            if (c == 0) duew0 += dV_te.dot(ti_w);
                            else duew1 += dV_te.dot(ti_w);
                        }
                        // TE factor derivatives
                        {
                            double dtcp_dk, dtdp_dk;
                            GetTEInfoDeriv(state.foil.x, k, dim, dtcp_dk, dtdp_dk);
                            if (std::abs(dtcp_dk) > 1e-20 || std::abs(dtdp_dk) > 1e-20) {
                                Matrix2d XjTE;
                                XjTE.col(0) = state.foil.x.col(N-1);
                                XjTE.col(1) = state.foil.x.col(0);
                                Vector2d v_src = PanelConstSourceVelocity(XjTE, xi_w, nullptr);
                                double dsigma = 0.5 * dtcp_dk *
                                    (state.isol.gamref(N-1, c) - state.isol.gamref(0, c));
                                Vector2d a_v, b_v;
                                PanelLinVortexVelocity(XjTE, xi_w, nullptr, false, a_v, b_v);
                                double dgam = 0.5 * dtdp_dk *
                                    (state.isol.gamref(0, c) - state.isol.gamref(N-1, c));
                                Vector2d dV_te = v_src * dsigma + (a_v + b_v) * dgam;
                                if (c == 0) duew0 += dV_te.dot(ti_w);
                                else duew1 += dV_te.dot(ti_w);
                            }
                        }
                    }

                    dRue_dx(N + i, xidx) -= (cosa * duew0 + sina * duew1);
                }
                // Wake geometry correction (uewiref chain through wake.x/wake.t)
                // is now handled entirely by ComputeWakeGeomVJP backward sweep.
            }
        }
    }

    // --- (2) d(ue_m * mass)/dx_foil ---
    // ue_m = diag(sgnue) * ue_sigma * sigma_m
    // sigma_m is fixed (foil.s doesn't change with foil.x)
    // So d(ue_m * mass)/dx = diag(sgnue) * d(ue_sigma)/dx * phi
    // where phi = sigma_m * mass
    //
    // ue_sigma structure:
    //   foil rows: Bp = -AIC^{-1} * B
    //   wake rows: Dw = Cgam*Bp + Csig, with Dw(0) = Bp(N-1)
    //
    // For d(Bp*phi)/dx_k: from AIC*(-Bp_full) = B, differentiate:
    //   d(Bp*phi)/dx_k = -Ainv * (dB_phi/dx_k + dAIC_Bpphi/dx_k)
    // where dB_phi = d(B*phi)/dx_k and dAIC_Bpphi = dAIC/dx_k * (Bp*phi)_full
    //
    // For d(Dw*phi)/dx_k = d(Cgam)/dx_k * (Bp*phi) + Cgam * d(Bp*phi)/dx_k
    //                     + d(Csig*phi)/dx_k
    //
    // d(Cgam)/dx_k and d(Csig*phi)/dx_k involve velocity influence derivatives
    // at wake evaluation points from foil panels.

    int Nw = state.wake.N;
    if (Nw > 0) {  // Re-enabled
        // Compute phi = sigma_m * mass [N+Nw-2]
        VectorXd ue_vec = state.glob.U.row(3).transpose();
        VectorXd ds_vec = state.glob.U.row(1).transpose();
        VectorXd mass = ds_vec.cwiseProduct(ue_vec);
        VectorXd phi = state.vsol.sigma_m * mass;

        // Bp_phi_full = AIC^{-1} * (B * phi) extended to N+1
        VectorXd Bphi = VectorXd::Zero(N + 1);
        // B is [N+1 x Npanel], but only first N rows are nonzero
        // Bp * phi = -(AIC^{-1} * B * phi)(0:N-1)
        // We need (B * phi) computed from the stored B...
        // Actually, we don't store B. Recompute B*phi directly.

        // B(i,j) contributions: j < N-1 are const source stream from foil panels
        // j >= N-1 are lin source stream from wake panels
        for (int i = 0; i < N; ++i) {
            Vector2d xi_eval = state.foil.x.col(i);
            // Foil source panels
            for (int j = 0; j < N - 1; ++j) {
                Matrix2d Xj;
                Xj.col(0) = state.foil.x.col(j);
                Xj.col(1) = state.foil.x.col(j + 1);
                Bphi(i) += PanelConstSourceStream(Xj, xi_eval) * phi(j);
            }
            // Wake source panels (frozen geometry, won't contribute to dx_foil)
            for (int j = 0; j < Nw - 1; ++j) {
                Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
                Matrix2d Xj_l, Xj_r;
                Xj_l.col(0) = state.wake.x.col(j);
                Xj_l.col(1) = Xm;
                Vector2d X_right;
                if (j == Nw - 2) {
                    X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
                } else {
                    X_right = state.wake.x.col(j + 1);
                }
                Xj_r.col(0) = Xm;
                Xj_r.col(1) = X_right;

                double a_l, b_l, a_r, b_r;
                PanelLinSourceStream(Xj_l, xi_eval, a_l, b_l);
                PanelLinSourceStream(Xj_r, xi_eval, a_r, b_r);

                // Accumulate into Bphi using same assembly as CalcUeM step 2
                if (j > 0) {
                    Bphi(i) += (0.5 * a_l + b_l) * phi(N - 1 + j);
                    Bphi(i) += 0.5 * a_l * phi(N - 1 + j - 1);
                } else {
                    Bphi(i) += b_l * phi(N - 1 + j);
                }
                Bphi(i) += (a_r + 0.5 * b_r) * phi(N - 1 + j);
                if (j < Nw - 2) {
                    Bphi(i) += 0.5 * b_r * phi(N - 1 + j + 1);
                } else {
                    Bphi(i) += 0.5 * b_r * phi(N - 1 + j);
                }
            }
        }

        // Bp_phi_full = -Ainv * Bphi (the full N+1 vector)
        VectorXd Bp_phi_full = -Ainv * Bphi;
        VectorXd Bp_phi = Bp_phi_full.head(N);  // [N]

        // Precompute Cgam * Bp_phi for wake contribution
        // Cgam(i,:) * Bp_phi is the wake velocity from foil vortex sheets
        // Cgam was computed in CalcUeM but not stored separately.
        // Recompute Cgam * Bp_phi:
        VectorXd Cgam_Bpphi = VectorXd::Zero(Nw);
        for (int i = 0; i < Nw; ++i) {
            Vector2d xi_w = state.wake.x.col(i);
            Vector2d ti_w = state.wake.t.col(i);
            MatrixXd V_G;
            InviscidVelocity(state.foil.x, Bp_phi, 0.0, 0.0, xi_w, true, &V_G);
            // Cgam_Bpphi(i) = sum_j V_G(:,j) . ti * Bp_phi(j) — but InviscidVelocity
            // already multiplied by Bp_phi as the "gamma" vector.
            // Actually InviscidVelocity computes V = sum_j (a_j*G_j + b_j*G_{j+1})
            // We want Cgam(i,:) * Bp_phi = tangent component of velocity from Bp_phi vortices
            Vector2d V = InviscidVelocity(state.foil.x, Bp_phi, 0.0, 0.0,
                                          xi_w, false, nullptr);
            Cgam_Bpphi(i) = V.dot(ti_w);
        }

        // Csig * phi: wake velocity from source panels
        VectorXd Csig_phi = VectorXd::Zero(Nw);
        // This involves foil source panels and wake source panels.
        // Only foil panel contributions depend on foil.x.
        // We need Csig_phi for the Dw contribution.
        // Csig was computed in CalcUeM. Recompute Csig*phi:
        for (int i = 0; i < Nw; ++i) {
            Vector2d xi_w = state.wake.x.col(i);
            Vector2d ti_w = state.wake.t.col(i);
            // Foil source panels
            int jstart = (i == 0) ? 1 : 0;
            int jend = (i == 0) ? N - 2 : N - 1;
            for (int j = jstart; j < jend; ++j) {
                Matrix2d Xj;
                Xj.col(0) = state.foil.x.col(j);
                Xj.col(1) = state.foil.x.col(j + 1);
                Vector2d v = PanelConstSourceVelocity(Xj, xi_w, &ti_w);
                Csig_phi(i) += v(0) * phi(j);
            }
            // Wake source panels (frozen, skip for dx computation)
            // ... (include for correct Csig_phi but not needed for derivatives)
        }

        // sgnue vector
        VectorXd sgue(Nsys);
        sgue.head(N) = state.isol.sgnue;
        sgue.tail(Nw).setOnes();

        // Now for each design variable x_k, compute d(ue_m * mass)/dx_k
        for (int k = 0; k < N; ++k) {
            for (int dim = 0; dim < 2; ++dim) {
                int xidx = 2 * k + dim;

                // --- d(B*phi)/dx_k: source stream derivatives ---
                // B(i,j) = PanelConstSourceStream(Xj_j, xi_i) for j<N-1
                // Perturbing node k affects panels j=k-1,k (Xj endpoints)
                // and row i=k (xi evaluation point)
                VectorXd dBphi = VectorXd::Zero(N + 1);

                // Xj endpoint contributions (panels adjacent to k)
                for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
                    int xj_idx = (k == j) ? dim : 2 + dim;
                    for (int i = 0; i < N; ++i) {
                        if (nogap && i == N - 1) continue;  // no-gap: row N-1 is extrapolation, not streamfn
                        Vector2d xi_eval = state.foil.x.col(i);
                        Matrix2d Xj;
                        Xj.col(0) = state.foil.x.col(j);
                        Xj.col(1) = state.foil.x.col(j + 1);
                        auto sd = PanelConstSourceStreamDeriv(Xj, xi_eval);
                        dBphi(i) += sd.val_Xj(xj_idx) * phi(j);
                    }
                }

                // xi evaluation point contributions (row k, ALL panels in B)
                // Skip for no-gap row N-1 (extrapolation, not streamfunction)
                if (!nogap || k != N - 1)
                {
                    Vector2d xi_eval = state.foil.x.col(k);

                    // Foil const source panels
                    for (int j = 0; j < N - 1; ++j) {
                        Matrix2d Xj_src;
                        Xj_src.col(0) = state.foil.x.col(j);
                        Xj_src.col(1) = state.foil.x.col(j + 1);
                        auto sd = PanelConstSourceStreamDeriv(Xj_src, xi_eval);
                        dBphi(k) += sd.val_xi(dim) * phi(j);
                    }

                    // Wake linear source panels: Xj are wake geometry (frozen),
                    // but xi = foil.x.col(k) changes. Need LinSourceStream xi derivs.
                    for (int j = 0; j < Nw - 1; ++j) {
                        Vector2d Xm_w = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
                        Matrix2d Xj_wl, Xj_wr;
                        Xj_wl.col(0) = state.wake.x.col(j);
                        Xj_wl.col(1) = Xm_w;
                        Vector2d X_right_w;
                        if (j == Nw - 2) {
                            X_right_w = 2.0 * state.wake.x.col(j + 1) - Xm_w;
                        } else {
                            X_right_w = state.wake.x.col(j + 1);
                        }
                        Xj_wr.col(0) = Xm_w;
                        Xj_wr.col(1) = X_right_w;

                        auto ld = PanelLinSourceStreamDeriv(Xj_wl, xi_eval);
                        auto rd = PanelLinSourceStreamDeriv(Xj_wr, xi_eval);

                        // Assembly mirrors CalcUeM step 2 exactly:
                        // Left half-panel [wake(j), Xm]:
                        //   j>0: B(i, N-1+j) += 0.5*a + b; B(i, N-1+j-1) += 0.5*a
                        //   j==0: B(i, N-1+j) += b
                        if (j > 0) {
                            dBphi(k) += (0.5 * ld.a_xi(dim) + ld.b_xi(dim)) * phi(N - 1 + j);
                            dBphi(k) += 0.5 * ld.a_xi(dim) * phi(N - 1 + j - 1);
                        } else {
                            dBphi(k) += ld.b_xi(dim) * phi(N - 1 + j);
                        }
                        // Right half-panel [Xm, X_right]:
                        //   B(i, N-1+j) += a + 0.5*b
                        //   j<Nw-2: B(i, N-1+j+1) += 0.5*b
                        //   j==Nw-2: B(i, N-1+j) += 0.5*b  (ghost fold-back)
                        dBphi(k) += (rd.a_xi(dim) + 0.5 * rd.b_xi(dim)) * phi(N - 1 + j);
                        if (j < Nw - 2) {
                            dBphi(k) += 0.5 * rd.b_xi(dim) * phi(N - 1 + j + 1);
                        } else {
                            dBphi(k) += 0.5 * rd.b_xi(dim) * phi(N - 1 + j);
                        }
                    }
                }

                // --- d(AIC)/dx_k * Bp_phi_full ---
                // Same structure as dAIC/dx_k * gamref, but with Bp_phi_full
                VectorXd dAIC_Bpphi = VectorXd::Zero(N + 1);

                // Xj endpoint contributions
                for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
                    int xj_idx = (k == j) ? dim : 2 + dim;
                    for (int i = 0; i < N; ++i) {
                        if (nogap && i == N - 1) continue;
                        Vector2d xi_eval = state.foil.x.col(i);
                        Matrix2d Xj;
                        Xj.col(0) = state.foil.x.col(j);
                        Xj.col(1) = state.foil.x.col(j + 1);
                        auto vd = PanelLinVortexStreamDeriv(Xj, xi_eval);
                        dAIC_Bpphi(i) += vd.a_Xj(xj_idx) * Bp_phi_full(j) +
                                         vd.b_Xj(xj_idx) * Bp_phi_full(j + 1);
                    }
                }

                // xi evaluation point (row k)
                if (!nogap || k != N - 1) {
                    Vector2d xi_eval = state.foil.x.col(k);
                    for (int j = 0; j < N - 1; ++j) {
                        Matrix2d Xj;
                        Xj.col(0) = state.foil.x.col(j);
                        Xj.col(1) = state.foil.x.col(j + 1);
                        auto vd = PanelLinVortexStreamDeriv(Xj, xi_eval);
                        dAIC_Bpphi(k) += vd.a_xi(dim) * Bp_phi_full(j) +
                                         vd.b_xi(dim) * Bp_phi_full(j + 1);
                    }
                    // TE panels at row k
                    {
                        Matrix2d XjTE;
                        XjTE.col(0) = state.foil.x.col(N - 1);
                        XjTE.col(1) = state.foil.x.col(0);
                        double src_coeff = 0.5 * te.tcp * (Bp_phi_full(N - 1) - Bp_phi_full(0));
                        auto sd = PanelConstSourceStreamDeriv(XjTE, xi_eval);
                        dAIC_Bpphi(k) += sd.val_xi(dim) * src_coeff;
                        double vor_coeff = 0.5 * te.tdp * (Bp_phi_full(0) - Bp_phi_full(N - 1));
                        auto vd = PanelLinVortexStreamDeriv(XjTE, xi_eval);
                        dAIC_Bpphi(k) += (vd.a_xi(dim) + vd.b_xi(dim)) * vor_coeff;
                    }
                }

                // TE panel Xj endpoint changes
                if (k == 0 || k == N - 1) {
                    int te_idx = (k == N - 1) ? dim : 2 + dim;
                    Matrix2d XjTE;
                    XjTE.col(0) = state.foil.x.col(N - 1);
                    XjTE.col(1) = state.foil.x.col(0);
                    double src_coeff = 0.5 * te.tcp * (Bp_phi_full(N - 1) - Bp_phi_full(0));
                    double vor_coeff = 0.5 * te.tdp * (Bp_phi_full(0) - Bp_phi_full(N - 1));
                    for (int i = 0; i < N; ++i) {
                        if (nogap && i == N - 1) continue;
                        Vector2d xi_eval = state.foil.x.col(i);
                        auto sd = PanelConstSourceStreamDeriv(XjTE, xi_eval);
                        dAIC_Bpphi(i) += sd.val_Xj(te_idx) * src_coeff;
                        auto vd = PanelLinVortexStreamDeriv(XjTE, xi_eval);
                        dAIC_Bpphi(i) += (vd.a_Xj(te_idx) + vd.b_Xj(te_idx)) * vor_coeff;
                    }
                }

                // --- TE factor derivatives for dAIC_Bpphi ---
                {
                    double dtcp_dk, dtdp_dk;
                    GetTEInfoDeriv(state.foil.x, k, dim, dtcp_dk, dtdp_dk);
                    if (std::abs(dtcp_dk) > 1e-20 || std::abs(dtdp_dk) > 1e-20) {
                        Matrix2d XjTE;
                        XjTE.col(0) = state.foil.x.col(N - 1);
                        XjTE.col(1) = state.foil.x.col(0);
                        for (int i = 0; i < N; ++i) {
                            if (nogap && i == N - 1) continue;
                            Vector2d xi_eval = state.foil.x.col(i);
                            double a_src = PanelConstSourceStream(XjTE, xi_eval);
                            dAIC_Bpphi(i) += a_src * 0.5 * dtcp_dk *
                                (Bp_phi_full(N - 1) - Bp_phi_full(0));
                            double a_v, b_v;
                            PanelLinVortexStream(XjTE, xi_eval, a_v, b_v);
                            dAIC_Bpphi(i) += (a_v + b_v) * 0.5 * dtdp_dk *
                                (Bp_phi_full(0) - Bp_phi_full(N - 1));
                        }
                    }
                }

                // d(Bp*phi)/dx_k = -Ainv * (dBphi + dAIC_Bpphi)
                VectorXd dBpphi_full = -Ainv * (dBphi + dAIC_Bpphi);
                VectorXd dBpphi = dBpphi_full.head(N);

                // --- d(Cgam * Bp_phi)/dx_k ---
                // Cgam(i,j) = tangent velocity at wake point i from vortex j
                // Perturbing foil node k affects the vortex velocity influence
                // at wake points through panels adjacent to k and TE panels.
                VectorXd dCgam_Bpphi = VectorXd::Zero(Nw);
                for (int i = 0; i < Nw; ++i) {
                    Vector2d xi_w = state.wake.x.col(i);
                    Vector2d ti_w = state.wake.t.col(i);

                    // Panels adjacent to k
                    for (int j = std::max(0, k - 1); j <= std::min(N - 2, k); ++j) {
                        int xj_idx = (k == j) ? dim : 2 + dim;
                        Matrix2d Xj;
                        Xj.col(0) = state.foil.x.col(j);
                        Xj.col(1) = state.foil.x.col(j + 1);
                        auto vd = PanelLinVortexVelocityDeriv(Xj, xi_w);
                        // dCgam(i,j)/dx_k = d(V_j . ti)/dx_k
                        Vector2d da = vd.a_Xj.col(xj_idx);
                        Vector2d db = vd.b_Xj.col(xj_idx);
                        dCgam_Bpphi(i) += da.dot(ti_w) * Bp_phi(j) +
                                          db.dot(ti_w) * Bp_phi(j + 1);
                    }

                    // TE panels
                    if (k == 0 || k == N - 1) {
                        int te_idx = (k == N - 1) ? dim : 2 + dim;
                        Matrix2d XjTE;
                        XjTE.col(0) = state.foil.x.col(N - 1);
                        XjTE.col(1) = state.foil.x.col(0);

                        // TE source velocity
                        auto svd = PanelConstSourceVelocityDeriv(XjTE, xi_w);
                        double sigma_te = 0.5 * te.tcp * (Bp_phi(N - 1) - Bp_phi(0));
                        dCgam_Bpphi(i) += svd.val_Xj.col(te_idx).dot(ti_w) * sigma_te;
                        // Actually Cgam uses InviscidVelocity which includes TE.
                        // The TE contribution to Cgam is through the TE source and vortex panels.
                        // But Cgam multiplied by Bp_phi is already the velocity from Bp_phi "gammas".
                        // The InviscidVelocity function includes TE panels internally.
                        // So d(Cgam*Bp_phi)/dx_k needs d(InviscidVelocity(foil.x, Bp_phi, ...))/dx_k
                        // which includes TE panel derivatives.

                        // TE vortex velocity
                        auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi_w);
                        Vector2d fvor = vvd.a_Xj.col(te_idx) + vvd.b_Xj.col(te_idx);
                        double gam_te = 0.5 * te.tdp * (Bp_phi(0) - Bp_phi(N - 1));
                        dCgam_Bpphi(i) += fvor.dot(ti_w) * gam_te;
                    }

                    // TE factor derivatives for Cgam
                    {
                        double dtcp_dk, dtdp_dk;
                        GetTEInfoDeriv(state.foil.x, k, dim, dtcp_dk, dtdp_dk);
                        if (std::abs(dtcp_dk) > 1e-20 || std::abs(dtdp_dk) > 1e-20) {
                            Matrix2d XjTE;
                            XjTE.col(0) = state.foil.x.col(N - 1);
                            XjTE.col(1) = state.foil.x.col(0);
                            Vector2d v_src = PanelConstSourceVelocity(XjTE, xi_w, nullptr);
                            dCgam_Bpphi(i) += v_src.dot(ti_w) * 0.5 * dtcp_dk *
                                (Bp_phi(N - 1) - Bp_phi(0));
                            Vector2d a_vor, b_vor;
                            PanelLinVortexVelocity(XjTE, xi_w, nullptr, false, a_vor, b_vor);
                            Vector2d fv = a_vor + b_vor;
                            dCgam_Bpphi(i) += fv.dot(ti_w) * 0.5 * dtdp_dk *
                                (Bp_phi(0) - Bp_phi(N - 1));
                        }
                    }
                }

                // --- d(Csig*phi)/dx_k ---
                // Csig has two parts:
                // (A) Const source panels on airfoil: Csig(i,j) for j<N-1
                //     These depend on foil.x through panel endpoints.
                // (B) Self-influence at wake j=0 (TE point):
                //     Csig(0,0) += (0.5/pi)*(log(dl/d2)+1)
                //     Csig(0,N-2) += (0.5/pi)*(log(du/d2)+1)
                //     where dl = ||foil.x(1)-foil.x(0)||, du = ||foil.x(N-1)-foil.x(N-2)||
                //     These depend on foil.x through dl and du.
                // (C) Non-self wake j=0 (i!=0): uses PanelLinSourceVelocity
                //     at wake Xm/Xr points — frozen, no foil.x dependence.
                // (D) All other wake j>0: use wake geometry only — frozen.

                VectorXd dCsig_phi = VectorXd::Zero(Nw);

                // (A) Foil const source panels
                for (int i = 0; i < Nw; ++i) {
                    Vector2d xi_w = state.wake.x.col(i);
                    Vector2d ti_w = state.wake.t.col(i);

                    int jstart = (i == 0) ? 1 : 0;
                    int jend = (i == 0) ? N - 2 : N - 1;

                    for (int j = std::max(jstart, k - 1);
                         j <= std::min(jend - 1, k); ++j) {
                        int xj_idx = (k == j) ? dim : 2 + dim;
                        Matrix2d Xj_src;
                        Xj_src.col(0) = state.foil.x.col(j);
                        Xj_src.col(1) = state.foil.x.col(j + 1);
                        auto svd = PanelConstSourceVelocityDeriv(Xj_src, xi_w);
                        dCsig_phi(i) += svd.val_Xj.col(xj_idx).dot(ti_w) * phi(j);
                    }
                }

                // (B) Self-influence at wake j=0 (i==0, the TE wake point)
                // Csig(0, 0) += (0.5/pi) * (log(dl/d2) + 1)
                // Csig(0, N-2) += (0.5/pi) * (log(du/d2) + 1)
                // dl = ||foil.x(1) - foil.x(0)||
                // du = ||foil.x(N-1) - foil.x(N-2)||
                // d2 = distance from Xr to Xm at wake j=0 — frozen wake geometry
                {
                    // dl depends on nodes 0 and 1
                    Vector2d dx_l = state.foil.x.col(1) - state.foil.x.col(0);
                    double dl = dx_l.norm();
                    // du depends on nodes N-1 and N-2
                    Vector2d dx_u = state.foil.x.col(N - 1) - state.foil.x.col(N - 2);
                    double du = dx_u.norm();

                    // d(dl)/d(x_k_dim)
                    double ddl = 0.0;
                    if (dl > 1e-12) {
                        if (k == 0) ddl = -dx_l(dim) / dl;
                        else if (k == 1) ddl = dx_l(dim) / dl;
                    }
                    // d(du)/d(x_k_dim)
                    double ddu = 0.0;
                    if (du > 1e-12) {
                        if (k == N - 2) ddu = -dx_u(dim) / du;
                        else if (k == N - 1) ddu = dx_u(dim) / du;
                    }

                    // d(Csig(0,0))/dx_k = (0.5/pi) * (1/dl) * ddl
                    // because d(log(dl/d2))/dx = (1/dl)*ddl (d2 is frozen)
                    if (dl > 1e-12 && ddl != 0.0) {
                        dCsig_phi(0) += (0.5 / kPi) * (1.0 / dl) * ddl * phi(0);
                    }
                    // d(Csig(0,N-2))/dx_k = (0.5/pi) * (1/du) * ddu
                    if (du > 1e-12 && ddu != 0.0) {
                        dCsig_phi(0) += (0.5 / kPi) * (1.0 / du) * ddu * phi(N - 2);
                    }
                }

                // --- Assemble d(ue_sigma * phi)/dx_k ---
                // Foil part: dBpphi
                // Wake part: Cgam * dBpphi + dCgam_Bpphi + dCsig_phi
                // with Dw(0) = Bp(N-1) override
                VectorXd due_sigma_phi = VectorXd::Zero(Nsys);
                due_sigma_phi.head(N) = dBpphi;

                VectorXd dDw = VectorXd::Zero(Nw);
                // Compute Cgam * dBpphi
                for (int i = 0; i < Nw; ++i) {
                    Vector2d xi_w = state.wake.x.col(i);
                    Vector2d ti_w = state.wake.t.col(i);
                    Vector2d V = InviscidVelocity(state.foil.x, dBpphi, 0.0, 0.0,
                                                  xi_w, false, nullptr);
                    dDw(i) = V.dot(ti_w);
                }
                dDw += dCgam_Bpphi + dCsig_phi;
                // Override: Dw(0) = Bp(N-1)
                dDw(0) = dBpphi(N - 1);

                due_sigma_phi.tail(Nw) = dDw;

                // dR_ue/dx has -d(ue_m*mass)/dx = -diag(sgnue)*d(ue_sigma*phi)/dx
                for (int i = 0; i < Nsys; ++i) {
                    dRue_dx(i, xidx) -= sgue(i) * due_sigma_phi(i);
                }
            }
        }
    }

    return dRue_dx;
}

// ============================================================================
// ComputeWakeGeomVJP: VJP through wake geometry dependence on foil.x
//
// The 4th equation R_ue at wake stations depends on wake.x and wake.t through:
//   uewiref[i,c] = V_c(wake.x[i]) · wake.t[i]
//   Cgam(i,j) and Csig(i,j): velocity influences at wake.x[i] projected onto wake.t[i]
//
// ComputeDRueDxFoil treats wake.x and wake.t as constants.  The missing terms are:
//   ∂R_ue/∂x_foil|_full = ∂R_ue/∂x_foil|_{wake fixed}     [already computed]
//                        + ∂R_ue/∂wake.x · dwake.x/dx_foil  [this function]
//                        + ∂R_ue/∂wake.t · dwake.t/dx_foil  [this function]
//
// For the adjoint, we compute the VJP: λ_ue^T · (missing terms).
// ============================================================================
static VectorXd ComputeWakeGeomVJPImpl(FoilState& state,
                                        const VectorXd& lambda_ue,
                                        MatrixXd* mu_x_out,
                                        MatrixXd* mu_t_out) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    double Vinf = state.oper.Vinf;
    double cosa = Cosd(state.oper.alpha);
    double sina = Sind(state.oper.alpha);

    VectorXd vjp = VectorXd::Zero(2 * N);

    if (Nw <= 0) return vjp;

    // Extract wake portion of lambda_ue
    VectorXd lam_w = lambda_ue.tail(Nw);

    // =========================================================================
    // Step 1: Compute wake adjoints μ_x[i] and μ_t[i]
    //
    // For each wake point i, compute how R_ue (weighted by lambda_ue) depends
    // on wake.x[i] and wake.t[i].
    //
    // R_ue at wake point i has two main dependencies on wake geometry:
    // (A) ueinv contribution: uewiref[i,c] = V_c(wake.x[i]) · wake.t[i]
    // (B) ue_m contribution: through Cgam and Csig at wake points
    // =========================================================================

    // mu_x[i] = adjoint w.r.t. wake.x[i] (2-vector per wake point)
    // mu_t[i] = adjoint w.r.t. wake.t[i] (2-vector per wake point)
    MatrixXd mu_x = MatrixXd::Zero(2, Nw);
    MatrixXd mu_t = MatrixXd::Zero(2, Nw);

    // --- (A) ueinv contribution ---
    // R_ue(N+i) has term -ueinv(N+i) = -(uewiref(i,0)*cos(a) + uewiref(i,1)*sin(a))
    // uewiref(i,c) = V_c(wake.x[i]) · wake.t[i]
    //
    // ∂R_ue(N+i)/∂wake.x[i] = -(dV_c/dxi(wake.x[i]) · wake.t[i]) * [cos(a), sin(a)]
    // ∂R_ue(N+i)/∂wake.t[i] = -V_c(wake.x[i]) * [cos(a), sin(a)]
    //
    // Skip i=0: first wake point uses continuity (uewiref overridden by last airfoil)

    for (int i = 0; i < Nw; ++i) {
        double lam_i = lam_w(i);
        if (std::abs(lam_i) < 1e-20) continue;

        // The i=0 wake point: ueinv(N+0) = ueinv(N-1) (continuity override)
        // So ∂R_ue(N+0)/∂wake geometry through ueinv = 0
        // But Cgam/Csig contributions to ue_m at i=0 are also overridden (Dw(0)=Bp(N-1))
        // So i=0 has NO wake geometry dependence at all
        if (i == 0) continue;

        Vector2d xi_w = state.wake.x.col(i);
        Vector2d ti_w = state.wake.t.col(i);

        // Compute velocity and velocity gradient at wake point for each reference
        for (int c = 0; c < 2; ++c) {
            double wt = (c == 0) ? cosa : sina;  // weight for this reference

            // V_c(xi_w) is the velocity from gamref[:,c] + freestream at ref alpha
            Vector2d Vc = InviscidVelocity(state.foil.x, state.isol.gamref.col(c),
                                           Vinf, (c == 0) ? 0.0 : 90.0, xi_w);

            // ∂uewiref(i,c)/∂wake.t[i] = V_c
            // So ∂R_ue(N+i)/∂wake.t[i] contribution = -wt * V_c * lam_i
            mu_t.col(i) -= lam_i * wt * Vc;

            // ∂uewiref(i,c)/∂wake.x[i] = (dV_c/dxi) · t_i
            // Need dV/dxi: velocity gradient at evaluation point
            // Use panel influence derivatives with _xi terms
            TEInfo te = GetTEInfo(state.foil.x);

            // Accumulate dV/dxi for all panels
            Matrix2d dVdxi = Matrix2d::Zero();  // dV/dxi [2x2]

            for (int j = 0; j < N - 1; ++j) {
                Matrix2d Xj;
                Xj.col(0) = state.foil.x.col(j);
                Xj.col(1) = state.foil.x.col(j + 1);

                auto vd = PanelLinVortexVelocityDeriv(Xj, xi_w);
                double gj = state.isol.gamref(j, c);
                double gjp = state.isol.gamref(j + 1, c);

                // dV/dxi from this panel: a_xi * gj + b_xi * gjp
                // a_xi, b_xi are [2x2] matrices
                dVdxi += vd.a_xi * gj + vd.b_xi * gjp;
            }

            // TE source panel
            {
                Matrix2d XjTE;
                XjTE.col(0) = state.foil.x.col(N - 1);
                XjTE.col(1) = state.foil.x.col(0);

                auto svd = PanelConstSourceVelocityDeriv(XjTE, xi_w);
                double sigma_te = 0.5 * te.tcp * (state.isol.gamref(N - 1, c) -
                                                   state.isol.gamref(0, c));
                dVdxi += svd.val_xi * sigma_te;
            }

            // TE vortex panel
            {
                Matrix2d XjTE;
                XjTE.col(0) = state.foil.x.col(N - 1);
                XjTE.col(1) = state.foil.x.col(0);

                auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi_w);
                double gam_factor = 0.5 * te.tdp *
                    (state.isol.gamref(0, c) - state.isol.gamref(N - 1, c));
                Matrix2d fvor_xi = vvd.a_xi + vvd.b_xi;
                dVdxi += fvor_xi * gam_factor;
            }

            // ∂uewiref(i,c)/∂wake.x[i] = dVdxi^T · t_i (chain rule)
            // Actually: uewiref = V · t, so d(uewiref)/d(xi_k) = dV_k/dxi_k · t
            // dVdxi is [2x2]: dVdxi(component, xi_dim)
            // d(V·t)/d(xi_dim) = sum_k dV_k/d(xi_dim) * t_k = dVdxi^T * t
            Vector2d duewiref_dxi = dVdxi.transpose() * ti_w;

            // R_ue(N+i) has -uewiref * wt, so:
            mu_x.col(i) -= lam_i * wt * duewiref_dxi;
        }
    }

    // --- (B) ue_m contribution ---
    // ue_m at wake rows depends on Cgam and Csig, which are evaluated at wake points.
    // The full ue_m VJP through wake geometry is more complex: it involves
    // differentiating Cgam(i,j) and Csig(i,j) w.r.t. wake.x[i] and wake.t[i].
    //
    // Cgam(i,j) = (a_j(wake.x[i]) · wake.t[i]) (vortex panel j influence at wake point i)
    // d(Cgam(i,j))/d(wake.x[i]) = dV_j/dxi · t_i  (from a_xi/b_xi)
    // d(Cgam(i,j))/d(wake.t[i]) = V_j(wake.x[i])  (velocity vector from panel j)
    //
    // Similarly for Csig (source panel influences at wake points).
    //
    // The ue_m contribution is: -sgnue * ue_sigma * sigma_m * mass
    // At wake station i: R_ue(N+i) += -ue_sigma(N+i,:) * phi
    // where phi = sigma_m * mass, and ue_sigma(N+i,:) = Dw(i,:)
    // Dw = Cgam*Bp + Csig, with Dw(0) = Bp(N-1)
    //
    // For the adjoint through wake geometry, we need:
    //   lambda_ue(N+i) * sgnue(N+i) * d(Dw(i,:)*phi)/d(wake.x[i])
    //   lambda_ue(N+i) * sgnue(N+i) * d(Dw(i,:)*phi)/d(wake.t[i])
    //
    // d(Dw(i,:)*phi)/d(wake.x[i]) has contributions from both Cgam and Csig.
    // Csig at wake points involves wake source panels (frozen geometry) and
    // foil source panels (evaluation point changes). For non-self wake source
    // panels, the Xj endpoints are functions of wake.x, but we account for that
    // in the backward sweep.

    // Compute phi = sigma_m * mass (same as in ComputeDRueDxFoil)
    VectorXd ue_vec = state.glob.U.row(3).transpose();
    VectorXd ds_vec = state.glob.U.row(1).transpose();
    VectorXd mass = ds_vec.cwiseProduct(ue_vec);
    VectorXd phi = state.vsol.sigma_m * mass;

    // Compute Bp_phi = -(AIC^{-1} * B * phi)(0:N-1) for Cgam contribution
    // Actually we need Bp * phi directly. Bp = state.vsol.ue_sigma.topRows(N)
    VectorXd Bp_phi = state.vsol.ue_sigma.topRows(N) * phi;

    // sgnue for wake points is always +1
    // lambda_ue at wake station i is lam_w(i)

    TEInfo te = GetTEInfo(state.foil.x);

    for (int i = 1; i < Nw; ++i) {  // skip i=0 (Dw(0) = Bp(N-1), no wake geom dep)
        double lam_i = lam_w(i);
        if (std::abs(lam_i) < 1e-20) continue;

        Vector2d xi_w = state.wake.x.col(i);
        Vector2d ti_w = state.wake.t.col(i);

        // d(Cgam(i,:)*Bp_phi)/d(wake.x[i]) and d(wake.t[i])
        // Cgam(i,:)*Bp_phi = V_Bp(wake.x[i]) · wake.t[i]
        // where V_Bp is the velocity from treating Bp_phi as gamma strengths

        // Velocity from Bp_phi vortices
        Vector2d V_Bp = InviscidVelocity(state.foil.x, Bp_phi, 0.0, 0.0,
                                          xi_w, false, nullptr);

        // d(V_Bp · t)/d(t) = V_Bp
        mu_t.col(i) -= lam_i * V_Bp;

        // d(V_Bp · t)/d(xi) = dV_Bp/dxi^T * t
        Matrix2d dVBp_dxi = Matrix2d::Zero();
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);

            auto vd = PanelLinVortexVelocityDeriv(Xj, xi_w);
            dVBp_dxi += vd.a_xi * Bp_phi(j) + vd.b_xi * Bp_phi(j + 1);
        }
        // TE panels
        {
            Matrix2d XjTE;
            XjTE.col(0) = state.foil.x.col(N - 1);
            XjTE.col(1) = state.foil.x.col(0);

            auto svd = PanelConstSourceVelocityDeriv(XjTE, xi_w);
            double sigma_te = 0.5 * te.tcp * (Bp_phi(N - 1) - Bp_phi(0));
            dVBp_dxi += svd.val_xi * sigma_te;

            auto vvd = PanelLinVortexVelocityDeriv(XjTE, xi_w);
            double gam_te = 0.5 * te.tdp * (Bp_phi(0) - Bp_phi(N - 1));
            dVBp_dxi += (vvd.a_xi + vvd.b_xi) * gam_te;
        }
        mu_x.col(i) -= lam_i * (dVBp_dxi.transpose() * ti_w);

        // d(Csig(i,:)*phi)/d(wake.x[i]) and d(wake.t[i])
        // Csig has foil source panels evaluated at wake.x[i] projected onto wake.t[i]
        // d(Csig(i,j)*phi(j))/d(wake.x[i]) = d(v_src(j,wake.x[i]) · t) / d(xi) * phi(j)
        // d(Csig(i,j)*phi(j))/d(wake.t[i]) = v_src(j,wake.x[i]) * phi(j)
        int jstart = (i == 0) ? 1 : 0;
        int jend = (i == 0) ? N - 2 : N - 1;
        for (int j = jstart; j < jend; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);

            auto svd = PanelConstSourceVelocityDeriv(Xj, xi_w);

            // val_xi is [2x2]: d(velocity)/d(xi)
            // Csig(i,j) = v · t, so d(Csig)/d(xi) = (dv/dxi)^T * t
            Vector2d dCsig_dxi = svd.val_xi.transpose() * ti_w;
            mu_x.col(i) -= lam_i * dCsig_dxi * phi(j);

            // d(Csig)/d(t) = v_src
            mu_t.col(i) -= lam_i * svd.val * phi(j);
        }

        // Wake source panel contributions to Csig at non-self points
        // These involve wake Xj endpoints and wake.t projection —
        // we'll handle the wake.x dependence of wake source Xj endpoints
        // in the backward sweep (Step 3) since those are functions of other
        // wake.x values, not the evaluation point itself.
        // The evaluation-point xi sensitivity for wake source panels:
        for (int j = 0; j < Nw; ++j) {
            if (i == j) continue;  // self-influence has special formula
            int il = std::max(j - 1, 0);
            int im = j;
            int ir = std::min(j + 1, Nw - 1);

            Vector2d Xl = 0.5 * (state.wake.x.col(il) + state.wake.x.col(im));
            Vector2d Xm = state.wake.x.col(im);
            Vector2d Xr = 0.5 * (state.wake.x.col(im) + state.wake.x.col(ir));
            if (j == Nw - 1) Xr = 2.0 * Xm - Xl;

            if (j == 0) {
                // LinSourceVelocity from [Xm, Xr]
                Matrix2d Xj_src;
                Xj_src.col(0) = Xm;
                Xj_src.col(1) = Xr;
                auto lsvd = PanelLinSourceVelocityDeriv(Xj_src, xi_w);

                // Csig(i, N-1) += b·t, Csig(i, 0) += a·t, Csig(i, N-2) += a·t
                // d(a·t)/dxi = a_xi^T * t, d(b·t)/dxi = b_xi^T * t
                Vector2d da_dxi = lsvd.a_xi.transpose() * ti_w;
                Vector2d db_dxi = lsvd.b_xi.transpose() * ti_w;
                mu_x.col(i) -= lam_i * (da_dxi * (phi(0) + phi(N - 2)) +
                                         db_dxi * phi(N - 1));
                // d/dt
                mu_t.col(i) -= lam_i * (lsvd.a * (phi(0) + phi(N - 2)) +
                                         lsvd.b * phi(N - 1));
            } else if (j == Nw - 1) {
                // ConstSourceVelocity from [Xl, Xr]
                Matrix2d Xj_src;
                Xj_src.col(0) = Xl;
                Xj_src.col(1) = Xr;
                auto csvd = PanelConstSourceVelocityDeriv(Xj_src, xi_w);
                Vector2d dv_dxi = csvd.val_xi.transpose() * ti_w;
                mu_x.col(i) -= lam_i * dv_dxi * phi(N + Nw - 3);
                mu_t.col(i) -= lam_i * csvd.val * phi(N + Nw - 3);
            } else {
                // Two half-panels: [Xl,Xm] and [Xm,Xr]
                Matrix2d Xj_l, Xj_r;
                Xj_l.col(0) = Xl;
                Xj_l.col(1) = Xm;
                Xj_r.col(0) = Xm;
                Xj_r.col(1) = Xr;

                auto lsvd1 = PanelLinSourceVelocityDeriv(Xj_l, xi_w);
                auto lsvd2 = PanelLinSourceVelocityDeriv(Xj_r, xi_w);

                // Csig(i, N-1+j-1) += a1·t + 0.5*b1·t
                // Csig(i, N-1+j) += 0.5*a2·t + b2·t
                double phi_jm1 = phi(N - 1 + j - 1);
                double phi_j = phi(N - 1 + j);

                Vector2d da1_dxi = lsvd1.a_xi.transpose() * ti_w;
                Vector2d db1_dxi = lsvd1.b_xi.transpose() * ti_w;
                Vector2d da2_dxi = lsvd2.a_xi.transpose() * ti_w;
                Vector2d db2_dxi = lsvd2.b_xi.transpose() * ti_w;

                mu_x.col(i) -= lam_i * ((da1_dxi + 0.5 * db1_dxi) * phi_jm1 +
                                         (0.5 * da2_dxi + db2_dxi) * phi_j);
                mu_t.col(i) -= lam_i * ((lsvd1.a + 0.5 * lsvd1.b) * phi_jm1 +
                                         (0.5 * lsvd2.a + lsvd2.b) * phi_j);
            }
        }
    }

    // --- (C) Csig Xj endpoint sensitivity at wake stations ---
    // When wake.x[k] moves, wake source panels adjacent to k change their
    // Xj endpoints. This changes Csig at ALL other wake stations.
    // Only affects mu_x (not mu_t — t[i] at evaluation station i is fixed).
    //
    // Panel j (centered at wake node j) uses:
    //   Xl = 0.5*(wake.x[max(j-1,0)] + wake.x[j])   (left midpoint)
    //   Xm = wake.x[j]                                (center)
    //   Xr = 0.5*(wake.x[j] + wake.x[min(j+1,Nw-1)]) (right midpoint)
    //   if j==Nw-1: Xr = 2*Xm - Xl  (ghost extension)
    //
    // Csig(i, cols) = panel_linsource_velocity(half-panels, wake.x[i], wake.t[i])
    // Moving wake.x[k] changes the half-panel Xj endpoints through the
    // midpoint formulas, which changes Csig(i,:)*phi at all i ≠ j.
    for (int j = 0; j < Nw; ++j) {
        int jl = std::max(j - 1, 0);
        int jr = std::min(j + 1, Nw - 1);
        Vector2d Xl = 0.5 * (state.wake.x.col(jl) + state.wake.x.col(j));
        Vector2d Xm_j = state.wake.x.col(j);
        Vector2d Xr = 0.5 * (state.wake.x.col(j) + state.wake.x.col(jr));
        if (j == Nw - 1) Xr = 2.0 * Xm_j - Xl;

        // Self-influence (i==j): Csig has an analytical formula that
        // depends on d1=|Xm-Xl| and d2=|Xr-Xm|, which are functions of
        // the panel endpoints. d(Csig(j,:)*phi)/d(wake.x[k]) through d1,d2.
        // j=0: Dw(0,:) = Bp(N-1,:) overrides Csig at station 0, so no
        //   self-influence derivative needed.
        // j=Nw-1: Csig self-influence is 0 (ghost extension), no derivative.
        if (j > 0 && j < Nw - 1) {
            double lam_j = lam_w(j);
            if (std::abs(lam_j) > 1e-20) {
                double d1 = (Xm_j - Xl).norm();
                double d2 = (Xr - Xm_j).norm();
                // Csig(j,N-1+j-1) += (0.25/pi)*log(d1/d2) + 0.5/pi
                // Csig(j,N-1+j)   += (0.25/pi)*log(d1/d2) - 0.5/pi
                // Only the log(d1/d2) part depends on endpoints.
                // d(log(d1/d2))/d(wake.x[k]) = d(d1)/d(x[k])/d1 - d(d2)/d(x[k])/d2
                double phi_jm1 = phi(N - 1 + j - 1);
                double phi_jv = phi(N - 1 + j);
                // phi_sum weights the aa term: phi_jm1 + phi_jv
                double phi_aa = phi_jm1 + phi_jv;

                Vector2d e1 = (Xm_j - Xl) / std::max(d1, 1e-20);  // unit vec
                Vector2d e2 = (Xr - Xm_j) / std::max(d2, 1e-20);

                for (int dim = 0; dim < 2; ++dim) {
                    // d(d1)/d(wake.x[k]_dim):
                    //   d1 = |Xm - Xl|, Xm = wake.x[j], Xl = 0.5*(wake.x[j-1]+wake.x[j])
                    //   Xm - Xl = 0.5*(wake.x[j] - wake.x[j-1])
                    //   d(d1)/d(wake.x[j]_dim) = 0.5*e1(dim)
                    //   d(d1)/d(wake.x[j-1]_dim) = -0.5*e1(dim)
                    // d(d2)/d(wake.x[k]_dim):
                    //   d2 = |Xr - Xm|, Xr = 0.5*(wake.x[j]+wake.x[j+1]), Xm = wake.x[j]
                    //   Xr - Xm = 0.5*(wake.x[j+1] - wake.x[j])
                    //   d(d2)/d(wake.x[j+1]_dim) = 0.5*e2(dim)
                    //   d(d2)/d(wake.x[j]_dim) = -0.5*e2(dim)

                    double coeff = (0.25 / M_PI) * phi_aa;

                    // wake.x[j-1]: dd1/d(j-1) = -0.5*e1
                    double daa_jm1 = (-0.5 * e1(dim)) / d1;
                    mu_x(dim, j - 1) -= lam_j * coeff * daa_jm1;

                    // wake.x[j]: dd1/d(j) = +0.5*e1, dd2/d(j) = -0.5*e2
                    double daa_jj = (0.5 * e1(dim)) / d1 - (-0.5 * e2(dim)) / d2;
                    mu_x(dim, j) -= lam_j * coeff * daa_jj;

                    // wake.x[j+1]: dd2/d(j+1) = +0.5*e2
                    double daa_jp1 = -(0.5 * e2(dim)) / d2;
                    mu_x(dim, j + 1) -= lam_j * coeff * daa_jp1;
                }
            }
        }

        for (int i = 1; i < Nw; ++i) {
            if (i == j) continue;
            double lam_i = lam_w(i);
            if (std::abs(lam_i) < 1e-20) continue;
            Vector2d xi_w = state.wake.x.col(i);
            Vector2d ti_w = state.wake.t.col(i);

            if (j == 0) {
                // Panel 0: right half only [Xm, Xr]
                Matrix2d Xj_src;
                Xj_src.col(0) = Xm_j;        // = wake.x[0]
                Xj_src.col(1) = Xr;           // = 0.5*(wake.x[0]+wake.x[1])
                auto lsvd = PanelLinSourceVelocityDeriv(Xj_src, xi_w);

                // Csig(i,0) += a·t, Csig(i,N-2) += a·t, Csig(i,N-1) += b·t
                double phi_a = phi(0) + phi(N - 2);
                double phi_b = phi(N - 1);

                for (int dim = 0; dim < 2; ++dim) {
                    // d(Csig*phi)/d(Xj0_dim) · t + d(Csig*phi)/d(Xj1_dim) · t
                    double dW_Xj0 = lsvd.a_Xj.col(dim).dot(ti_w) * phi_a +
                                    lsvd.b_Xj.col(dim).dot(ti_w) * phi_b;
                    double dW_Xj1 = lsvd.a_Xj.col(2 + dim).dot(ti_w) * phi_a +
                                    lsvd.b_Xj.col(2 + dim).dot(ti_w) * phi_b;

                    // Xj0 = wake.x[0], Xj1 = 0.5*(wake.x[0]+wake.x[1])
                    mu_x(dim, 0) -= lam_i * (dW_Xj0 + 0.5 * dW_Xj1);
                    mu_x(dim, 1) -= lam_i * (0.5 * dW_Xj1);
                }
            } else if (j == Nw - 1) {
                // Last point: const source panel [Xl, Xr]
                Matrix2d Xj_src;
                Xj_src.col(0) = Xl;   // = 0.5*(wake.x[Nw-2]+wake.x[Nw-1])
                Xj_src.col(1) = Xr;   // = 2*wake.x[Nw-1] - Xl
                auto csvd = PanelConstSourceVelocityDeriv(Xj_src, xi_w);

                double phi_last = phi(N + Nw - 3);
                for (int dim = 0; dim < 2; ++dim) {
                    double dW_Xj0 = csvd.val_Xj.col(dim).dot(ti_w) * phi_last;
                    double dW_Xj1 = csvd.val_Xj.col(2 + dim).dot(ti_w) * phi_last;

                    // Xl = 0.5*(wake.x[Nw-2]+wake.x[Nw-1])
                    // Xr = 2*wake.x[Nw-1] - Xl = 1.5*wake.x[Nw-1] - 0.5*wake.x[Nw-2]
                    // d(Xj0)/d(wake.x[Nw-2]) = 0.5,  d(Xj0)/d(wake.x[Nw-1]) = 0.5
                    // d(Xj1)/d(wake.x[Nw-2]) = -0.5, d(Xj1)/d(wake.x[Nw-1]) = 1.5
                    mu_x(dim, Nw - 2) -= lam_i * (0.5 * dW_Xj0 - 0.5 * dW_Xj1);
                    mu_x(dim, Nw - 1) -= lam_i * (0.5 * dW_Xj0 + 1.5 * dW_Xj1);
                }
            } else {
                // Interior panel: left [Xl,Xm] and right [Xm,Xr] halves
                Matrix2d Xj_l, Xj_r;
                Xj_l.col(0) = Xl;
                Xj_l.col(1) = Xm_j;
                Xj_r.col(0) = Xm_j;
                Xj_r.col(1) = Xr;

                auto lsvd1 = PanelLinSourceVelocityDeriv(Xj_l, xi_w);
                auto lsvd2 = PanelLinSourceVelocityDeriv(Xj_r, xi_w);

                double phi_jm1 = phi(N - 1 + j - 1);
                double phi_j_val = phi(N - 1 + j);

                for (int dim = 0; dim < 2; ++dim) {
                    // Left half: Csig(i,N-1+j-1) += a1·t + 0.5*b1·t
                    //            Csig(i,N-1+j) += 0.5*a2·t + b2·t
                    double dW1_Xj0 = (lsvd1.a_Xj.col(dim).dot(ti_w) +
                                      0.5 * lsvd1.b_Xj.col(dim).dot(ti_w)) * phi_jm1;
                    double dW1_Xj1 = (lsvd1.a_Xj.col(2 + dim).dot(ti_w) +
                                      0.5 * lsvd1.b_Xj.col(2 + dim).dot(ti_w)) * phi_jm1;

                    double dW2_Xj0 = (0.5 * lsvd2.a_Xj.col(dim).dot(ti_w) +
                                      lsvd2.b_Xj.col(dim).dot(ti_w)) * phi_j_val;
                    double dW2_Xj1 = (0.5 * lsvd2.a_Xj.col(2 + dim).dot(ti_w) +
                                      lsvd2.b_Xj.col(2 + dim).dot(ti_w)) * phi_j_val;

                    double dW_Xj = dW1_Xj0 + dW1_Xj1 + dW2_Xj0 + dW2_Xj1;
                    (void)dW_Xj;

                    // Xl = 0.5*(wake.x[j-1]+wake.x[j])
                    //   d(Xl)/d(wake.x[j-1]) = 0.5, d(Xl)/d(wake.x[j]) = 0.5
                    // Xm = wake.x[j]
                    //   d(Xm)/d(wake.x[j]) = 1
                    // Xr = 0.5*(wake.x[j]+wake.x[j+1])
                    //   d(Xr)/d(wake.x[j]) = 0.5, d(Xr)/d(wake.x[j+1]) = 0.5

                    // Left half Xj_l = [Xl, Xm]:
                    //   d(Xj_l0)/d(wake.x[j-1]) = 0.5
                    //   d(Xj_l0)/d(wake.x[j]) = 0.5
                    //   d(Xj_l1)/d(wake.x[j]) = 1.0

                    // Right half Xj_r = [Xm, Xr]:
                    //   d(Xj_r0)/d(wake.x[j]) = 1.0
                    //   d(Xj_r1)/d(wake.x[j]) = 0.5
                    //   d(Xj_r1)/d(wake.x[j+1]) = 0.5

                    // Contribution to wake.x[j-1]:
                    double dW_jm1 = 0.5 * dW1_Xj0;  // from left-half Xj0
                    mu_x(dim, j - 1) -= lam_i * dW_jm1;

                    // Contribution to wake.x[j]:
                    double dW_jj = 0.5 * dW1_Xj0 + 1.0 * dW1_Xj1 +
                                   1.0 * dW2_Xj0 + 0.5 * dW2_Xj1;
                    mu_x(dim, j) -= lam_i * dW_jj;

                    // Contribution to wake.x[j+1]:
                    double dW_jp1 = 0.5 * dW2_Xj1;  // from right-half Xj1
                    mu_x(dim, j + 1) -= lam_i * dW_jp1;
                }
            }
        }
    }

    // =========================================================================
    // Step 1b: B matrix wake source panel sensitivity → mu_x
    //
    // ue_m at foil stations depends on Bp = -AIC^{-1}*B. The B matrix has
    // wake source panel contributions that depend on wake.x through the panel
    // endpoints. When wake.x changes, B changes → Bp changes → R_ue changes.
    //
    // VJP: v = sgnue_foil .* lambda_ue_foil + Cgam^T * lam_wake + e_{N-1}*lam_w(0),
    //       w = AIC^{-T} * [v; 0]
    //       mu_wake_x[k] += w^T * d(B*phi)/d(wake.x[k])
    // =========================================================================
    //
    // ue_m at foil stations depends on Bp = -AIC^{-1}*B. The B matrix has
    // wake source panel contributions that depend on wake.x through the panel
    // endpoints. When wake.x changes, B changes → Bp changes → R_ue changes.
    //
    // VJP: v = sgnue_foil .* lambda_ue_foil + Cgam^T * lam_wake,
    //       w = AIC^{-T} * [v; 0]
    //       mu_wake_x[k] += w^T * d(B*phi)/d(wake.x[k])
    // =========================================================================
    {
        // Build the combined adjoint for Bp sensitivity
        VectorXd v_ext = VectorXd::Zero(N + 1);
        VectorXd lambda_ue_foil = lambda_ue.head(N);
        for (int i = 0; i < N; ++i) {
            v_ext(i) = state.isol.sgnue(i) * lambda_ue_foil(i);
        }

        // Wake station contributions through Cgam*Bp:
        for (int i = 1; i < Nw; ++i) {
            double lam_i = lam_w(i);
            if (std::abs(lam_i) < 1e-20) continue;
            Vector2d xi_w = state.wake.x.col(i);
            Vector2d ti_w = state.wake.t.col(i);
            MatrixXd V_G;
            InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi_w,
                             true, &V_G);
            for (int jj = 0; jj < N; ++jj) {
                v_ext(jj) += lam_i * (V_G(0, jj) * ti_w(0) + V_G(1, jj) * ti_w(1));
            }
        }

        // Dw(0) = Bp(N-1) override: wake station 0 uses Bp(N-1,:) directly
        v_ext(N - 1) += lam_w(0);

        MatrixXd AIC_T_b = state.isol.AIC.transpose();
        VectorXd w = AIC_T_b.partialPivLu().solve(v_ext);

        // For each wake panel j, compute the full d(B*phi)/d(wake.x[k]) including
        // ALL B columns that the panel contributes to (N-1+j-1, N-1+j, N-1+j+1).
        for (int j = 0; j < Nw - 1; ++j) {
            Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
            Vector2d X_right;
            if (j == Nw - 2) {
                X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
            } else {
                X_right = state.wake.x.col(j + 1);
            }

            // Phi values for the columns this panel writes to
            double phi_jm1 = (j > 0) ? phi(N - 1 + j - 1) : 0.0;
            double phi_j = phi(N - 1 + j);
            double phi_jp1 = (j < Nw - 2) ? phi(N - 1 + j + 1) : 0.0;

            for (int i = 0; i < N; ++i) {
                if (std::abs(w(i)) < 1e-20) continue;
                Vector2d xi_eval = state.foil.x.col(i);

                Matrix2d Xj_l;
                Xj_l.col(0) = state.wake.x.col(j);
                Xj_l.col(1) = Xm;
                auto ld = PanelLinSourceStreamDeriv(Xj_l, xi_eval);

                Matrix2d Xj_r;
                Xj_r.col(0) = Xm;
                Xj_r.col(1) = X_right;
                auto rd = PanelLinSourceStreamDeriv(Xj_r, xi_eval);

                // Left half: a_l, b_l map to columns:
                //   j>0: col N-1+j gets (0.5*a_l+b_l), col N-1+j-1 gets 0.5*a_l
                //   j==0: col N-1+j gets b_l
                // Right half: a_r, b_r map to columns:
                //   col N-1+j gets (a_r + 0.5*b_r)
                //   j<Nw-2: col N-1+j+1 gets 0.5*b_r
                //   j==Nw-2: col N-1+j gets extra 0.5*b_r

                // Weighted sum: B(i,:)*phi contribution from this panel
                // W = left_col_jm1*phi_jm1 + left_col_j*phi_j + right_col_j*phi_j + right_col_jp1*phi_jp1

                for (int dim = 0; dim < 2; ++dim) {
                    // d(a_l)/d(wake.x[j]_dim) = a_Xj(dim) + 0.5*a_Xj(2+dim)
                    // d(b_l)/d(wake.x[j]_dim) = b_Xj(dim) + 0.5*b_Xj(2+dim)
                    double da_l_dj = ld.a_Xj(dim) + 0.5 * ld.a_Xj(2 + dim);
                    double db_l_dj = ld.b_Xj(dim) + 0.5 * ld.b_Xj(2 + dim);

                    // d(a_r)/d(wake.x[j]_dim) = 0.5*a_Xj(dim) (through Xm in col 0)
                    double da_r_dj = 0.5 * rd.a_Xj(dim);
                    double db_r_dj = 0.5 * rd.b_Xj(dim);
                    if (j == Nw - 2) {
                        da_r_dj += -0.5 * rd.a_Xj(2 + dim);
                        db_r_dj += -0.5 * rd.b_Xj(2 + dim);
                    }

                    // Total d(B*phi)/d(wake.x[j]_dim) from this panel
                    double dWj = 0.0;
                    if (j > 0) {
                        dWj += 0.5 * da_l_dj * phi_jm1;  // col N-1+j-1
                        dWj += (0.5 * da_l_dj + db_l_dj) * phi_j;  // col N-1+j from left
                    } else {
                        dWj += db_l_dj * phi_j;  // col N-1+j from left
                    }
                    dWj += (da_r_dj + 0.5 * db_r_dj) * phi_j;  // col N-1+j from right
                    if (j < Nw - 2) {
                        dWj += 0.5 * db_r_dj * phi_jp1;  // col N-1+j+1
                    } else {
                        dWj += 0.5 * db_r_dj * phi_j;  // ghost fold-back
                    }

                    mu_x(dim, j) += w(i) * dWj;

                    // d/d(wake.x[j+1]_dim)
                    double da_l_djp1 = 0.5 * ld.a_Xj(2 + dim);
                    double db_l_djp1 = 0.5 * ld.b_Xj(2 + dim);

                    double da_r_djp1 = 0.5 * rd.a_Xj(dim);
                    double db_r_djp1 = 0.5 * rd.b_Xj(dim);
                    if (j < Nw - 2) {
                        da_r_djp1 += rd.a_Xj(2 + dim);
                        db_r_djp1 += rd.b_Xj(2 + dim);
                    } else {
                        da_r_djp1 += 1.5 * rd.a_Xj(2 + dim);
                        db_r_djp1 += 1.5 * rd.b_Xj(2 + dim);
                    }

                    double dWjp1 = 0.0;
                    if (j > 0) {
                        dWjp1 += 0.5 * da_l_djp1 * phi_jm1;
                        dWjp1 += (0.5 * da_l_djp1 + db_l_djp1) * phi_j;
                    } else {
                        dWjp1 += db_l_djp1 * phi_j;
                    }
                    dWjp1 += (da_r_djp1 + 0.5 * db_r_djp1) * phi_j;
                    if (j < Nw - 2) {
                        dWjp1 += 0.5 * db_r_djp1 * phi_jp1;
                    } else {
                        dWjp1 += 0.5 * db_r_djp1 * phi_j;
                    }

                    mu_x(dim, j + 1) += w(i) * dWjp1;
                }
            }
        }
    }


    // Save mu_x/mu_t before backward sweep (for validation)
    if (mu_x_out) *mu_x_out = mu_x;
    if (mu_t_out) *mu_t_out = mu_t;

    // =========================================================================
    // Steps 2-3: Backward sweep through predictor-corrector wake tracing.
    //
    // Propagates mu_x[i] and mu_t[i] backward through the streamline tracing
    // to foil.x. At each step, velocity evaluation-point gradients and panel
    // geometry changes contribute to the foil.x VJP.
    // =========================================================================

    TEInfo te_wake = GetTEInfo(state.foil.x);
    double alpha = state.oper.alpha;

    // Recompute wake spacing (same as BuildWake)
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                         (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double chord_val = state.geom.chord;
    VectorXd sv(Nw);
    {
        int Np = Nw;
        double L = state.geom.wakelen * chord_val;
        if (Np <= 1) {
            sv(0) = 0;
        } else {
            int Ns = Np - 1;
            double dd = L / ds1;
            double a2 = Ns * (Ns - 1.0) * (Ns - 2.0) / 6.0;
            double b2 = Ns * (Ns - 1.0) / 2.0;
            double c2 = Ns - dd;
            double disc = std::max(b2 * b2 - 4 * a2 * c2, 0.0);
            double r = 1.0 + (-b2 + std::sqrt(disc)) / (2 * a2);
            for (int it = 0; it < 10; ++it) {
                double rN = std::pow(r, Ns);
                double R = rN - 1 - dd * (r - 1);
                double Rr = Ns * std::pow(r, Ns - 1) - dd;
                double dr = -R / Rr;
                if (std::abs(dr) < 1e-6) break;
                r += dr;
            }
            sv(0) = 0;
            double dxs = ds1;
            for (int ii = 1; ii < Np; ++ii) {
                sv(ii) = sv(ii - 1) + dxs;
                dxs *= r;
            }
        }
    }

    // Handle wake.t[Nw-1] = wake.t[Nw-2] (copy in forward pass)
    if (Nw >= 2) {
        mu_t.col(Nw - 2) += mu_t.col(Nw - 1);
        mu_t.col(Nw - 1).setZero();
    }

    // Backward sweep: i = Nw-1 down to 1
    for (int i = Nw - 1; i >= 1; --i) {
        double ds = sv(i) - sv(i - 1);
        Vector2d x0 = state.wake.x.col(i - 1);

        // Recompute forward quantities at this step
        Vector2d V0 = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, x0);
        double V0m = V0.norm();
        if (V0m < 1e-10) V0m = 1e-10;
        Vector2d t0 = V0 / V0m;

        Vector2d x1_pred = x0 + ds * t0;

        Vector2d V1 = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, x1_pred);
        double V1m = V1.norm();
        if (V1m < 1e-10) V1m = 1e-10;
        Vector2d t1 = V1 / V1m;

        Vector2d t_sum = 0.5 * (t0 + t1);
        double tmag = t_sum.norm();
        if (tmag < 1e-10) tmag = 1e-10;
        Vector2d t_avg = t_sum / tmag;

        // --- Adjoint propagation ---
        Vector2d mu_xi = mu_x.col(i);
        Vector2d mu_ti = mu_t.col(i - 1);

        // From wake.x[i] = x0 + ds * t_avg:
        mu_x.col(i - 1) += mu_xi;             // direct x0 dependence
        Vector2d mu_tavg = ds * mu_xi;         // through t_avg

        // From wake.t[i-1] = t_avg:
        mu_tavg += mu_ti;

        // NormalizeVJP: t_avg = t_sum / |t_sum|
        // Adjoint: mu_tsum = (I - t_avg*t_avg^T)/tmag * mu_tavg
        Vector2d mu_tsum =
            (mu_tavg - t_avg * t_avg.dot(mu_tavg)) / tmag;

        // t_sum = 0.5*(t0 + t1)
        Vector2d mu_t0 = 0.5 * mu_tsum;
        Vector2d mu_t1 = 0.5 * mu_tsum;

        // NormalizeVJP: t1 = V1 / V1m
        Vector2d mu_V1 = (mu_t1 - t1 * t1.dot(mu_t1)) / V1m;

        // V1 at x1_pred: evaluation point gradient
        Matrix2d G1 = VelocityGradient(state, x1_pred, te_wake);
        Vector2d mu_x1pred = G1.transpose() * mu_V1;

        // V1 panel geometry VJP
        PanelVelocityVJP(state, x1_pred, mu_V1, te_wake, vjp);

        // From x1_pred = x0 + ds * t0:
        mu_x.col(i - 1) += mu_x1pred;    // direct x0 dependence
        mu_t0 += ds * mu_x1pred;          // through t0

        // NormalizeVJP: t0 = V0 / V0m
        Vector2d mu_V0 = (mu_t0 - t0 * t0.dot(mu_t0)) / V0m;

        // V0 at x0: evaluation point gradient
        Matrix2d G0 = VelocityGradient(state, x0, te_wake);
        mu_x.col(i - 1) += G0.transpose() * mu_V0;

        // V0 panel geometry VJP
        PanelVelocityVJP(state, x0, mu_V0, te_wake, vjp);
    }

    // Handle wake.t[0]: computed separately as V(wake.x[0])/|V(wake.x[0])|
    {
        Vector2d x0 = state.wake.x.col(0);
        Vector2d Vt = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, x0);
        double Vtm = Vt.norm();
        if (Vtm < 1e-10) Vtm = 1e-10;
        Vector2d nt = Vt / Vtm;

        // NormalizeVJP: wake.t[0] = Vt / |Vt|
        Vector2d mu_Vt =
            (mu_t.col(0) - nt * nt.dot(mu_t.col(0))) / Vtm;

        // Evaluation point gradient
        Matrix2d Gt = VelocityGradient(state, x0, te_wake);
        mu_x.col(0) += Gt.transpose() * mu_Vt;

        // Panel geometry VJP
        PanelVelocityVJP(state, x0, mu_Vt, te_wake, vjp);
    }

    // TE initial point: wake.x[0] = 0.5*(foil.x[0]+foil.x[N-1]) + eps*t_te
    {
        Vector2d n_te = state.foil.x.col(N - 1) - state.foil.x.col(0);
        Vector2d t_te_raw(n_te(1), -n_te(0));
        double sign_te = (t_te_raw(0) < 0) ? -1.0 : 1.0;
        double eps_wake = 1e-5;

        Matrix2d dwx0_dxy1 = 0.5 * Matrix2d::Identity();
        dwx0_dxy1(0, 1) += eps_wake * chord_val * (-sign_te);
        dwx0_dxy1(1, 0) += eps_wake * chord_val * sign_te;

        Matrix2d dwx0_dxyN = 0.5 * Matrix2d::Identity();
        dwx0_dxyN(0, 1) += eps_wake * chord_val * sign_te;
        dwx0_dxyN(1, 0) += eps_wake * chord_val * (-sign_te);

        Vector2d mu0 = mu_x.col(0);
        vjp(0) += mu0.dot(dwx0_dxy1.col(0));
        vjp(1) += mu0.dot(dwx0_dxy1.col(1));
        vjp(2 * (N - 1)) += mu0.dot(dwx0_dxyN.col(0));
        vjp(2 * (N - 1) + 1) += mu0.dot(dwx0_dxyN.col(1));
    }

    return vjp;
}

// Public overloads
VectorXd ComputeWakeGeomVJP(FoilState& state, const VectorXd& lambda_ue) {
    return ComputeWakeGeomVJPImpl(state, lambda_ue, nullptr, nullptr);
}

VectorXd ComputeWakeGeomVJP(FoilState& state, const VectorXd& lambda_ue,
                             MatrixXd& mu_x_out, MatrixXd& mu_t_out) {
    return ComputeWakeGeomVJPImpl(state, lambda_ue, &mu_x_out, &mu_t_out);
}

// ============================================================================
// ComputeXiVJP: VJP through stagnation point location
//
// xi = |s - sstag| for foil points, xi = wake.s - sstag for wake points
// sstag depends on gamma through linear interpolation: sstag = f(gam(I0), gam(I1))
// gamma depends on foil.x through the AIC system.
//
// Returns lambda_bl^T * R_x * (∂xi/∂sstag · dsstag/dx_foil) [2*N]
// ============================================================================
VectorXd ComputeXiVJP(FoilState& state, const VectorXd& lambda_bl) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    int Nw = state.wake.N;
    double alpha = state.oper.alpha;

    VectorXd vjp = VectorXd::Zero(2 * N);

    // Compute c = lambda_bl^T * R_x * x_st
    // where x_st = ∂xi/∂sstag [Nsys]
    // For foil points: ∂xi(i)/∂sstag = -sgnue(i) (since xi = |s - sstag|)
    // For wake points: ∂xi(N+i)/∂sstag = -1 (since xi = wake.s - sstag)
    VectorXd dxi_dsstag = VectorXd::Zero(Nsys);
    for (int i = 0; i < N; ++i) {
        dxi_dsstag(i) = -state.isol.sgnue(i);
    }
    for (int i = 0; i < Nw; ++i) {
        dxi_dsstag(N + i) = -1.0;
    }

    // R_x is [3*Nsys x Nsys], dxi_dsstag is [Nsys]
    // c = lambda_bl^T * R_x * dxi_dsstag (scalar)
    double c = lambda_bl.dot(state.glob.R_x * dxi_dsstag);

    if (std::abs(c) < 1e-20) return vjp;

    // Now chain through sstag → gamma → foil.x
    // sstag is linearly interpolated between gammas at Istag[0] and Istag[1]
    // In the viscous solver (StagPointMove):
    //   ue0, ue1 are the viscous edge velocities at Istag nodes
    //   sstag = (ue1*s0 + ue0*s1) / (ue0 + ue1)
    //   sstag_ue = [ue1*(s1-s0)/(ue0+ue1)^2, -ue0*(s1-s0)/(ue0+ue1)^2]
    //
    // The ue values at stagnation are in turn functions of the BL state,
    // but for coordinate sensitivities, the stagnation location is
    // determined by where gamma = 0, which depends on gamma through the AIC.
    //
    // Use the initial stagnation formula (from StagPointFind):
    //   sstag depends on gam(I0) and gam(I1) through linear interpolation:
    //   sstag_g = [sstag_g(0), sstag_g(1)] = d(sstag)/d(gam) at Istag nodes
    //
    // VJP: d_sstag = c (scalar); adjoint in gamma space has 2 nonzero entries:
    //   d[Istag[0]] = c * sstag_g[0]
    //   d[Istag[1]] = c * sstag_g[1]
    //
    // Then solve AIC^T * mu = d and compute AIC VJP.

    int I0 = state.isol.Istag[0];
    int I1 = state.isol.Istag[1];

    VectorXd d = VectorXd::Zero(N + 1);
    d(I0) = c * state.isol.sstag_g(0);
    d(I1) = c * state.isol.sstag_g(1);

    // Solve AIC^T * mu = d
    MatrixXd AIC_T = state.isol.AIC.transpose();
    VectorXd mu = AIC_T.partialPivLu().solve(d);

    // Compute VJP through AIC residual
    vjp = AICResidualSensVJP(state, mu, alpha);

    return vjp;
}

// ============================================================================
// AlphaSensitivities: compute dJ/dalpha for cl, cd, cm
// ============================================================================
void AlphaSensitivities(FoilState& state, Sensitivity& sens) {
    int N = state.foil.N;
    double alpha = state.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    // dgam/dalpha = gamref[:,0]*(-sin a) + gamref[:,1]*(cos a), times deg2rad
    VectorXd dgam_dalpha = (-sina * state.isol.gamref.col(0) +
                             cosa * state.isol.gamref.col(1)) * kDegToRad;

    if (!state.oper.viscous) {
        // --- Inviscid path ---
        // Force rotation partials are already in post: cl_alpha, cdpi_alpha
        // cm_alpha = 0 (moment arm is alpha-independent when xref is fixed)
        //
        // dueinv/dalpha = sgnue .* dgam/dalpha
        // dJ/dalpha = J_alpha_rotation + J_ue * dueinv/dalpha
        VectorXd dueinv_dalpha = state.isol.sgnue.cwiseProduct(dgam_dalpha);

        sens.dcl_dalpha = state.post.cl_alpha +
                          state.post.cl_ue.dot(dueinv_dalpha);
        sens.dcm_dalpha = 0.0 + state.post.cm_ue.dot(dueinv_dalpha);
        sens.dcd_dalpha = state.post.cdpi_alpha +
                          state.post.cdpi_ue.dot(dueinv_dalpha);
    } else {
        // --- Viscous path ---
        // Need adjoint vectors from stored LU factorization
        if (!state.glob.lu_factorization) {
            throw std::runtime_error(
                "AlphaSensitivities: no LU factorization stored");
        }

        int Nsys = state.glob.Nsys;
        int Nw = state.wake.N;

        // Ensure force linearizations are computed
        CalcForceLin(state);
        CalcCdLin(state);

        // Build dJ/dV for each objective (interleaved column ordering)
        VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
        VectorXd dcm_dV = VectorXd::Zero(4 * Nsys);
        VectorXd dcd_dV = state.post.cd_V;
        for (int i = 0; i < N; ++i) {
            dcl_dV(4 * i + 3) = state.post.cl_ue(i);
            dcm_dV(4 * i + 3) = state.post.cm_ue(i);
        }

        // Solve adjoint: R_V^T * lambda = dJ/dV
        // lambda has block row ordering: [BL(3*Nsys) | ue(Nsys)]
        auto& lu = *state.glob.lu_factorization;
        VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
        VectorXd lambda_cm = lu.transpose().solve(dcm_dV);
        VectorXd lambda_cd = lu.transpose().solve(dcd_dV);

        // Compute dR_ue/dalpha [Nsys]
        // BL residuals don't depend on alpha directly: dR_BL/dalpha = 0
        // 4th equation: R_ue = ue - ueinv - ue_m*mass
        // dR_ue/dalpha = -dueinv/dalpha
        VectorXd dRue_dalpha = VectorXd::Zero(Nsys);

        // Airfoil: dueinv/dalpha = sgnue .* dgam/dalpha
        VectorXd dueinv_dalpha = state.isol.sgnue.cwiseProduct(dgam_dalpha);
        dRue_dalpha.head(N) = -dueinv_dalpha;

        // Wake: dueinv_wake/dalpha = uewiref * [-sin a, cos a]^T * deg2rad
        if (Nw > 0 && state.isol.uewiref.rows() > 0) {
            Vector2d dcs(-sina * kDegToRad, cosa * kDegToRad);
            VectorXd dueinv_wake_dalpha = state.isol.uewiref * dcs;
            // First wake point = last airfoil point (continuity)
            dueinv_wake_dalpha(0) = dueinv_dalpha(N - 1);
            dRue_dalpha.tail(Nw) = -dueinv_wake_dalpha;
        }

        // R_V has block row ordering: [BL(3*Nsys) | ue(Nsys)]
        // lambda from R_V^T solve has same ordering
        // dR/dalpha = [0_{3*Nsys}; dRue_dalpha_{Nsys}] (block ordered)
        // So lambda^T * dR/dalpha = lambda_ue^T * dRue_dalpha
        VectorXd lambda_cl_ue = lambda_cl.tail(Nsys);
        VectorXd lambda_cm_ue = lambda_cm.tail(Nsys);
        VectorXd lambda_cd_ue = lambda_cd.tail(Nsys);

        // Partial force derivatives w.r.t. alpha (rotation terms)
        // cl: cl_alpha from force integration rotation
        // cm: 0 (moment arm is alpha-independent)
        // cd: 0 (Squire-Young drag has no direct alpha dependence)
        double dcl_dalpha_partial = state.post.cl_alpha;
        double dcm_dalpha_partial = 0.0;
        double dcd_dalpha_partial = 0.0;

        // Total: dJ/dalpha = partial - lambda_ue^T * dRue/dalpha
        sens.dcl_dalpha = dcl_dalpha_partial - lambda_cl_ue.dot(dRue_dalpha);
        sens.dcm_dalpha = dcm_dalpha_partial - lambda_cm_ue.dot(dRue_dalpha);
        sens.dcd_dalpha = dcd_dalpha_partial - lambda_cd_ue.dot(dRue_dalpha);
    }
}

// ============================================================================
// ComputeWgapVJP: VJP through wake gap (wgap) dependence on foil.x
//
// The BL residual at wake stations subtracts wgap from delta*:
//   U.row(1) -= Aux, where Aux = wgap at wake stations
//
// wgap depends on foil.x through hTE and dtdx (GetTEInfo), which depend
// on nodes {0, 1, N-2, N-1} only.
//
// VJP: lambda_bl^T * dR_BL/d(wgap) * d(wgap)/d(foil.x)
// ============================================================================
VectorXd ComputeWgapVJP(FoilState& state, const VectorXd& lambda_bl) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    int Nw = state.wake.N;

    VectorXd vjp = VectorXd::Zero(2 * N);
    if (Nw <= 0) return vjp;

    // Recompute wgap parameters
    TEInfo te = GetTEInfo(state.foil.x);
    double hTE = te.hTE;
    double flen = 2.5;
    double dtdx = std::min(std::max(te.dtdx, -3.0 / flen), 3.0 / flen);
    bool dtdx_clipped = (te.dtdx < -3.0 / flen || te.dtdx > 3.0 / flen);
    double Lw = flen * hTE;
    double xi0 = state.isol.xi(N);
    double alpha = 2.0 + flen * dtdx;

    // Step 1: Compute mu_wgap(w) = -lambda_bl^T * R_V_col(4*(N+w)+1)
    // dR_BL/d(wgap_w) = -dR_BL/d(ds_w) since U.row(1) -= Aux
    const auto& RV = state.glob.R_V;
    VectorXd mu_wgap = VectorXd::Zero(Nw);
    for (int w = 0; w < Nw; ++w) {
        int col = 4 * (N + w) + 1;  // ds column for wake station w
        double dot = 0.0;
        for (SparseMatrixd::InnerIterator it(RV, col); it; ++it) {
            if (it.row() < 3 * Nsys) {
                dot += lambda_bl(it.row()) * it.value();
            }
        }
        mu_wgap(w) = -dot;
    }

    // Step 2: Assemble VJP through d(wgap)/d(foil.x)
    // Only nodes {0, 1, N-2, N-1} affect hTE and dtdx
    int te_nodes[] = {0, 1, N - 2, N - 1};
    for (int k : te_nodes) {
        for (int dim = 0; dim < 2; ++dim) {
            double dhTE_k, ddtdx_k;
            GetTEGapDeriv(state.foil.x, k, dim, dhTE_k, ddtdx_k);

            if (dtdx_clipped) ddtdx_k = 0.0;

            double sum = 0.0;
            for (int w = 0; w < Nw; ++w) {
                double xi_wake = state.isol.xi(N + w);
                double xib = (xi_wake - xi0) / std::max(Lw, 1e-10);

                if (xib >= 0.0 && xib <= 1.0) {
                    double omxib = 1.0 - xib;

                    // dwgap/dhTE = (1-xib) * (1 + xib + 2*alpha*xib^2)
                    double dwgap_dhTE =
                        omxib * (1.0 + xib + 2.0 * alpha * xib * xib);

                    // dwgap/ddtdx = hTE * flen * xib * (1-xib)^2
                    double dwgap_ddtdx = hTE * flen * xib * omxib * omxib;

                    double dwgap_k =
                        dwgap_dhTE * dhTE_k + dwgap_ddtdx * ddtdx_k;
                    sum += mu_wgap(w) * dwgap_k;
                }
            }
            vjp(2 * k + dim) += sum;
        }
    }

    return vjp;
}

// ============================================================================
// ViscousSensitivities: compute dJ/dx using the viscous adjoint method
// ============================================================================
Sensitivity ViscousSensitivities(FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;

    if (!state.glob.lu_factorization) {
        throw std::runtime_error(
            "ViscousSensitivities: no LU factorization stored");
    }

    // Step 1: Compute force partial derivatives
    CalcForceLin(state);
    CalcCdLin(state);

    // Step 2: Form dJ/dV [4*Nsys] for each objective
    // V = [th, ds, sa, ue] at each node, stored column-major as 4*Nsys vector
    // cl depends on ue at airfoil nodes: dcl/dV has entries at indices 4*i+3
    VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
    VectorXd dcm_dV = VectorXd::Zero(4 * Nsys);
    VectorXd dcd_dV = state.post.cd_V;  // Already computed in CalcCdLin

    for (int i = 0; i < N; ++i) {
        dcl_dV(4 * i + 3) = state.post.cl_ue(i);
        dcm_dV(4 * i + 3) = state.post.cm_ue(i);
    }

    // Step 3: Solve adjoint equations: R_V^T * lambda = -(dJ/dV)
    // Note sign: R_V * dV = -R, so the adjoint of the Newton system
    // at convergence (R=0) gives us R_V^T * lambda = (dJ/dV)^T
    // But for the total derivative: dJ/dx = dJ/dx_partial + lambda^T * dR/dx
    // where R_V^T * lambda = (dJ/dV)^T
    auto& lu = *state.glob.lu_factorization;

    // Eigen SparseLU supports transpose solve
    VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
    VectorXd lambda_cm = lu.transpose().solve(dcm_dV);
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);

    // Step 4: Assemble dR/dx_foil [4*Nsys x 2*N]
    // Last Nsys rows: 4th equation
    MatrixXd dRue_dx = ComputeDRueDxFoil(state);

    // Compute lambda^T * dR/dx for each objective
    // Split lambda into BL part (first 3*Nsys) and ue part (last Nsys)
    Sensitivity sens;
    sens.dcl_dx = MatrixXd::Zero(2, N);
    sens.dcm_dx = MatrixXd::Zero(2, N);
    sens.dcd_dx = MatrixXd::Zero(2, N);
    sens.dcdpi_dx = MatrixXd::Zero(2, N);

    auto computeVJP = [&](const VectorXd& lambda,
                          const char* label) -> VectorXd {
        VectorXd lambda_bl = lambda.head(3 * Nsys);
        VectorXd lambda_ue = lambda.tail(Nsys);

        // Arc-length chain: foil.x → foil.s → {xi, sigma_m}
        VectorXd vjp_bl = ComputeArcLengthVJP(state, lambda_bl, lambda_ue);

        // 4th equation contribution: lambda_ue^T * dRue/dx (wake-fixed terms)
        VectorXd vjp_ue = dRue_dx.transpose() * lambda_ue;

        // Wake geometry chain rule: lambda_ue^T * (dR_ue/dwake · dwake/dx_foil)
        VectorXd vjp_wake = ComputeWakeGeomVJP(state, lambda_ue);

        // Wake gap chain rule: lambda_bl^T * dR_BL/dwgap * dwgap/dx_foil
        VectorXd vjp_wgap = ComputeWgapVJP(state, lambda_bl);

        // Diagnostic: print VJP components at test nodes
        int diag_nodes[] = {20, 50, 80};
        for (int kk : diag_nodes) {
            if (kk >= N) continue;
            for (int dd = 0; dd < 2; ++dd) {
                int idx = 2 * kk + dd;
                std::fprintf(stderr,
                    "VJP_%s n=%d d=%d: bl=%.4e ue=%.4e wake=%.4e wgap=%.4e "
                    "partial=%.4e\n",
                    label, kk, dd,
                    vjp_bl(idx), vjp_ue(idx), vjp_wake(idx), vjp_wgap(idx),
                    (std::strcmp(label, "cl") == 0
                         ? state.post.dcl_dx_partial(dd, kk)
                     : std::strcmp(label, "cd") == 0
                         ? state.post.dcd_dx_partial(dd, kk)
                         : 0.0));
            }
        }

        // Xi chain rule: foil.x → AIC → gamma → sstag → xi → R_BL
        VectorXd vjp_xi = ComputeXiVJP(state, lambda_bl);

        return vjp_bl + vjp_ue + vjp_wake + vjp_wgap + vjp_xi;
    };

    // Debug: dump lambda vectors for external validation
    {
        auto dump = [](const char* name, const VectorXd& v) {
            FILE* fp = std::fopen(name, "wb");
            if (fp) {
                std::fwrite(v.data(), sizeof(double), v.size(), fp);
                std::fclose(fp);
                std::fprintf(stderr, "Dumped %s (%d doubles)\n",
                             name, (int)v.size());
            }
        };
        dump("/tmp/lambda_cl.bin", lambda_cl);
        dump("/tmp/lambda_cd.bin", lambda_cd);
    }

    VectorXd vjp_cl = computeVJP(lambda_cl, "cl");
    VectorXd vjp_cm = computeVJP(lambda_cm, "cm");
    VectorXd vjp_cd = computeVJP(lambda_cd, "cd");

    // Step 5: Total derivative = partial - VJP
    // dJ/dx = ∂J/∂x - λ^T * ∂R/∂x  (implicit function theorem)
    for (int k = 0; k < N; ++k) {
        sens.dcl_dx(0, k) = state.post.dcl_dx_partial(0, k) - vjp_cl(2 * k);
        sens.dcl_dx(1, k) = state.post.dcl_dx_partial(1, k) - vjp_cl(2 * k + 1);
        sens.dcm_dx(0, k) = state.post.dcm_dx_partial(0, k) - vjp_cm(2 * k);
        sens.dcm_dx(1, k) = state.post.dcm_dx_partial(1, k) - vjp_cm(2 * k + 1);
        sens.dcd_dx(0, k) = state.post.dcd_dx_partial(0, k) - vjp_cd(2 * k);
        sens.dcd_dx(1, k) = state.post.dcd_dx_partial(1, k) - vjp_cd(2 * k + 1);
    }

    return sens;
}

// ============================================================================
// Mach and Reynolds Sensitivity Infrastructure
// ============================================================================

// Derivative of thermodynamic parameters w.r.t. a single "param perturbation"
// For Mach: dMinf=1, drho0 and dmu0 computed from init_thermo chain rule
// For Re: dMinf=0, drho0=0, dmu0=-mu0/Re
struct ParamDeriv {
    double dMinf;
    double drho0;
    double dmu0;
};

// Compute ParamDeriv for Mach sensitivity (differentiate InitThermo w.r.t. Ma)
static ParamDeriv InitThermoDeriv_Ma(const Param& param, double Ma, double Re,
                                     double rhoinf, double Vinf, double chord) {
    ParamDeriv pd;
    pd.dMinf = 1.0;  // dMinf/dMa = 1

    double gam = param.gam;
    double gmi = gam - 1.0;

    if (Ma > 0) {
        // drho0/dMa: rho0 = rhoinf * (1 + 0.5*gmi*Ma^2)^(1/gmi)
        // drho0/dMa = rhoinf * (1/gmi) * (1 + 0.5*gmi*Ma^2)^(1/gmi - 1) * gmi * Ma
        //           = rhoinf * Ma * (1 + 0.5*gmi*Ma^2)^(1/gmi - 1)
        double fac = 1.0 + 0.5 * gmi * Ma * Ma;
        pd.drho0 = rhoinf * Ma * std::pow(fac, 1.0 / gmi - 1.0);

        // dmu0/dMa: mu0 = muinf / finf
        // muinf = rhoinf * Vinf * chord / Re (independent of Ma)
        // finf = Tr^1.5 * (1 + Tsrat) / (Tr + Tsrat)
        // Tr = 1 - 0.5 * Vinf^2 / H0
        // H0 = (1 + 0.5*gmi*Ma^2) * Vinf^2 / (gmi * Ma^2)
        double muinf = rhoinf * Vinf * chord / Re;

        double H0 = (1.0 + 0.5 * gmi * Ma * Ma) * Vinf * Vinf / (gmi * Ma * Ma);
        // dH0/dMa = d/dMa[ Vinf^2/(gmi*Ma^2) + 0.5*Vinf^2/gmi ]
        //         = Vinf^2/(gmi) * d/dMa[1/Ma^2] = Vinf^2/(gmi) * (-2/Ma^3)
        //         = -2*Vinf^2 / (gmi * Ma^3)
        double dH0 = -2.0 * Vinf * Vinf / (gmi * Ma * Ma * Ma);

        double Tr = 1.0 - 0.5 * Vinf * Vinf / H0;
        // dTr/dMa = 0.5 * Vinf^2 / H0^2 * dH0
        double dTr = 0.5 * Vinf * Vinf / (H0 * H0) * dH0;

        double Tsrat = param.Tsrat;
        double finf = std::pow(Tr, 1.5) * (1.0 + Tsrat) / (Tr + Tsrat);
        // dfinf/dTr = 1.5*Tr^0.5*(1+Tsrat)/(Tr+Tsrat) - Tr^1.5*(1+Tsrat)/(Tr+Tsrat)^2
        //           = finf * (1.5/Tr - 1/(Tr+Tsrat))
        double dfinf_dTr = finf * (1.5 / Tr - 1.0 / (Tr + Tsrat));
        double dfinf = dfinf_dTr * dTr;

        // mu0 = muinf / finf => dmu0/dMa = -muinf / finf^2 * dfinf
        pd.dmu0 = -muinf / (finf * finf) * dfinf;
    } else {
        // Ma=0: rho0 = rhoinf, mu0 = muinf (both constant)
        pd.drho0 = 0.0;
        pd.dmu0 = 0.0;
    }

    return pd;
}

// ============================================================================
// Leaf-level dMa companions (Level 0)
// ============================================================================

// dM2/dMa: M2 = ue^2 * Minf^2 / Vinf^2
// dM2/dMa = 2*M2/Minf * dMinf (since dMinf/dMa = 1 for Ma sensitivities)
static double GetMach2_dMa(const Vector4d& U, const Param& param,
                           const ParamDeriv& pd) {
    if (param.Minf <= 0) return 0.0;
    auto [M2, M2_U] = GetMach2(U, param);
    return 2.0 * M2 / param.Minf * pd.dMinf;
}

// dRet/dMa: Ret = rho0 * |ue| * th / mu0
// dRet/dMa = Ret/rho0 * drho0 - Ret/mu0 * dmu0
static double GetRet_dMa(const Vector4d& U, const Param& param,
                         const ParamDeriv& pd) {
    auto [Ret, Ret_U] = GetRet(U, param);
    double dRet = 0.0;
    if (param.rho0 > 0) dRet += Ret / param.rho0 * pd.drho0;
    if (param.mu0 > 0) dRet -= Ret / param.mu0 * pd.dmu0;
    return dRet;
}

// ============================================================================
// Level 1: GetHk_dMa
// ============================================================================
static double GetHk_dMa(const Vector4d& U, const Param& param,
                        const ParamDeriv& pd) {
    if (param.Minf <= 0) return 0.0;

    auto [H, H_U] = GetH(U);
    auto [M2, M2_U] = GetMach2(U, param);
    double dM2 = GetMach2_dMa(U, param, pd);

    // Hk = (H - 0.29*M2) / (1 + 0.113*M2)
    double den = 1.0 + 0.113 * M2;
    double Hk = (H - 0.29 * M2) / den;

    // dHk/dMa = (-0.29*dM2)/den - Hk/den * 0.113*dM2
    return (-0.29 * dM2 - Hk * 0.113 * dM2) / den;
}

// ============================================================================
// Level 2: GetHs_dMa, GetHss_dMa, GetDe_dMa, GetDamp_dMa
// ============================================================================
static double GetHs_dMa(const Vector4d& U, const Param& param,
                        const ParamDeriv& pd) {
    auto [Hk, Hk_U] = GetHk(U, param);
    double dHk = GetHk_dMa(U, param, pd);

    // Limit Hk (matching primal)
    if (param.wake && Hk < 1.00005) { dHk = 0.0; Hk = 1.00005; }
    if (!param.wake && Hk < 1.05) { dHk = 0.0; Hk = 1.05; }

    double dHs;

    if (param.turb) {
        auto [Ret, Ret_U] = GetRet(U, param);
        double dRet = GetRet_dMa(U, param, pd);

        double Ho = 4.0;
        double dHo = 0.0;
        if (Ret > 400.0) {
            Ho = 3.0 + 400.0 / Ret;
            dHo = -400.0 / (Ret * Ret) * dRet;
        }

        double Reb = Ret;
        double dReb = dRet;
        if (Ret < 200.0) { Reb = 200.0; dReb = 0.0; }

        double Hs_val;
        if (Hk < Ho) {
            double Hr = (Ho - Hk) / (Ho - 1.0);
            double dHr = (dHo - dHk) / (Ho - 1.0) -
                         (Ho - Hk) / ((Ho - 1.0) * (Ho - 1.0)) * dHo;

            double aa = (2.0 - 1.5 - 4.0 / Reb) * Hr * Hr;
            double daa = (4.0 / (Reb * Reb) * dReb) * Hr * Hr +
                         (2.0 - 1.5 - 4.0 / Reb) * 2.0 * Hr * dHr;

            Hs_val = 1.5 + 4.0 / Reb + aa * 1.5 / (Hk + 0.5);
            dHs = -4.0 / (Reb * Reb) * dReb + daa * 1.5 / (Hk + 0.5) -
                  aa * 1.5 / ((Hk + 0.5) * (Hk + 0.5)) * dHk;
        } else {
            double lrb = std::log(Reb);
            double dlrb = dReb / Reb;

            double aa = Hk - Ho + 4.0 / lrb;
            double daa = dHk - dHo - 4.0 / (lrb * lrb) * dlrb;

            double bb = 0.007 * lrb / (aa * aa) + 0.015 / Hk;
            double dbb = 0.007 * (dlrb / (aa * aa) -
                                  2.0 * lrb / (aa * aa * aa) * daa) -
                         0.015 / (Hk * Hk) * dHk;

            Hs_val = 1.5 + 4.0 / Reb + (Hk - Ho) * (Hk - Ho) * bb;
            dHs = -4.0 / (Reb * Reb) * dReb +
                  2.0 * (Hk - Ho) * (dHk - dHo) * bb +
                  (Hk - Ho) * (Hk - Ho) * dbb;
        }

        // Mach correction: Hs = (Hs + 0.028*M2) / (1 + 0.014*M2)
        auto [M2, M2_U2] = GetMach2(U, param);
        double dM2 = GetMach2_dMa(U, param, pd);
        double den = 1.0 + 0.014 * M2;
        Hs_val = (Hs_val + 0.028 * M2) / den;
        dHs = (dHs + 0.028 * dM2) / den - Hs_val / den * 0.014 * dM2;
    } else {
        // Laminar: Hs depends only on Hk
        double a = Hk - 4.35;
        double Hs_Hk;

        if (Hk < 4.35) {
            double num = 0.0111 * a * a - 0.0278 * a * a * a;
            Hs_Hk = (0.0111 * 2.0 * a - 0.0278 * 3.0 * a * a) / (Hk + 1.0) -
                    num / ((Hk + 1.0) * (Hk + 1.0)) -
                    0.0002 * 2.0 * a * Hk * (Hk + a);
        } else {
            Hs_Hk = 0.015 * 2.0 * a / Hk - 0.015 * a * a / (Hk * Hk);
        }
        dHs = Hs_Hk * dHk;
    }

    return dHs;
}

static double GetHss_dMa(const Vector4d& U, const Param& param,
                         const ParamDeriv& pd) {
    // Hss = M2 * (0.064/(Hk-0.8) + 0.251)
    auto [M2, M2_U] = GetMach2(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);
    double dM2 = GetMach2_dMa(U, param, pd);
    double dHk = GetHk_dMa(U, param, pd);

    double num = 0.064 / (Hk - 0.8) + 0.251;
    double dnum = -0.064 / ((Hk - 0.8) * (Hk - 0.8)) * dHk;

    return dM2 * num + M2 * dnum;
}

static double GetDe_dMa(const Vector4d& U, const Param& param,
                        const ParamDeriv& pd) {
    // de depends on Hk through aa = 3.15 + 1.72/(Hk-1)
    auto [Hk, Hk_U] = GetHk(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double th = std::max(U(0), 1e-12);
    double ds = U(1);

    double aa = 3.15 + 1.72 / (Hk - 1.0);
    double daa = -1.72 / ((Hk - 1.0) * (Hk - 1.0)) * dHk;

    double de = th * aa + ds;
    double dde = th * daa;

    // Cap check
    if (de > 12.0 * th) {
        dde = 0.0;  // de = 12*th, no param dependence
    }

    return dde;
}

static double GetDamp_dMa(const Vector4d& U, const Param& param,
                          const ParamDeriv& pd) {
    // Amplification rate dn/dx
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);
    double th = U(0);

    // Limit Hk
    if (Hk < 1.05) { Hk = 1.05; dHk = 0.0; }

    double Hmi = 1.0 / (Hk - 1.0);
    double dHmi = -Hmi * Hmi * dHk;

    double aa = 2.492 * std::pow(Hmi, 0.43);
    double daa = 0.43 * aa / Hmi * dHmi;

    double bb = std::tanh(14.0 * Hmi - 9.24);
    double dbb = (1.0 - bb * bb) * 14.0 * dHmi;

    double lrc = aa + 0.7 * (bb + 1.0);
    double dlrc = daa + 0.7 * dbb;

    constexpr double lten = 2.302585092994046;
    double lr = std::log(Ret) / lten;
    double dlr = (1.0 / Ret) * dRet / lten;

    double ddamp = 0.0;

    constexpr double dl = 0.1;
    if (lr >= lrc - dl) {
        double rn = (lr - (lrc - dl)) / (2.0 * dl);
        double drn = (dlr - dlrc) / (2.0 * dl);

        double rf, drf_coeff;
        if (rn >= 1.0) {
            rf = 1.0;
            drf_coeff = 0.0;
        } else {
            rf = 3.0 * rn * rn - 2.0 * rn * rn * rn;
            drf_coeff = 6.0 * rn - 6.0 * rn * rn;
        }
        double drf = drf_coeff * drn;

        double ar = 3.87 * Hmi - 2.52;
        double dar = 3.87 * dHmi;

        double ex = std::exp(-ar * ar);
        double dex = ex * (-2.0 * ar * dar);

        double da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
        double dda = 0.028 * dHk - 0.0345 * dex;

        double af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi +
                    3.0 * Hmi * Hmi * Hmi + 0.1 * std::exp(-20.0 * Hmi);
        double daf =
            (2.7 - 11.0 * Hmi + 9.0 * Hmi * Hmi - 2.0 * std::exp(-20.0 * Hmi)) *
            dHmi;

        ddamp = (drf * af * da + rf * daf * da + rf * af * dda) / th;
        // damp_U(0) -= damp / th is for dth, not for dMa
    }

    // Extra amplification: depends on param.ncrit which doesn't change with Ma
    // ed = eex * 0.001 / th where eex depends only on U(2) and ncrit
    // No Ma contribution from ed.

    return ddamp;
}

// ============================================================================
// Level 3: GetUs_dMa, GetCf_dMa, GetCDi_Lam_dMa, GetCteq_dMa
// ============================================================================
static double GetUs_dMa(const Vector4d& U, const Param& param,
                        const ParamDeriv& pd) {
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [H, H_U] = GetH(U);
    double dHs = GetHs_dMa(U, param, pd);
    double dHk = GetHk_dMa(U, param, pd);
    // dH = 0 (H = ds/th, no param dependence)

    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) { Hk = Hk_lim; dHk = 0.0; }

    double beta = param.GB;
    double bi = 1.0 / beta;

    double Us = 0.5 * Hs * (1.0 - bi * (Hk - 1.0) / H);
    double dUs = 0.5 * dHs * (1.0 - bi * (Hk - 1.0) / H) +
                 0.5 * Hs * (-bi * dHk / H);
    // dH = 0, so the H term in Us_U drops out

    if (!param.wake && Us > 0.98) {
        dUs = 0.0;
    }

    return dUs;
}

static double GetCf_dMa(const Vector4d& U, const Param& param,
                        const ParamDeriv& pd) {
    if (param.wake) return 0.0;

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);

    double dcf;

    if (param.turb) {
        auto [M2, M2_U] = GetMach2(U, param);
        double dM2 = GetMach2_dMa(U, param, pd);

        double Fc = std::sqrt(1.0 + 0.5 * (param.gam - 1.0) * M2);
        double dFc = 0.5 / Fc * 0.5 * (param.gam - 1.0) * dM2;

        double aa = -1.33 * Hk;
        double daa = -1.33 * dHk;
        if (aa < -17.0) {
            double old_aa = aa;
            aa = -20.0 + 3.0 * std::exp((old_aa + 17.0) / 3.0);
            daa = (aa + 20.0) / 3.0 * daa;
        }

        double bb = std::log(Ret / Fc);
        double dbb = dRet / Ret - dFc / Fc;
        if (bb < 3.0) { bb = 3.0; dbb = 0.0; }
        bb /= std::log(10.0);
        dbb /= std::log(10.0);

        double cc = -1.74 - 0.31 * Hk;
        double dcc = -0.31 * dHk;

        double dd = std::tanh(4.0 - Hk / 0.875);
        double ddd = (1.0 - dd * dd) * (-dHk / 0.875);

        double cf0 = 0.3 * std::exp(aa) * std::pow(bb, cc);
        double dcf0 = cf0 * daa +
                      0.3 * std::exp(aa) * cc * std::pow(bb, cc - 1.0) * dbb +
                      cf0 * std::log(bb) * dcc;

        double cf = (cf0 + 1.1e-4 * (dd - 1.0)) / Fc;
        dcf = (dcf0 + 1.1e-4 * ddd) / Fc - cf / Fc * dFc;
    } else {
        // Laminar cf
        double num, num_Hk;
        if (Hk < 5.5) {
            double term = std::pow(5.5 - Hk, 3) / (Hk + 1.0);
            num = 0.0727 * term - 0.07;
            num_Hk =
                0.0727 * (3.0 * std::pow(5.5 - Hk, 2) / (Hk + 1.0) * (-1.0) -
                          std::pow(5.5 - Hk, 3) / std::pow(Hk + 1.0, 2));
        } else {
            double invterm = 1.0 - 1.0 / (Hk - 4.5);
            num = 0.015 * invterm * invterm - 0.07;
            num_Hk = 0.015 * 2.0 * invterm / std::pow(Hk - 4.5, 2);
        }
        dcf = num_Hk * dHk / Ret - num / (Ret * Ret) * dRet;
    }

    return dcf;
}

static double GetCDi_Lam_dMa(const Vector4d& U, const Param& param,
                              const ParamDeriv& pd) {
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);

    double num, num_Hk;
    if (Hk < 4.0) {
        num = 0.00205 * std::pow(4.0 - Hk, 5.5) + 0.207;
        num_Hk = 0.00205 * 5.5 * std::pow(4.0 - Hk, 4.5) * (-1.0);
    } else {
        double Hk1 = Hk - 4.0;
        double denom = 1.0 + 0.02 * Hk1 * Hk1;
        num = -0.0016 * Hk1 * Hk1 / denom + 0.207;
        num_Hk = -0.0016 * (2.0 * Hk1 / denom -
                            Hk1 * Hk1 / (denom * denom) * 0.02 * 2.0 * Hk1);
    }

    return num_Hk * dHk / Ret - num / (Ret * Ret) * dRet;
}

static double GetCteq_dMa(const Vector4d& U, const Param& param,
                          const ParamDeriv& pd) {
    double CC = 0.5 / (param.GA * param.GA * param.GB);
    double C = param.GC;

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [H, H_U] = GetH(U);
    auto [Ret, Ret_U] = GetRet(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double dHs = GetHs_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);
    double dUs = GetUs_dMa(U, param, pd);
    // dH = 0

    double Hkc;
    double dHkc;
    if (param.wake) {
        if (Hk < 1.00005) { Hk = 1.00005; dHk = 0.0; }
        Hkc = Hk - 1.0;
        dHkc = dHk;
    } else {
        if (Hk < 1.05) { Hk = 1.05; dHk = 0.0; }
        Hkc = Hk - 1.0 - C / Ret;
        dHkc = dHk + C / (Ret * Ret) * dRet;
        if (Hkc < 0.01) { Hkc = 0.01; dHkc = 0.0; }
    }

    double num = CC * Hs * (Hk - 1.0) * Hkc * Hkc;
    double dnum = CC * (dHs * (Hk - 1.0) * Hkc * Hkc + Hs * dHk * Hkc * Hkc +
                        Hs * (Hk - 1.0) * 2.0 * Hkc * dHkc);

    double den = (1.0 - Us) * H * Hk * Hk;
    double dden = (-dUs) * H * Hk * Hk + (1.0 - Us) * H * 2.0 * Hk * dHk;
    // dH = 0

    double cteq = std::sqrt(num / den);
    return 0.5 / cteq * (dnum / den - num / (den * den) * dden);
}

// ============================================================================
// Level 4: GetCDi_TurbWall_dMa, GetCDi_Outer_dMa, GetCDi_LamStress_dMa,
//          GetCDi_LamWake_dMa
// ============================================================================
static double GetCDi_TurbWall_dMa(const Vector4d& U, const Param& param,
                                   const ParamDeriv& pd) {
    if (param.wake) return 0.0;

    auto [cf, cf_U] = GetCf(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dcf = GetCf_dMa(U, param, pd);
    double dHk = GetHk_dMa(U, param, pd);
    double dHs = GetHs_dMa(U, param, pd);
    double dUs = GetUs_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);

    double lr = std::log(Ret);
    double dlr = dRet / Ret;

    double Hmin = 1.0 + 2.1 / lr;
    double dHmin = -2.1 / (lr * lr) * dlr;

    double aa = std::tanh((Hk - 1.0) / (Hmin - 1.0));
    double fac = 0.5 + 0.5 * aa;
    double dfac = 0.5 * (1.0 - aa * aa) *
                  (dHk / (Hmin - 1.0) -
                   (Hk - 1.0) / ((Hmin - 1.0) * (Hmin - 1.0)) * dHmin);

    double cDi = 0.5 * cf * Us * (2.0 / Hs) * fac;
    return dcf * Us / Hs * fac + cf * dUs / Hs * fac -
           cDi / Hs * dHs + cf * Us / Hs * dfac;
}

static double GetCDi_Outer_dMa(const Vector4d& U, const Param& param,
                                const ParamDeriv& pd) {
    if (!param.turb) return 0.0;

    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    double dHs = GetHs_dMa(U, param, pd);
    double dUs = GetUs_dMa(U, param, pd);

    double ct = U(2) * U(2);  // sa^2, no param dependence

    return ct * (-dUs) * 2.0 / Hs - ct * (0.995 - Us) * 2.0 / (Hs * Hs) * dHs;
}

static double GetCDi_LamStress_dMa(const Vector4d& U, const Param& param,
                                    const ParamDeriv& pd) {
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dHs = GetHs_dMa(U, param, pd);
    double dUs = GetUs_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);

    double HsRet = Hs * Ret;
    double dHsRet = dHs * Ret + Hs * dRet;

    double num = 0.15 * (0.995 - Us) * (0.995 - Us) * 2.0;
    double dnum = 0.15 * 2.0 * (0.995 - Us) * (-1.0) * 2.0 * dUs;

    return dnum / HsRet - num / (HsRet * HsRet) * dHsRet;
}

static double GetCDi_LamWake_dMa(const Vector4d& U, const Param& paramin,
                                  const ParamDeriv& pd) {
    // Force laminar mode
    Param param = paramin;
    param.turb = false;

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double dHk = GetHk_dMa(U, param, pd);
    double dHs = GetHs_dMa(U, param, pd);
    double dRet = GetRet_dMa(U, param, pd);

    double HsRet = Hs * Ret;
    double dHsRet = dHs * Ret + Hs * dRet;

    if (Hk < 1.05) { Hk = 1.05; dHk = 0.0; }

    double Hk_inv = 1.0 / Hk;
    double one_m = 1.0 - Hk_inv;
    double num = 2.0 * 1.1 * one_m * one_m * Hk_inv;
    double num_Hk = 2.0 * 1.1 *
                    (2.0 * one_m * (Hk_inv * Hk_inv) * Hk_inv +
                     one_m * one_m * (-Hk_inv * Hk_inv));

    return num_Hk * dHk / HsRet - num / (HsRet * HsRet) * dHsRet;
}

// ============================================================================
// Level 5: GetCDi_dMa, GetUpw_dMa, GetCfxt_dMa, GetCDixt_dMa
// ============================================================================
static double GetCDi_dMa(const Vector4d& U, const Param& param,
                          const ParamDeriv& pd) {
    if (param.turb) {
        double dcDi_turb = 0.0;

        if (!param.wake) {
            dcDi_turb += GetCDi_TurbWall_dMa(U, param, pd);
        }
        dcDi_turb += GetCDi_Outer_dMa(U, param, pd);
        dcDi_turb += GetCDi_LamStress_dMa(U, param, pd);

        // Check against laminar floor using the public GetCDi
        // GetCDi internally computes turbulent sum and compares with laminar floor.
        // To determine which branch wins, compare laminar value with turbulent sum.
        // Reconstruct turbulent sum value:
        // We can't call the static helpers, but we can evaluate GetCDi and check
        // if the laminar branch derivative should be used.

        // Compute laminar floor value
        Param param_lam = param;
        double cDi_lam;
        if (param.wake) {
            // LamWake: force laminar mode for sub-functions
            param_lam.turb = false;
            auto [Hk, Hk_U] = GetHk(U, param_lam);
            auto [Hs, Hs_U] = GetHs(U, param_lam);
            auto [Ret, Ret_U] = GetRet(U, param_lam);
            if (Hk < 1.05) Hk = 1.05;
            double Hk_inv = 1.0 / Hk;
            double one_m = 1.0 - Hk_inv;
            double num = 2.0 * 1.1 * one_m * one_m * Hk_inv;
            cDi_lam = num / (Hs * Ret);
        } else {
            param_lam.turb = false;
            auto [Hk, Hk_U] = GetHk(U, param_lam);
            auto [Ret, Ret_U] = GetRet(U, param_lam);
            double num;
            if (Hk < 4.0) {
                num = 0.00205 * std::pow(4.0 - Hk, 5.5) + 0.207;
            } else {
                double Hk1 = Hk - 4.0;
                double denom = 1.0 + 0.02 * Hk1 * Hk1;
                num = -0.0016 * Hk1 * Hk1 / denom + 0.207;
            }
            cDi_lam = num / Ret;
        }

        // Compute turbulent sum (without laminar floor)
        // The turbulent value is GetCDi - cDi_lam when cDi_lam is the floor
        // Actually, we need the sum before the max check.
        // Use GetCDi to get the final value and compare with cDi_lam.
        auto [cDi_total, cDi_total_U] = GetCDi(U, param);
        // If wake, cDi_total includes the 2x factor, and cDi_lam doesn't.
        double cDi_lam_adj = param.wake ? 2.0 * cDi_lam : cDi_lam;

        // If the laminar floor dominates, use laminar derivative
        // The total value equals cDi_lam_adj when floor wins
        double dcDi;
        if (std::abs(cDi_total - cDi_lam_adj) < 1e-15 * std::max(1.0, std::abs(cDi_total))) {
            // Laminar floor wins
            dcDi = param.wake ? GetCDi_LamWake_dMa(U, param, pd)
                              : GetCDi_Lam_dMa(U, param, pd);
            if (param.wake) dcDi *= 2.0;
        } else {
            dcDi = dcDi_turb;
            if (param.wake) dcDi *= 2.0;
        }

        return dcDi;
    } else {
        return GetCDi_Lam_dMa(U, param, pd);
    }
}

// GetUpw_dMa: upwinding depends only on Hk which depends on Ma through M2
static double GetUpw_dMa(const Vector4d& U1, const Vector4d& U2,
                          const Param& param, const ParamDeriv& pd) {
    auto [Hk1, Hk1_U] = GetHk(U1, param);
    auto [Hk2, Hk2_U] = GetHk(U2, param);
    double dHk1 = GetHk_dMa(U1, param, pd);
    double dHk2 = GetHk_dMa(U2, param, pd);

    double C = param.wake ? 1.0 : 5.0;
    double Huc = C * 1.0 / (Hk2 * Hk2);
    double dHuc = -2.0 * Huc / Hk2 * dHk2;

    double aa = (Hk2 - 1.0) / (Hk1 - 1.0);
    double sga = (aa >= 0) ? 1.0 : -1.0;
    double la = std::log(sga * aa);
    double dla = -dHk1 / (Hk1 - 1.0) + dHk2 / (Hk2 - 1.0);

    double Hls = la * la;
    double dHls = 2.0 * la * dla;
    if (Hls > 15.0) { Hls = 15.0; dHls = 0.0; }

    return -0.5 * std::exp(-Hls * Huc) * (-dHls * Huc - Hls * dHuc);
}

static double GetCfxt_dMa(const Vector4d& U, double x, const Param& param,
                           const ParamDeriv& pd) {
    double dcf = GetCf_dMa(U, param, pd);
    double th = std::max(U(0), 1e-12);
    return dcf * x / th;
}

static double GetCDixt_dMa(const Vector4d& U, double x, const Param& param,
                            const ParamDeriv& pd) {
    double dcDi = GetCDi_dMa(U, param, pd);
    double th = std::max(U(0), 1e-12);
    return dcDi * x / th;
}

// ============================================================================
// Level 6: GetUq_dMa
// ============================================================================
static double GetUq_dMa(double ds, double /*cf*/, double dcf_dMa,
                        double Hk, double dHk_dMa,
                        double Ret, double dRet_dMa,
                        const Param& param) {
    double beta = param.GB;
    double A = param.GA;
    double C = param.GC;

    if (param.wake) { A *= param.Dlr; C = 0.0; }

    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) { Hk = Hk_lim; dHk_dMa = 0.0; }

    double Hkc = Hk - 1.0 - C / std::max(Ret, 1e-10);
    double dHkc = dHk_dMa + C / (Ret * Ret) * dRet_dMa;
    if (Hkc < 0.01) { Hkc = 0.01; dHkc = 0.0; }

    double dut = 0.5 * dcf_dMa -
                 2.0 * (Hkc / (A * Hk)) *
                     (dHkc / (A * Hk) - Hkc / (A * Hk * Hk) * dHk_dMa);

    double ds_safe = std::max(ds, 1e-12);
    // uq = ut / (beta * ds_safe), ds has no param dependence here (it's averaged from U)
    return dut / (beta * ds_safe);
}

// ============================================================================
// Top level: ResidualStation_dMa — computes 3-component dR/dMa at one station
// ============================================================================
static Vector3d ResidualStation_dMa(const Param& param, const Vector2d& x,
                                     const Matrix42d& Uin, const Vector2d& Aux,
                                     const ParamDeriv& pd) {
    // Copy and modify for wake gap (same as primal)
    Matrix42d U = Uin;
    U.row(1) -= Aux.transpose();

    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);

    auto [uk1, uk1_U] = GetUk(U1(3), param);
    auto [uk2, uk2_U] = GetUk(U2(3), param);
    // GetUk is a stub: uk = ue, no Ma dependence => duk_dMa = 0
    // So uelog_dMa = 0

    // Similarity handling
    bool simi = param.simi;

    double dx = simi ? 0.5 * (x(0) + x(1)) : (x(1) - x(0));

    // Upwinding
    auto [upw, upw_U] = GetUpw(U1, U2, param);
    double dupw = GetUpw_dMa(U1, U2, param, pd);

    // --- Momentum equation dRmom/dMa ---
    // Rmom = thlog + (2+H+Hw-Ms)*uelog - 0.5*xlog*cfxt
    // thlog, uelog, xlog have no Ma dependence
    // H has no Ma dependence (H = ds/th)
    // Hw has no Ma dependence
    // Ms (Mach^2) depends on Ma
    // cfxt depends on Ma
    double dMs1 = GetMach2_dMa(U1, param, pd);
    double dMs2 = GetMach2_dMa(U2, param, pd);
    double dMs = 0.5 * (dMs1 + dMs2);  // simple average

    double uelog;
    if (simi) {
        uelog = 1.0;
    } else {
        uelog = std::log(std::max(uk2, 1e-12) / std::max(uk1, 1e-12));
    }

    double xlog;
    if (simi) {
        xlog = 1.0;
    } else {
        xlog = std::log(std::max(x(1), 1e-12) / std::max(x(0), 1e-12));
    }

    // cfxt = 0.25*cfxt1 + 0.5*cfxtm + 0.25*cfxt2
    double dcfxt1 = GetCfxt_dMa(U1, x(0), param, pd);
    double xm = 0.5 * (x(0) + x(1));
    Vector4d Um = 0.5 * (U1 + U2);
    double dcfxtm = GetCfxt_dMa(Um, xm, param, pd);
    double dcfxt2 = GetCfxt_dMa(U2, x(1), param, pd);
    double dcfxt = 0.25 * dcfxt1 + 0.5 * dcfxtm + 0.25 * dcfxt2;

    double dRmom = -dMs * uelog - 0.5 * xlog * dcfxt;

    // --- Shape parameter equation dRshape/dMa ---
    // Rshape = Hslog + (2*Hss/Hs + 1 - H - Hw)*uelog + xlog*(0.5*cfxtu - cDixt)
    auto [Hs1, Hs1_U] = GetHs(U1, param);
    auto [Hs2, Hs2_U] = GetHs(U2, param);
    double dHs1 = GetHs_dMa(U1, param, pd);
    double dHs2 = GetHs_dMa(U2, param, pd);

    // Hslog = log(Hs2/Hs1)
    double dHslog = simi ? 0.0 : (dHs2 / std::max(Hs2, 1e-12) - dHs1 / std::max(Hs1, 1e-12));

    double Hs = 0.5 * (Hs1 + Hs2);
    double dHs = 0.5 * (dHs1 + dHs2);

    auto [Hss1, Hss1_U] = GetHss(U1, param);
    auto [Hss2, Hss2_U] = GetHss(U2, param);
    double dHss1 = GetHss_dMa(U1, param, pd);
    double dHss2 = GetHss_dMa(U2, param, pd);
    double Hss = 0.5 * (Hss1 + Hss2);
    double dHss = 0.5 * (dHss1 + dHss2);

    // Upwinded cDixt
    auto [cDixt1, cDixt1_U, cDixt1_x] = GetCDixt(U1, x(0), param);
    auto [cDixt2, cDixt2_U, cDixt2_x] = GetCDixt(U2, x(1), param);
    double dcDixt1 = GetCDixt_dMa(U1, x(0), param, pd);
    double dcDixt2 = GetCDixt_dMa(U2, x(1), param, pd);
    double dcDixt = (1.0 - upw) * dcDixt1 + upw * dcDixt2 +
                    dupw * (cDixt2 - cDixt1);

    // Upwinded cfxtu
    auto [cfxt1_val, cfxt1_val_U, cfxt1_val_x] = GetCfxt(U1, x(0), param);
    auto [cfxt2_val, cfxt2_val_U, cfxt2_val_x] = GetCfxt(U2, x(1), param);
    double dcfxtu1 = GetCfxt_dMa(U1, x(0), param, pd);
    double dcfxtu2 = GetCfxt_dMa(U2, x(1), param, pd);
    double dcfxtu = (1.0 - upw) * dcfxtu1 + upw * dcfxtu2 +
                    dupw * (cfxt2_val - cfxt1_val);

    double dRshape = dHslog +
                     (2.0 * dHss / Hs - 2.0 * Hss / (Hs * Hs) * dHs) * uelog +
                     xlog * (0.5 * dcfxtu - dcDixt);

    // --- Lag equation dRlag/dMa ---
    double dRlag;

    if (param.turb) {
        auto [de1, de1_U] = GetDe(U1, param);
        auto [de2, de2_U] = GetDe(U2, param);
        double dde1 = GetDe_dMa(U1, param, pd);
        double dde2 = GetDe_dMa(U2, param, pd);
        double de = 0.5 * (de1 + de2);
        double dde = 0.5 * (dde1 + dde2);

        auto [Us1, Us1_U] = GetUs(U1, param);
        auto [Us2, Us2_U] = GetUs(U2, param);
        double dUs1 = GetUs_dMa(U1, param, pd);
        double dUs2 = GetUs_dMa(U2, param, pd);
        double Us_upw = 0.5 * (Us1 + Us2);
        double dUs_upw = 0.5 * (dUs1 + dUs2);

        auto [Hk1, Hk1_U] = GetHk(U1, param);
        auto [Hk2, Hk2_U] = GetHk(U2, param);
        double dHk1 = GetHk_dMa(U1, param, pd);
        double dHk2 = GetHk_dMa(U2, param, pd);
        auto [Hk_upw_bind, Hk_upw_U] = Upwind(upw, upw_U, Hk1, Hk1_U, Hk2, Hk2_U);
        double Hk_upw = Hk_upw_bind;
        double dHk_upw = (1.0 - upw) * dHk1 + upw * dHk2 + dupw * (Hk2 - Hk1);

        auto [Ret1, Ret1_U] = GetRet(U1, param);
        auto [Ret2, Ret2_U] = GetRet(U2, param);
        double dRet1 = GetRet_dMa(U1, param, pd);
        double dRet2 = GetRet_dMa(U2, param, pd);
        double Ret_avg = 0.5 * (Ret1 + Ret2);
        double dRet_avg = 0.5 * (dRet1 + dRet2);

        auto [cf1, cf1_U] = GetCf(U1, param);
        auto [cf2, cf2_U] = GetCf(U2, param);
        double dcf1 = GetCf_dMa(U1, param, pd);
        double dcf2 = GetCf_dMa(U2, param, pd);
        double cf_upw = (1.0 - upw) * cf1 + upw * cf2;
        double dcf_upw = (1.0 - upw) * dcf1 + upw * dcf2 + dupw * (cf2 - cf1);

        // ds averaged (no param dep)
        double dsa_avg = 0.5 * (U1(1) + U2(1));

        // uq
        double uq_val = [&]() {
            auto [uq, uq_U] = GetUq(dsa_avg, Vector8d::Zero(), cf_upw, Vector8d::Zero(),
                                     Hk_upw, Vector8d::Zero(), Ret_avg, Vector8d::Zero(), param);
            return uq;
        }();
        double duq = GetUq_dMa(dsa_avg, cf_upw, dcf_upw, Hk_upw, dHk_upw,
                                Ret_avg, dRet_avg, param);

        auto [cteq1, cteq1_U] = GetCteq(U1, param);
        auto [cteq2, cteq2_U] = GetCteq(U2, param);
        double dcteq1 = GetCteq_dMa(U1, param, pd);
        double dcteq2 = GetCteq_dMa(U2, param, pd);
        double cteq_upw = (1.0 - upw) * cteq1 + upw * cteq2;
        double dcteq_upw = (1.0 - upw) * dcteq1 + upw * dcteq2 +
                           dupw * (cteq2 - cteq1);

        // saa upwinded (sa values, no param dep => dsaa = dupw*(sa2-sa1))
        double sa1 = std::max(U1(2), 1e-12);
        double sa2 = std::max(U2(2), 1e-12);
        double saa = (1.0 - upw) * sa1 + upw * sa2;
        double dsaa = dupw * (sa2 - sa1);

        double Klag = param.SlagK;
        double beta = param.GB;
        double Clag = Klag / beta / (1.0 + Us_upw);
        double dClag = -Clag / (1.0 + Us_upw) * dUs_upw;

        double ald = param.wake ? param.Dlr : 1.0;
        double Cuq = param.Cuq;

        // salog has no Ma dependency
        double salog = std::log(std::max(sa2, 1e-10) / std::max(sa1, 1e-10));

        // Rlag = Clag*(cteq - ald*saa)*dx - 2*de*salog + 2*de*(uq*dx - uelog)*Cuq
        dRlag = dClag * (cteq_upw - ald * saa) * dx +
                Clag * (dcteq_upw - ald * dsaa) * dx -
                2.0 * dde * salog +
                2.0 * dde * (uq_val * dx - uelog) * Cuq +
                2.0 * de * (duq * dx) * Cuq;
        // uelog_dMa = 0 since GetUk is stub
    } else {
        if (simi) {
            dRlag = 0.0;
        } else {
            double ddamp1 = GetDamp_dMa(U1, param, pd);
            double ddamp2 = GetDamp_dMa(U2, param, pd);
            double ddamp = 0.5 * (ddamp1 + ddamp2);
            // Rlag = sa2 - sa1 - damp*dx => dRlag/dMa = -ddamp*dx
            dRlag = -ddamp * dx;
        }
    }

    return Vector3d(dRmom, dRshape, dRlag);
}

// Helper: d(cttr)/dMa (transition initial ctau)
static double GetCttr_dMa(const Vector4d& Ut, const Param& paramin,
                           const ParamDeriv& pd) {
    // cttr = CtauC * exp(-CtauE/(Hk-1)) * cteq
    Param p = paramin;
    p.wake = false;

    auto [Hk, Hk_U] = GetHk(Ut, p);
    double dHk = GetHk_dMa(Ut, p, pd);
    auto [cteq, cteq_U] = GetCteq(Ut, p);
    double dcteq = GetCteq_dMa(Ut, p, pd);

    if (Hk < 1.05) { Hk = 1.05; dHk = 0.0; }

    double E = p.CtauE;
    double c = p.CtauC * std::exp(-E / (Hk - 1.0));
    double dc = c * E / ((Hk - 1.0) * (Hk - 1.0)) * dHk;

    return dc * cteq + c * dcteq;
}

// ============================================================================
// ResidualTransition_dMa: fully analytical dR/dMa at a transition station
//
// R = R_lam(param_lam, x1, xt, U1, Utl) + R_turb(param_turb, xt, x2, Utt, U2)
// where xt satisfies Rxt(xt, U, x, param) = ncrit - sa1 - 0.5*(xt-x1)*(damp1+dampt) = 0
//
// dR/dMa = ∂R_lam/∂Ma + ∂R_turb/∂Ma
//        + (Rl_x(:,1) + Rt_x(:,0)) * dxt/dMa
//        + Rl_Utl * dUtl/dMa + Rt_Utt * dUtt/dMa
// ============================================================================
static Vector3d ResidualTransition_dMa(FoilState& /*state*/, const Param& param_base,
                                        const Vector2d& x, const Matrix42d& U,
                                        const Vector2d& Aux,
                                        const ParamDeriv& pd) {
    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);
    double x1 = x(0), x2 = x(1);
    double dx = x2 - x1;
    double ncrit = param_base.ncrit;

    // Reconstruct converged transition location (same Newton as primal)
    double xt = x1 + 0.5 * dx;
    for (int iN = 0; iN < 20; ++iN) {
        double w2t = (xt - x1) / dx, w1t = 1.0 - w2t;
        Vector4d Utt = w1t * U1 + w2t * U2;
        Utt(2) = ncrit;
        Vector4d Ut_xt_t = (U2 - U1) / dx;
        Ut_xt_t(2) = 0.0;

        Param pm = param_base; pm.turb = false;
        auto [d1, d1u] = GetDamp(U1, pm);
        auto [dt, dtu] = GetDamp(Utt, pm);
        double Rxt = ncrit - U1(2) - 0.5 * (xt - x1) * (d1 + dt);
        double Rxt_xt = -0.5 * (d1 + dt) - 0.5 * (xt - x1) * dtu.dot(Ut_xt_t);
        double step = -Rxt / Rxt_xt;
        double dmax = 0.2 * dx * (1.1 - static_cast<double>(iN) / 20);
        if (std::abs(step) > dmax) step = step * dmax / std::abs(step);
        if (std::abs(Rxt) < 1e-10) break;
        if (iN < 19) xt += step;
    }

    // At converged xt, set up interpolated state
    double w2 = (xt - x1) / dx, w1 = 1.0 - w2;
    Vector4d Ut = w1 * U1 + w2 * U2;
    Vector4d Ut_xt = (U2 - U1) / dx;
    Ut(2) = ncrit;
    Ut_xt(2) = 0.0;

    Param param_lam = param_base; param_lam.turb = false;
    auto [damp1, damp1_U1] = GetDamp(U1, param_lam);
    auto [dampt, dampt_Ut] = GetDamp(Ut, param_lam);
    dampt_Ut(2) = 0.0;

    double Rxt_xt = -0.5 * (damp1 + dampt) -
                    0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);

    // dxt/dMa from implicit function theorem on Rxt = 0
    // ∂Rxt/∂Ma|_{xt fixed} = -0.5*(xt-x1)*(∂damp1/∂Ma + ∂dampt/∂Ma)
    double ddamp1_dMa = GetDamp_dMa(U1, param_lam, pd);
    double ddampt_dMa = GetDamp_dMa(Ut, param_lam, pd);
    double dRxt_dMa = -0.5 * (xt - x1) * (ddamp1_dMa + ddampt_dMa);
    double dxt_dMa = -dRxt_dMa / Rxt_xt;

    // dUtl/dMa = Ut_xt * dxt_dMa (with sa = 0, since sa = ncrit always)
    Vector4d dUtl_dMa = Ut_xt * dxt_dMa;
    dUtl_dMa(2) = 0.0;

    // dUtt/dMa: position shift + cttr param change
    Param param_turb = param_base; param_turb.turb = true;
    auto [cttr, cttr_Ut] = GetCttr(Ut, param_turb);
    double dcttr_direct = GetCttr_dMa(Ut, param_turb, pd);
    Vector4d dUtt_dMa = Ut_xt * dxt_dMa;
    dUtt_dMa(2) = cttr_Ut.dot(Ut_xt * dxt_dMa) + dcttr_direct;

    // Build transition states
    Vector4d Utl = Ut; Utl(2) = ncrit;
    Vector4d Utt = Ut; Utt(2) = cttr;

    // Laminar half residual and direct Ma derivative
    param_lam.turb = false;
    Matrix42d U_lam; U_lam.col(0) = U1; U_lam.col(1) = Utl;
    Vector2d x_lam(x1, xt);
    auto Rl = ResidualStation(param_lam, x_lam, U_lam, Aux);
    Vector3d dRl_direct = ResidualStation_dMa(param_lam, x_lam, U_lam, Aux, pd);

    // Turbulent half residual and direct Ma derivative
    param_turb.turb = true;
    Matrix42d U_turb; U_turb.col(0) = Utt; U_turb.col(1) = U2;
    Vector2d x_turb(xt, x2);
    auto Rt = ResidualStation(param_turb, x_turb, U_turb, Aux);
    Vector3d dRt_direct = ResidualStation_dMa(param_turb, x_turb, U_turb, Aux, pd);

    // Chain rule assembly
    Matrix34d Rl_Utl = Rl.lin_U.block<3, 4>(0, 4);
    Matrix34d Rt_Utt = Rt.lin_U.block<3, 4>(0, 0);

    return dRl_direct + dRt_direct +
           (Rl.lin_x.col(1) + Rt.lin_x.col(0)) * dxt_dMa +
           Rl_Utl * dUtl_dMa + Rt_Utt * dUtt_dMa;
}

// ============================================================================
// ComputeDRblDMa: assemble dR_BL/dMa vector [3*Nsys]
// Mirrors BuildGlobSys exactly for station ordering and stagnation handling.
// ============================================================================
static VectorXd ComputeDRblDMa(FoilState& state, const ParamDeriv& pd) {
    int Nsys = state.glob.Nsys;
    VectorXd dR(3 * Nsys);
    dR.setZero();

    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        if (Is.empty()) continue;
        int Ns = static_cast<int>(Is.size());

        Param param = BuildParam(state, si);

        // Local xi and Aux, matching BuildGlobSys exactly
        VectorXd xi(Ns);
        for (int idx = 0; idx < Ns; ++idx) xi(idx) = state.isol.xi(Is[idx]);

        VectorXd Aux_surf = VectorXd::Zero(Ns);
        if (si == 2 && state.vsol.wgap.size() > 0) {
            for (int idx = 0;
                 idx < Ns && idx < static_cast<int>(state.vsol.wgap.size());
                 ++idx) {
                Aux_surf(idx) = state.vsol.wgap(idx);
            }
        }

        int i0 = (si < 2 && xi(0) < 1e-8 * xi(Ns - 1)) ? 1 : 0;

        if (si < 2) {
            // Stagnation: similarity residual at extrapolated Ust
            int ip0 = i0, ip1 = i0 + 1;
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip0]);
            Uin.col(1) = state.glob.U.col(Is[ip1]);
            Vector2d xin(xi(ip0), xi(ip1));

            auto stag = StagnationState(Uin, xin);  // No param dependence

            param.turb = false;
            param.simi = true;
            Matrix42d Ust_pair;
            Ust_pair.col(0) = stag.Ust;
            Ust_pair.col(1) = stag.Ust;
            Vector2d xst_pair(stag.xst, stag.xst);
            Vector2d Aux_pair(Aux_surf(ip0), Aux_surf(ip0));

            Vector3d dRs = ResidualStation_dMa(param, xst_pair, Ust_pair,
                                                Aux_pair, pd);
            dR.segment<3>(3 * Is[ip1]) = dRs;
            param.simi = false;

            // Constraint U[Is[i0]] = Ust at first node: no param dependence → dR = 0
        } else {
            // Wake init: WakeSys
            int il = state.vsol.Is[0].back();
            int iu = state.vsol.Is[1].back();
            int iw = Is[0];

            double th_l = state.glob.U(0, il);
            double th_u = state.glob.U(0, iu);
            double wl = th_l / (th_l + th_u);
            double wu = th_u / (th_l + th_u);

            // R1,R2 = theta/ds sums: no param dep
            // R3 = sa_w - weighted ctau: ctl/ctu depend on param if laminar
            double dctl = 0.0;
            if (il < static_cast<int>(state.vsol.turb.size()) &&
                !state.vsol.turb[il]) {
                Param p = param; p.wake = false;
                StationParam(state, p, il); p.wake = false;
                dctl = GetCttr_dMa(state.glob.U.col(il), p, pd);
            }
            double dctu = 0.0;
            if (iu < static_cast<int>(state.vsol.turb.size()) &&
                !state.vsol.turb[iu]) {
                Param p = param; p.wake = false;
                StationParam(state, p, iu); p.wake = false;
                dctu = GetCttr_dMa(state.glob.U.col(iu), p, pd);
            }

            dR(3 * iw + 2) = -(wl * dctl + wu * dctu);
            param.turb = true;
            param.wake = true;
        }

        // March through remaining stations
        int idx_start = (si < 2) ? i0 + 2 : 1;
        for (int idx = idx_start; idx < Ns; ++idx) {
            int ip = idx - 1;
            int i_glob = Is[idx];
            int ip_glob = Is[ip];

            Vector2d xin(xi(ip), xi(idx));
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(ip_glob);
            Uin.col(1) = state.glob.U.col(i_glob);
            Vector2d Aux_in(Aux_surf(ip), Aux_surf(idx));

            bool tran = (ip_glob < static_cast<int>(state.vsol.turb.size()) &&
                         i_glob < static_cast<int>(state.vsol.turb.size())) &&
                        (state.vsol.turb[ip_glob] != state.vsol.turb[i_glob]);

            Vector3d dRs;
            if (tran) {
                dRs = ResidualTransition_dMa(state, param, xin, Uin, Aux_in, pd);
                param.turb = true;
            } else {
                StationParam(state, param, i_glob);
                dRs = ResidualStation_dMa(param, xin, Uin, Aux_in, pd);
            }
            dR.segment<3>(3 * i_glob) = dRs;
        }
    }

    return dR;
}

// ============================================================================
// ComputeDcpDMa: d(cp)/d(Ma) at each panel node (KT correction differentiation)
// For the force partials w.r.t. Ma
// ============================================================================
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
// ============================================================================
void MachSensitivities(FoilState& state, Sensitivity& sens) {
    int N = state.foil.N;
    double Ma = state.oper.Ma;
    double chord = state.geom.chord;
    double alpha = state.oper.alpha;

    if (!state.oper.viscous) {
        // --- Inviscid path ---
        // Only the KT correction on cp depends on Ma
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
    } else {
        // --- Viscous path ---
        // dJ/dMa = ∂J/∂Ma_partial − λ^T · ∂R/∂Ma
        // ∂J/∂Ma_partial: from KT correction (same as inviscid but with viscous ue)
        // ∂R/∂Ma: [dR_BL/dMa; dR_ue/dMa] where dR_ue/dMa = 0

        if (!state.glob.lu_factorization) {
            throw std::runtime_error(
                "MachSensitivities: no LU factorization stored");
        }

        int Nsys = state.glob.Nsys;

        CalcForceLin(state);
        CalcCdLin(state);

        // Compute force partials w.r.t. Ma (KT correction)
        VectorXd ue;
        if (state.glob.U.cols() > 0) {
            ue = state.glob.U.row(3).transpose();
        } else {
            ue = GetUeInv(state);
        }

        double dcl_partial = 0.0, dcm_partial = 0.0;
        double cosa = Cosd(alpha);
        double sina = Sind(alpha);
        Vector2d xref = state.geom.xref;

        // Also compute d(cl_ue)/dMa for the cross-term
        // (not needed: cl_ue * due/dMa = 0 since ue equation doesn't depend on Ma)
        // The partial is just the direct cp-Ma dependency at fixed ue

        for (int i0 = 1; i0 <= N; ++i0) {
            int i = (i0 == N) ? 0 : i0;
            int ip = (i0 == N) ? N - 1 : i0 - 1;

            Vector2d x1 = state.foil.x.col(ip);
            Vector2d x2 = state.foil.x.col(i);
            Vector2d dxv = x2 - x1;
            Vector2d dx1 = x1 - xref;
            Vector2d dx2 = x2 - xref;

            double dx_lift = -dxv(0) * cosa - dxv(1) * sina;

            double dcp_ip, d2cp_ip;
            double dcp_i, d2cp_i;
            ComputeCpDerivsMa(ue(ip), state.param, Ma, dcp_ip, d2cp_ip);
            ComputeCpDerivsMa(ue(i), state.param, Ma, dcp_i, d2cp_i);

            double dcpbar = 0.5 * (dcp_ip + dcp_i);
            dcl_partial += dx_lift * dcpbar;

            double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
            double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
            dcm_partial += dcp_ip / 3.0 * dx1nds + dcp_ip / 6.0 * dx2nds +
                           dcp_i / 6.0 * dx1nds + dcp_i / 3.0 * dx2nds;
        }
        dcl_partial /= chord;
        dcm_partial /= (chord * chord);
        double dcd_partial = 0.0;  // Squire-Young cd has no direct Ma partial

        // Build dJ/dV for adjoint solve (same as alpha sensitivities)
        VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
        VectorXd dcm_dV = VectorXd::Zero(4 * Nsys);
        VectorXd dcd_dV = state.post.cd_V;

        // Need to include d(cl_ue)/dMa * ue_i contribution
        // But wait: the ue values in dJ/dV are the standard cl_ue linearizations.
        // The Ma partial from force is separate from the adjoint term.
        // The adjoint term handles: lambda^T * dR/dMa
        // dR/dMa = [dR_BL/dMa; 0] (ue equation doesn't depend on Ma)

        for (int i = 0; i < N; ++i) {
            dcl_dV(4 * i + 3) = state.post.cl_ue(i);
            dcm_dV(4 * i + 3) = state.post.cm_ue(i);
        }

        // However, we also need to account for d(cl_ue)/dMa (the KT correction
        // changes how cl relates to ue). This enters as:
        // dcl/dMa = dcl/dMa|_ue_fixed + (dcl_ue/dMa) * ue + cl_ue * (due/dMa)
        // Since due/dMa = 0 (ue equation doesn't depend on Ma), the second term
        // is zero. The third term: (dcl_ue/dMa) * ue is already part of dcl_partial.
        // Actually dcl_partial already captures all force-Ma partial at fixed ue.

        auto& lu = *state.glob.lu_factorization;
        VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
        VectorXd lambda_cm = lu.transpose().solve(dcm_dV);
        VectorXd lambda_cd = lu.transpose().solve(dcd_dV);

        // Compute dR_BL/dMa [3*Nsys]
        ParamDeriv pd_ma = InitThermoDeriv_Ma(
            state.param, Ma, state.oper.Re, state.oper.rho,
            state.oper.Vinf, chord);
        VectorXd dRbl_dMa = ComputeDRblDMa(state, pd_ma);

        // dR/dMa = [dRbl_dMa; 0_{Nsys}] (block ordered)
        // lambda^T * dR/dMa = lambda_bl^T * dRbl_dMa
        VectorXd lambda_cl_bl = lambda_cl.head(3 * Nsys);
        VectorXd lambda_cm_bl = lambda_cm.head(3 * Nsys);
        VectorXd lambda_cd_bl = lambda_cd.head(3 * Nsys);

        sens.dcl_dMa = dcl_partial - lambda_cl_bl.dot(dRbl_dMa);
        sens.dcm_dMa = dcm_partial - lambda_cm_bl.dot(dRbl_dMa);
        sens.dcd_dMa = dcd_partial - lambda_cd_bl.dot(dRbl_dMa);
    }
}

// ============================================================================
// ReSensitivities: compute dJ/dRe for cl, cd, cm
// ============================================================================
void ReSensitivities(FoilState& state, Sensitivity& sens) {
    // Re enters only through mu0: dmu0/dRe = -mu0/Re
    // For inviscid: no Re dependence
    if (!state.oper.viscous) {
        sens.dcl_dRe = 0.0;
        sens.dcd_dRe = 0.0;
        sens.dcm_dRe = 0.0;
        return;
    }

    if (!state.glob.lu_factorization) {
        throw std::runtime_error(
            "ReSensitivities: no LU factorization stored");
    }

    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    double Re = state.oper.Re;

    CalcForceLin(state);
    CalcCdLin(state);

    // ParamDeriv for Re: dMinf=0, drho0=0, dmu0=-mu0/Re
    ParamDeriv pd_re;
    pd_re.dMinf = 0.0;
    pd_re.drho0 = 0.0;
    pd_re.dmu0 = -state.param.mu0 / Re;

    // dR_BL/dRe uses the same _dMa companion functions with Re-specific ParamDeriv
    VectorXd dRbl_dRe = ComputeDRblDMa(state, pd_re);

    // Force partials: no direct Re dependence in forces
    // (cp from ue has no Re; Squire-Young cd depends on BL state, not Re directly)
    double dcl_partial = 0.0;
    double dcm_partial = 0.0;
    double dcd_partial = 0.0;

    // Adjoint solve
    VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
    VectorXd dcm_dV = VectorXd::Zero(4 * Nsys);
    VectorXd dcd_dV = state.post.cd_V;
    for (int i = 0; i < N; ++i) {
        dcl_dV(4 * i + 3) = state.post.cl_ue(i);
        dcm_dV(4 * i + 3) = state.post.cm_ue(i);
    }

    auto& lu = *state.glob.lu_factorization;
    VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
    VectorXd lambda_cm = lu.transpose().solve(dcm_dV);
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);

    VectorXd lambda_cl_bl = lambda_cl.head(3 * Nsys);
    VectorXd lambda_cm_bl = lambda_cm.head(3 * Nsys);
    VectorXd lambda_cd_bl = lambda_cd.head(3 * Nsys);

    sens.dcl_dRe = dcl_partial - lambda_cl_bl.dot(dRbl_dRe);
    sens.dcm_dRe = dcm_partial - lambda_cm_bl.dot(dRbl_dRe);
    sens.dcd_dRe = dcd_partial - lambda_cd_bl.dot(dRbl_dRe);
}

}  // namespace cppfoil
