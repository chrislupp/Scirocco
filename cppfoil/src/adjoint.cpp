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
// In cppfoil, foil.s (arc-length parametrization) is computed once during
// paneling and is NOT updated when panel coordinates are perturbed. Therefore
// xi = |s - sstag| does not depend directly on foil.x. The geometric
// coordinate sensitivity enters ONLY through the 4th equation (ueinv, ue_m).
// ============================================================================
SparseMatrixd ComputeDxiDxFoil(const FoilState& state) {
    int Nsys = state.glob.Nsys;
    int N = state.foil.N;
    SparseMatrixd dxi_dx(Nsys, 2 * N);
    return dxi_dx;
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
    if (Nw > 0) {
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
                        Vector2d xi_eval = state.foil.x.col(i);
                        Matrix2d Xj;
                        Xj.col(0) = state.foil.x.col(j);
                        Xj.col(1) = state.foil.x.col(j + 1);
                        auto sd = PanelConstSourceStreamDeriv(Xj, xi_eval);
                        dBphi(i) += sd.val_Xj(xj_idx) * phi(j);
                    }
                }

                // xi evaluation point contributions (row k, ALL panels in B)
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
    // First 3*Nsys rows: BL residuals depend on xi which depends on x
    SparseMatrixd dxi_dx = ComputeDxiDxFoil(state);
    // dR_BL/dx = R_x * dxi/dx [3*Nsys x 2*N]

    // Last Nsys rows: 4th equation
    MatrixXd dRue_dx = ComputeDRueDxFoil(state);

    // Compute lambda^T * dR/dx for each objective
    // Split lambda into BL part (first 3*Nsys) and ue part (last Nsys)
    Sensitivity sens;
    sens.dcl_dx = MatrixXd::Zero(2, N);
    sens.dcm_dx = MatrixXd::Zero(2, N);
    sens.dcd_dx = MatrixXd::Zero(2, N);
    sens.dcdpi_dx = MatrixXd::Zero(2, N);

    auto computeVJP = [&](const VectorXd& lambda) -> VectorXd {
        VectorXd lambda_bl = lambda.head(3 * Nsys);
        VectorXd lambda_ue = lambda.tail(Nsys);

        // BL contribution: lambda_bl^T * R_x * dxi/dx
        VectorXd vjp_bl = (lambda_bl.transpose() * state.glob.R_x * dxi_dx)
                              .transpose();

        // 4th equation contribution: lambda_ue^T * dRue/dx
        VectorXd vjp_ue = dRue_dx.transpose() * lambda_ue;

        return vjp_bl + vjp_ue;
    };

    VectorXd vjp_cl = computeVJP(lambda_cl);
    VectorXd vjp_cm = computeVJP(lambda_cm);
    VectorXd vjp_cd = computeVJP(lambda_cd);

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
