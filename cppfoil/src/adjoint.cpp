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

        double dcp_due_ip = -2.0 * ue(ip) / (Vinf * Vinf);
        double dcp_due_i = -2.0 * ue(i) / (Vinf * Vinf);

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

}  // namespace cppfoil
