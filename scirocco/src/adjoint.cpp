/*
    Scirocco: C++ Airfoil Analysis Library

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
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "scirocco/core/adjoint.h"

#include <array>
#include <cmath>

#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/inviscid_deriv.h"
#include "scirocco/core/panel_influence.h"
#include "scirocco/core/panel_influence_deriv.h"
#include "scirocco/core/viscous_geom_deriv.h"
#include "scirocco/core/viscous_residual_deriv.h"
#include "scirocco/core/wake.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

static VectorXd DRVDAlphaImpl(FoilState& state, bool flow, bool wake);

/**
 * @brief compute partial derivatives of cl, cm, cdpi w.r.t. foil.x
 * holding ue (gamma) fixed. Also computes cm_ue and cdpi_ue linearizations.
 * @param state foil state to update with force linearization fields
 */
void CalcForceLin(FoilState& state) {
    double chord = state.geom.chord;
    Vector2d xref = state.geom.xref;
    double Vinf = state.param.Vinf > 0 ? state.param.Vinf : 1.0;
    double alpha = state.oper.alpha;
    int N = state.foil.N;

    /* Get edge velocity (use inviscid ue = gamma for inviscid case) */
    VectorXd ue;
    if (state.oper.viscous && state.glob.U.cols() > 0) {
        ue = state.glob.U.row(3).transpose();
    } else {
        ue = GetUeInv(state);
    }

    /* Compute cp at each node (same as CalcForce) */
    VectorXd cp(N);
    for (int i = 0; i < N; ++i) {
        double cp_inc = 1.0 - (ue(i) / Vinf) * (ue(i) / Vinf);
        if (state.param.KTl > 0) {
            double denom = state.param.KTb + state.param.KTl *
                                                 (1.0 + state.param.KTb) * 0.5 *
                                                 cp_inc;
            if (std::abs(denom) > 1e-10) {
                cp(i) = cp_inc / denom;
            } else {
                cp(i) = cp_inc;
            }
        } else {
            cp(i) = cp_inc;
        }
    }

    /* Initialize partial derivative matrices */
    state.post.dcl_dx_partial = MatrixXd::Zero(2, N);
    state.post.dcm_dx_partial = MatrixXd::Zero(2, N);
    state.post.dcdpi_dx_partial = MatrixXd::Zero(2, N);
    state.post.cm_ue = VectorXd::Zero(N);
    state.post.cdpi_ue = VectorXd::Zero(N);

    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    /* The force integration loops over panel pairs (ip, i) matching CalcForce
     * For each panel:
     *   dx = -(x2(0)-x1(0))*cosa - (x2(1)-x1(1))*sina   (lift direction)
     *   dz = (x2(1)-x1(1))*cosa - (x2(0)-x1(0))*sina    (drag direction)
     *   cl += dx * cpbar / chord
     *   cm += moment terms / chord^2
     *   cdpi += dz * cpbar / chord
     *
     * Partial w.r.t. coordinates (holding cp fixed):
     *   dcl/dx1_x = cpbar * cosa / chord  (from d(-dxv(0))/dx1_x = 1 => no,
     *   dxv=x2-x1) dxv = x2 - x1 => d(dxv)/d(x1) = [-1, 0; 0, -1], d(dxv)/d(x2)
     *   = [1, 0; 0, 1] dx = -dxv(0)*cosa - dxv(1)*sina d(dx)/d(x1_x) = cosa,
     *   d(dx)/d(x1_z) = sina d(dx)/d(x2_x) = -cosa, d(dx)/d(x2_z) = -sina */

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

        /* --- cl partial w.r.t. coordinates ---
         * dcl/d(x1_x) = cpbar * cosa / chord  (from ddx/dx1_x = cosa)
         * dcl/d(x1_z) = cpbar * sina / chord
         * dcl/d(x2_x) = cpbar * (-cosa) / chord
         * dcl/d(x2_z) = cpbar * (-sina) / chord */
        state.post.dcl_dx_partial(0, ip) += cpbar * cosa / chord;
        state.post.dcl_dx_partial(1, ip) += cpbar * sina / chord;
        state.post.dcl_dx_partial(0, i) += cpbar * (-cosa) / chord;
        state.post.dcl_dx_partial(1, i) += cpbar * (-sina) / chord;

        /* --- cdpi partial w.r.t. coordinates ---
         * dz = dxv(1)*cosa - dxv(0)*sina
         * d(dz)/d(x1_x) = sina, d(dz)/d(x1_z) = -cosa
         * d(dz)/d(x2_x) = -sina, d(dz)/d(x2_z) = cosa */
        state.post.dcdpi_dx_partial(0, ip) += cpbar * sina / chord;
        state.post.dcdpi_dx_partial(1, ip) += cpbar * (-cosa) / chord;
        state.post.dcdpi_dx_partial(0, i) += cpbar * (-sina) / chord;
        state.post.dcdpi_dx_partial(1, i) += cpbar * cosa / chord;

        /* --- cm partial w.r.t. coordinates ---
         * cm += (cp1 * dx1nds/3 + cp1 * dx2nds/6 + cp2 * dx1nds/6 + cp2 *
         * dx2nds/3) / chord^2 where dx1nds = dxv(0)*dx1(0) + dxv(1)*dx1(1)
         *       dx2nds = dxv(0)*dx2(0) + dxv(1)*dx2(1)
         *
         * This depends on both dxv and dx1, dx2. The derivatives are:
         * d(dx1nds)/d(x1_x) = -dx1(0) + dxv(0)*(-1) = -dx1(0) - dxv(0)
         * [d(dxv)/dx1_x = -1, d(dx1)/dx1_x = 1... wait] dxv = x2 - x1, dx1 = x1
         * - xref, dx2 = x2 - xref d(dxv(0))/d(x1_x) = -1, d(dx1(0))/d(x1_x) = 1
         * dx1nds = dxv(0)*dx1(0) + dxv(1)*dx1(1)
         * d(dx1nds)/d(x1_x) = (-1)*dx1(0) + dxv(0)*(1) = dxv(0) - dx1(0) */
        double inv_chord2 = 1.0 / (chord * chord);

        /* d(dx1nds)/d(x1_x) = dxv(0) - dx1(0)  ... actually let me be more
         * careful dx1nds = dxv . dx1 d(dx1nds)/d(x1_x) = d(dxv)/d(x1_x) . dx1 +
         * dxv . d(dx1)/d(x1_x)
         *                   = (-1,0) . dx1 + dxv . (1,0) = -dx1(0) + dxv(0) */
        double ddx1nds_dx1x = -dx1(0) + dxv(0);
        double ddx1nds_dx1z = -dx1(1) + dxv(1);
        double ddx1nds_dx2x =
            dx1(0); /* d(dxv)/d(x2_x) = (1,0), d(dx1)/d(x2_x) = 0 */
        double ddx1nds_dx2z = dx1(1);

        /* d(dx2nds)/d(x1_x) = (-1,0).dx2 + dxv.(0,0) = -dx2(0) */
        double ddx2nds_dx1x = -dx2(0);
        double ddx2nds_dx1z = -dx2(1);
        /* d(dx2nds)/d(x2_x) = (1,0).dx2 + dxv.(1,0) = dx2(0) + dxv(0) */
        double ddx2nds_dx2x = dx2(0) + dxv(0);
        double ddx2nds_dx2z = dx2(1) + dxv(1);

        double w1 = cp1 / 3.0 + cp2 / 6.0;
        double w2 = cp1 / 6.0 + cp2 / 3.0;

        state.post.dcm_dx_partial(0, ip) +=
            (w1 * ddx1nds_dx1x + w2 * ddx2nds_dx1x) * inv_chord2;
        state.post.dcm_dx_partial(1, ip) +=
            (w1 * ddx1nds_dx1z + w2 * ddx2nds_dx1z) * inv_chord2;
        state.post.dcm_dx_partial(0, i) +=
            (w1 * ddx1nds_dx2x + w2 * ddx2nds_dx2x) * inv_chord2;
        state.post.dcm_dx_partial(1, i) +=
            (w1 * ddx1nds_dx2z + w2 * ddx2nds_dx2z) * inv_chord2;

        /* --- cm_ue linearization ---
         * cm depends on cp which depends on ue: dcp/due = -2*ue/(Vinf^2) */
        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);

        auto [cp_val_ip, cp_U_ip] = GetCpFromUe(ue(ip), state.param);
        auto [cp_val_i, cp_U_i] = GetCpFromUe(ue(i), state.param);
        double dcp_due_ip = cp_U_ip(3);
        double dcp_due_i = cp_U_i(3);

        state.post.cm_ue(ip) +=
            (dcp_due_ip / 3.0 * dx1nds + dcp_due_ip / 6.0 * dx2nds) *
            inv_chord2;
        state.post.cm_ue(i) +=
            (dcp_due_i / 6.0 * dx1nds + dcp_due_i / 3.0 * dx2nds) * inv_chord2;

        /* --- cdpi_ue linearization --- */
        state.post.cdpi_ue(ip) += dz * 0.5 * dcp_due_ip / chord;
        state.post.cdpi_ue(i) += dz * 0.5 * dcp_due_i / chord;
    }

    /* ========================================================================
     * Viscous output partials: dJ/dV at fixed foil.x, V = U flattened [4*Nsys].
     * - cl, cm, cdpi: only depend on U[3, i] for i ∈ [0, N) (foil edge
     * velocity)
     *   so cl_V[4*i + 3] = cl_ue[i] (similar cm, cdpi). All other entries zero.
     * - cd: Squire-Young from end-of-wake state Uw only; cd_V is sparse at iw.
     * - dcd/dx_partial = 0 because Squire-Young has no direct foil.x
     * dependence.
     * ========================================================================
     */
    if (state.oper.viscous && state.glob.U.cols() > 0) {
        int Nsys = static_cast<int>(state.glob.U.cols());
        state.post.cl_V = VectorXd::Zero(4 * Nsys);
        state.post.cm_V = VectorXd::Zero(4 * Nsys);
        state.post.cd_V = VectorXd::Zero(4 * Nsys);
        for (int i = 0; i < N; ++i) {
            state.post.cl_V(4 * i + 3) = state.post.cl_ue(i);
            state.post.cm_V(4 * i + 3) = state.post.cm_ue(i);
        }

        /* Squire-Young cd from end-of-wake state. */
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

                /* cd = 2 * Uw[0] * q^p
                 * dcd/dUw[k] = 2*δ_k0 * q^p + 2*Uw[0] * d(q^p)/dUw[k]
                 * d(q^p)/dUw[k] = q^p * (dp/dUw[k] * lnq + (p/q) * dq/dUw[k])
                 *               = q^p * (0.5*H_U[k]*lnq +
                 *               (p/q)*(1/Vinf)*uk_U[k]) */
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

/**
 * @brief compute dJ/dx using the adjoint method
 * @param state foil state containing the converged inviscid solution
 * @return Sensitivity struct populated with dcl_dx, dcm_dx, dcdpi_dx, dcd_dx
 */
Sensitivity InviscidSensitivities(FoilState& state) {
    int N = state.foil.N;
    if (N <= 0 || state.isol.AIC.rows() == 0) {
        throw std::runtime_error("InviscidSensitivities: no inviscid solution");
    }

    /* Step 1: Compute force partial derivatives */
    CalcForceLin(state);

    double alpha = state.oper.alpha;

    /* Step 2: Form dJ/dgamma for each objective
     * For inviscid flow, ue(i) = sgnue(i) * gamma(i)
     * cl_ue is dcl/due [N], and due/dgamma(i) = sgnue(i)
     * So dcl/dgamma = cl_ue .* sgnue (elementwise)
     * Similarly for cm, cdpi */

    /* The AIC system is [N+1 x N+1] (includes streamfunction unknown and Kutta)
     * gamma is [N], extended to [N+1] with streamfunction as last entry
     * dJ/dgamma_extended = [dJ/dgamma; 0] (streamfunction doesn't affect forces
     * directly) */

    VectorXd dcl_dgam = VectorXd::Zero(N + 1);
    VectorXd dcm_dgam = VectorXd::Zero(N + 1);
    VectorXd dcdpi_dgam = VectorXd::Zero(N + 1);

    for (int i = 0; i < N; ++i) {
        double sg = state.isol.sgnue(i);
        dcl_dgam(i) = state.post.cl_ue(i) * sg;
        dcm_dgam(i) = state.post.cm_ue(i) * sg;
        dcdpi_dgam(i) = state.post.cdpi_ue(i) * sg;
    }

    /* Step 3: Solve adjoint equations: AIC^T * lambda = dJ/dgamma
     * Note: the residual is R = A*gamma - rhs = 0
     * Adjoint: A^T * lambda = (dJ/dgamma)^T
     * But since R = A*gam - rhs, dR/dgamma = A, so adjoint is A^T * lambda =
     * dJ/dgamma. The sign convention: dJ/dx_total = dJ/dx_partial - lambda^T *
     * dR/dx Since dR = d(A*gam - rhs), we need to subtract. Actually, for R =
     * A*gam - rhs = 0, implicit function theorem gives:
     *   dgam/dx = -A^{-1} * dR/dx  (where dR/dx = dA/dx*gam - drhs/dx)
     *   dJ/dx = dJ/dx_partial + dJ/dgam * dgam/dx
     *         = dJ/dx_partial - dJ/dgam * A^{-1} * dR/dx
     *         = dJ/dx_partial - lambda^T * dR/dx
     * where A^T * lambda = dJ/dgam
     * So: dJ/dx_total = dJ/dx_partial - lambda^T * dR/dx_foil */

    const MatrixXd& AIC = state.isol.AIC;
    MatrixXd AIC_T = AIC.transpose();

    /* Solve for adjoint vectors */
    auto lu = AIC_T.partialPivLu();

    VectorXd lambda_cl = lu.solve(dcl_dgam);
    VectorXd lambda_cm = lu.solve(dcm_dgam);
    VectorXd lambda_cdpi = lu.solve(dcdpi_dgam);

    /* Step 4: Compute VJP: lambda^T * dR/dx_foil */
    VectorXd vjp_cl = AICResidualSensVJP(state, lambda_cl, alpha);
    VectorXd vjp_cm = AICResidualSensVJP(state, lambda_cm, alpha);
    VectorXd vjp_cdpi = AICResidualSensVJP(state, lambda_cdpi, alpha);

    /* Step 5: Total sensitivity = partial - VJP */
    Sensitivity sens;
    sens.dcl_dx = MatrixXd::Zero(2, N);
    sens.dcm_dx = MatrixXd::Zero(2, N);
    sens.dcdpi_dx = MatrixXd::Zero(2, N);
    sens.dcd_dx =
        MatrixXd::Zero(2, N); /* Zero for inviscid (no viscous drag) */

    for (int k = 0; k < N; ++k) {
        sens.dcl_dx(0, k) = state.post.dcl_dx_partial(0, k) - vjp_cl(2 * k);
        sens.dcl_dx(1, k) = state.post.dcl_dx_partial(1, k) - vjp_cl(2 * k + 1);
        sens.dcm_dx(0, k) = state.post.dcm_dx_partial(0, k) - vjp_cm(2 * k);
        sens.dcm_dx(1, k) = state.post.dcm_dx_partial(1, k) - vjp_cm(2 * k + 1);
        sens.dcdpi_dx(0, k) =
            state.post.dcdpi_dx_partial(0, k) - vjp_cdpi(2 * k);
        sens.dcdpi_dx(1, k) =
            state.post.dcdpi_dx_partial(1, k) - vjp_cdpi(2 * k + 1);
    }

    sens.dalpha_dx = MatrixXd::Zero(2, N);
    AlphaSensitivities(state, sens);
    MachSensitivities(state, sens);
    if (state.oper.givencl) ApplyLiftConstraint(sens);
    return sens;
}

void AlphaSensitivities(FoilState& state, Sensitivity& sens) {
    double alpha = state.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    /* dgam/dalpha = gamref[:,0]*(-sin a) + gamref[:,1]*(cos a), times deg2rad
     */
    VectorXd dgam_dalpha =
        (-sina * state.isol.gamref.col(0) + cosa * state.isol.gamref.col(1)) *
        kDegToRad;

    /* dueinv/dalpha = sgnue .* dgam/dalpha */
    VectorXd dueinv_dalpha = state.isol.sgnue.cwiseProduct(dgam_dalpha);

    sens.dcl_dalpha = state.post.cl_alpha + state.post.cl_ue.dot(dueinv_dalpha);
    sens.dcm_dalpha = 0.0 + state.post.cm_ue.dot(dueinv_dalpha);
    sens.dcd_dalpha = 0.0;  // Python mfoil reports zero inviscid cd.
    sens.dcdpi_dalpha =
        state.post.cdpi_alpha + state.post.cdpi_ue.dot(dueinv_dalpha);
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

    double b = param.KTb; /* sqrt(1-Ma^2) */
    double l = param.KTl; /* Ma^2/(1+b)^2 */

    /* db/dMa = -Ma/b */
    double db = -Ma / b;
    /* dl/dMa = 2*Ma/(1+b)^2 - 2*Ma^2/(1+b)^3 * db
     *        = 2*Ma/(1+b)^2 + 2*Ma^3/((1+b)^3 * b) */
    double opb = 1.0 + b;
    double dl =
        2.0 * Ma / (opb * opb) + 2.0 * Ma * Ma * Ma / (opb * opb * opb * b);

    double den = b + 0.5 * l * (1.0 + b) * cp_inc;
    double den_cp = 0.5 * l * (1.0 + b);

    /* cp = cp_inc / den */
    double cp = cp_inc / den;

    /* dden/dMa = db + 0.5*(dl*(1+b) + l*db)*cp_inc */
    double dden_dMa = db + 0.5 * (dl * opb + l * db) * cp_inc;
    /* dden_cp/dMa = 0.5*(dl*(1+b) + l*db) */
    double dden_cp_dMa = 0.5 * (dl * opb + l * db);

    /* dcp/dMa = -cp_inc / den^2 * dden_dMa  (cp_inc doesn't depend on Ma) */
    dcp_dMa = -cp / den * dden_dMa;

    /* d(dcp/due)/dMa: dcp/due = dcp_due_inc * (1 - cp*den_cp) / den
     * d(dcp/due)/dMa = dcp_due_inc * d/dMa[(1 - cp*den_cp)/den] */
    double dfac_dMa = (-(dcp_dMa * den_cp + cp * dden_cp_dMa) * den -
                       (1.0 - cp * den_cp) * dden_dMa) /
                      (den * den);
    d2cp_duedMa = dcp_due_inc * dfac_dMa;
}

/**
 * @brief MachSensitivities: compute dJ/dMa for cl, cd, cm
 * @param state foil state containing the converged inviscid solution
 * @param sens Sensitivity struct to populate with dcl_dMa, dcm_dMa, dcd_dMa
 */
void MachSensitivities(FoilState& state, Sensitivity& sens) {
    int N = state.foil.N;
    double Ma = state.oper.Ma;
    double chord = state.geom.chord;
    double alpha = state.oper.alpha;

    if (Ma <= 0) {
        sens.dcl_dMa = 0.0;
        sens.dcd_dMa = 0.0;
        sens.dcm_dMa = 0.0;
        sens.dcdpi_dMa = 0.0;
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

        double dcp_ip, d2cp_ip;
        double dcp_i, d2cp_i;
        ComputeCpDerivsMa(ue(ip), state.param, Ma, dcp_ip, d2cp_ip);
        ComputeCpDerivsMa(ue(i), state.param, Ma, dcp_i, d2cp_i);

        double dcpbar = 0.5 * (dcp_ip + dcp_i);
        dcl += dx_lift * dcpbar;
        dcdpi += (dxv(1) * cosa - dxv(0) * sina) * dcpbar;

        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
        dcm += dcp_ip / 3.0 * dx1nds + dcp_ip / 6.0 * dx2nds +
               dcp_i / 6.0 * dx1nds + dcp_i / 3.0 * dx2nds;
    }

    sens.dcl_dMa = dcl / chord;
    sens.dcm_dMa = dcm / (chord * chord);
    sens.dcd_dMa = 0.0;  // cdpi is a separate output from cd.
    sens.dcdpi_dMa = dcdpi / chord;
}

/**
 * @brief ViscousSensitivities — adjoint-based dJ/dx_foil for J ∈ {cl, cd, cm}.
 *
 * Steps:
 *   1. CalcForceLin → populates cl_V, cm_V, cd_V (∂J/∂V) and ∂J/∂x_partial.
 *   2. StoreAdjointFactorization → R_V LU at the converged state.
 *   3. For each J: λ_J = R_V^{-T} ∂J/∂V via SolveViscousAdjoint.
 *   4. VJP λ_J^T · ∂R_V/∂x_foil
 *        = λ_bl^T · DRblDFoilX(state, dfoil_s_dx)             [3*Nsys x 2N]
 *        + ComputeRueVJP(state, λ_ue, dfoil_s_dx)             [2N]
 *   5. dJ/dx_foil = ∂J/∂x_partial - vjp_J.
 * @param state foil state containing the converged viscous solution
 * @return Sensitivity struct with dcl_dx, dcm_dx, dcd_dx, and scalar
 * sensitivities
 */
Sensitivity ViscousSensitivities(FoilState& state) {
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;

    if (!state.oper.viscous || state.glob.U.cols() == 0 || !state.glob.conv) {
        throw std::runtime_error(
            "ViscousSensitivities: viscous primal not converged");
    }

    /* 1. Force partials */
    CalcForceLin(state);

    /* 2. Adjoint factorization */
    StoreAdjointFactorization(state);
    if (!state.glob.lu_factorization) {
        throw std::runtime_error(
            "ViscousSensitivities: adjoint LU factorization failed");
    }

    /* 3. Solve λ_J = R_V^{-T} ∂J/∂V for J ∈ {cl, cm, cd} */
    VectorXd lambda_cl = SolveViscousAdjoint(state, state.post.cl_V);
    VectorXd lambda_cm = SolveViscousAdjoint(state, state.post.cm_V);
    VectorXd lambda_cd = SolveViscousAdjoint(state, state.post.cd_V);
    VectorXd cdpi_V = VectorXd::Zero(4 * Nsys);
    for (int i = 0; i < N; ++i) cdpi_V(4 * i + 3) = state.post.cdpi_ue(i);
    VectorXd lambda_cdpi = SolveViscousAdjoint(state, cdpi_V);

    /* 4. VJPs — DRbl as dense matrix * lambda_bl (lambda_bl size 3*Nsys),
     *    plus ComputeRueVJP for the wake-side chain. */
    MatrixXd dfoil_s_dx = DFoilSDFoilX(state.foil);
    MatrixXd dRbl_dx = DRblDFoilX(state, dfoil_s_dx); /* [3*Nsys x 2N] */
    /* The R_ue chain's geometric Jacobians are the same for every output, so
     * build them once rather than once per VJP. */
    const RueVJPContext rue_ctx = BuildRueVJPContext(state, dfoil_s_dx);

    auto vjp = [&](const VectorXd& lambda) -> VectorXd {
        VectorXd lambda_bl = lambda.head(3 * Nsys);
        VectorXd lambda_ue = lambda.tail(Nsys);
        VectorXd out = (lambda_bl.transpose() * dRbl_dx).transpose(); /* [2N] */
        out += ComputeRueVJP(state, lambda_ue, rue_ctx);
        return out;
    };
    VectorXd vjp_cl = vjp(lambda_cl);
    VectorXd vjp_cm = vjp(lambda_cm);
    VectorXd vjp_cd = vjp(lambda_cd);
    VectorXd vjp_cdpi = vjp(lambda_cdpi);

    /* 5. Total dJ/dx = partial - vjp */
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
        sens.dcdpi_dx(0, k) =
            state.post.dcdpi_dx_partial(0, k) - vjp_cdpi(2 * k);
        sens.dcdpi_dx(1, k) =
            state.post.dcdpi_dx_partial(1, k) - vjp_cdpi(2 * k + 1);
        /* dcd/dx_partial is zero for Squire-Young, so this is just -vjp_cd. */
        sens.dcd_dx(0, k) = state.post.dcd_dx_partial(0, k) - vjp_cd(2 * k);
        sens.dcd_dx(1, k) = state.post.dcd_dx_partial(1, k) - vjp_cd(2 * k + 1);
    }

    /* 6. Scalar parameter sensitivities: dJ/dp = ∂J/∂p_partial - λ^T·∂R_V/∂p.
     *
     * DRVDAlpha includes the circulation, wake, and mass-coupling chains.
     * DRVDRe and DRVDMa differentiate the boundary-layer closures and
     * thermodynamic parameters analytically at fixed U.
     *
     * Explicit output dependence is the lift rotation for alpha and the
     * Karman-Tsien pressure/velocity correction for Mach. Reynolds number
     * enters the outputs only through the converged state. */
    {
        VectorXd dRV_dalpha = DRVDAlphaImpl(
            state, true, !state.oper.givencl || state.oper.redowake);
        sens.dcl_dalpha = state.post.cl_alpha - lambda_cl.dot(dRV_dalpha);
        sens.dcm_dalpha = 0.0 - lambda_cm.dot(dRV_dalpha);
        sens.dcd_dalpha = 0.0 - lambda_cd.dot(dRV_dalpha);
        sens.dcdpi_dalpha = state.post.cdpi_alpha - lambda_cdpi.dot(dRV_dalpha);
    }

    /* Ma: ∂J/∂Ma_partial analytically. NO FD anywhere.
     *   cl, cm, cdpi: panel integration with ∂cp/∂Ma via Karman-Tsien
     *     correction at fixed viscous ue. Uses ComputeCpDerivsMa (already
     *     analytic for inviscid MachSensitivities at line 443).
     *   cd: viscous Squire-Young from end-of-wake state. cd =
     *   2·Uw[0]·(uk/Vinf)^p,
     *     p = (5+H)/2. H is Ma-independent at fixed U; uk has Ma dep via
     *     Karman-Tsien lambda. */
    {
        VectorXd dRV_dMa = DRVDMa(state);
        double Ma = state.oper.Ma;
        Vector4d dJ_dMa = Vector4d::Zero();
        if (Ma > 0) {
            /* cl/cm/cdpi: panel integration with ∂cp/∂Ma at fixed ue. */
            int N = state.foil.N;
            double chord = state.geom.chord;
            double alpha = state.oper.alpha;
            double cosa = Cosd(alpha);
            double sina = Sind(alpha);
            Vector2d xref = state.geom.xref;
            VectorXd ue = state.glob.U.row(3).head(N).transpose();

            double dcl = 0, dcm = 0, dcdpi = 0;
            for (int i0 = 1; i0 <= N; ++i0) {
                int i = (i0 == N) ? 0 : i0;
                int ip = (i0 == N) ? N - 1 : i0 - 1;
                Vector2d x1 = state.foil.x.col(ip);
                Vector2d x2 = state.foil.x.col(i);
                Vector2d dxv = x2 - x1;
                Vector2d dx1 = x1 - xref;
                Vector2d dx2 = x2 - xref;
                double dx_lift = -dxv(0) * cosa - dxv(1) * sina;
                double dcp_ip, d2cp_ip, dcp_i, d2cp_i;
                ComputeCpDerivsMa(ue(ip), state.param, Ma, dcp_ip, d2cp_ip);
                ComputeCpDerivsMa(ue(i), state.param, Ma, dcp_i, d2cp_i);
                double dcpbar = 0.5 * (dcp_ip + dcp_i);
                dcl += dx_lift * dcpbar;
                dcdpi += (dxv(1) * cosa - dxv(0) * sina) * dcpbar;
                double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
                double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
                dcm += dcp_ip / 3.0 * dx1nds + dcp_ip / 6.0 * dx2nds +
                       dcp_i / 6.0 * dx1nds + dcp_i / 3.0 * dx2nds;
            }
            dJ_dMa(0) = dcl / chord;
            dJ_dMa(2) = dcm / (chord * chord);
            dJ_dMa(3) = dcdpi / chord;

            /* cd: viscous Squire-Young. Use the same uk_dMa from
             * ComputeThermoMaDeriv chain. */
            if (state.oper.viscous && !state.vsol.Is[2].empty()) {
                int iw = state.vsol.Is[2].back();
                Vector4d Uw = state.glob.U.col(iw);
                auto [H, unused_H_U] = GetH(Uw);
                (void)unused_H_U;
                auto [uk, unused_uk_U] = GetUk(Uw(3), state.param);
                (void)unused_uk_U;
                double Vinf = state.param.Vinf;
                double KTb = state.param.KTb;
                double KTl = state.param.KTl;
                double opb = 1.0 + KTb;
                double KTb_dMa_v = -Ma / KTb;
                double KTl_dMa_v =
                    2.0 * Ma / (opb * opb) -
                    2.0 * Ma * Ma / (opb * opb * opb) * KTb_dMa_v;
                double ue_v = Uw(3);
                double xKT = ue_v / Vinf;
                double den_uk = 1.0 - KTl * xKT * xKT;
                double uk_dMa =
                    ue_v * (xKT * xKT - 1.0) / (den_uk * den_uk) * KTl_dMa_v;
                double q = uk / Vinf;
                double p = 0.5 * (5.0 + H);
                double qp = std::pow(q, p);
                /* ∂cd/∂Ma = 2·Uw[0] · p · q^(p-1) · (uk_dMa/Vinf) */
                dJ_dMa(1) = 2.0 * Uw(0) * p * qp / q * (uk_dMa / Vinf);
            }
        }
        sens.dcl_dMa = dJ_dMa(0) - lambda_cl.dot(dRV_dMa);
        sens.dcd_dMa = dJ_dMa(1) - lambda_cd.dot(dRV_dMa);
        sens.dcm_dMa = dJ_dMa(2) - lambda_cm.dot(dRV_dMa);
        sens.dcdpi_dMa = dJ_dMa(3) - lambda_cdpi.dot(dRV_dMa);
    }

    /* Re: no explicit J(Re) dependence — only chain through R_V. */
    {
        VectorXd dRV_dRe = DRVDRe(state);
        sens.dcl_dRe = -lambda_cl.dot(dRV_dRe);
        sens.dcd_dRe = -lambda_cd.dot(dRV_dRe);
        sens.dcm_dRe = -lambda_cm.dot(dRV_dRe);
        sens.dcdpi_dRe = -lambda_cdpi.dot(dRV_dRe);
    }

    sens.dalpha_dx = MatrixXd::Zero(2, N);
    if (state.oper.givencl) {
        Vector4d seed = Vector4d::Zero();
        if (!state.oper.redowake) {
            const VectorXd rw = DRVDAlphaImpl(state, false, true);
            seed << -lambda_cl.dot(rw), -lambda_cd.dot(rw), -lambda_cm.dot(rw),
                -lambda_cdpi.dot(rw);
        }
        ApplyLiftConstraint(sens, seed);
    }

    return sens;
}

void ApplyLiftConstraint(Sensitivity& sens, const Vector4d& seed) {
    const Vector4d slope(sens.dcl_dalpha, sens.dcd_dalpha, sens.dcm_dalpha,
                         sens.dcdpi_dalpha);
    if (!slope.allFinite() || std::abs(slope(0)) < 1e-10)
        throw std::runtime_error(
            "Lift constraint has singular angle sensitivity");
    sens.dalpha_dx = -sens.dcl_dx / slope(0);
    sens.dalpha_dMa = -sens.dcl_dMa / slope(0);
    sens.dalpha_dRe = -sens.dcl_dRe / slope(0);
    sens.dalpha_dcltarget = 1.0 / slope(0);
    sens.dalpha_dalpha = -seed(0) / slope(0);
    std::array<MatrixXd*, 4> shape{&sens.dcl_dx, &sens.dcd_dx, &sens.dcm_dx,
                                   &sens.dcdpi_dx};
    std::array<double*, 4> ma{&sens.dcl_dMa, &sens.dcd_dMa, &sens.dcm_dMa,
                              &sens.dcdpi_dMa};
    std::array<double*, 4> re{&sens.dcl_dRe, &sens.dcd_dRe, &sens.dcm_dRe,
                              &sens.dcdpi_dRe};
    std::array<double*, 4> alpha{&sens.dcl_dalpha, &sens.dcd_dalpha,
                                 &sens.dcm_dalpha, &sens.dcdpi_dalpha};
    std::array<double*, 4> target{&sens.dcl_dcltarget, &sens.dcd_dcltarget,
                                  &sens.dcm_dcltarget, &sens.dcdpi_dcltarget};
    for (int j = 0; j < 4; ++j) {
        *shape[j] += slope(j) * sens.dalpha_dx;
        *ma[j] += slope(j) * sens.dalpha_dMa;
        *re[j] += slope(j) * sens.dalpha_dRe;
        *alpha[j] = seed(j) + slope(j) * sens.dalpha_dalpha;
        *target[j] = slope(j) / slope(0);
    }
    sens.dcl_dx.setZero();
    sens.dcl_dalpha = sens.dcl_dMa = sens.dcl_dRe = 0;
    sens.dcl_dcltarget = 1;
}

VectorXd DRVDAlpha(FoilState& state) {
    return DRVDAlphaImpl(state, true, true);
}

static VectorXd DRVDAlphaImpl(FoilState& state, bool flow, bool wake) {
    /* Analytic forward-mode chain for ∂R_V/∂α (NO FD).
     *
     * R_V depends on α through several intermediates:
     *   gamma(α)  ──→ ueinv (foil side, wake[0] via continuity)
     *               ──→ xi via sstag(gamma) ──→ R_BL via R_x · dxi/dα
     *   wake.x(α) ──→ ueinv (wake j>0 via V·t)
     *               ──→ ue_sigma (via B wake-source cols, Cgam, Csig)
     *   wake.t(α) ──→ ueinv (V·dt term)
     *               ──→ ue_sigma (via Cgam, Csig wake.t projection)
     *
     * dphi/dα = 0 (sigma_m depends on foil.s/wake.s, both α-indep here; mass
     * = U[1]·U[3] is fixed). So d(ue_sigma · phi)/dα has chain only via
     * ue_sigma. */
    int N = state.foil.N;
    int Nsys = state.glob.Nsys;
    double alpha = state.oper.alpha;
    int Nw = state.wake.N;

    VectorXd dgam_dalpha = flow ? DGammaDAlpha(state) : VectorXd::Zero(N + 1);
    VectorXd dxi_dalpha = DXiDAlpha(state);
    MatrixXd dwake_x_dalpha =
        wake ? DWakeXDAlpha(state) : MatrixXd::Zero(2, Nw); /* [2 x Nw] */
    MatrixXd dwake_t_dalpha =
        wake ? DWakeTDAlpha(state) : MatrixXd::Zero(2, Nw); /* [2 x Nw] */

    VectorXd dR = VectorXd::Zero(4 * Nsys);

    /* ---- R_BL part: R_x · dxi/dα ---- */
    dR.head(3 * Nsys) = state.glob.R_x * dxi_dalpha;

    /* ---- R_ue part: -dueinv/dα - sgue · d(ue_sigma · phi)/dα ----
     *
     * Step 1: dueinv/dα as above. */
    VectorXd dueinv = VectorXd::Zero(Nsys);
    for (int i = 0; i < N; ++i) {
        dueinv(i) = state.isol.sgnue(i) * dgam_dalpha(i);
    }
    dueinv(N) =
        state.isol.sgnue(N - 1) * dgam_dalpha(N - 1); /* wake[0] continuity */
    for (int j = 1; j < Nw; ++j) {
        Vector2d xi_w = state.wake.x.col(j);
        Vector2d t_w = state.wake.t.col(j);
        VelDeriv vd_xi = InviscidVelocityWithDeriv(
            state, xi_w, MatrixXd::Zero(N + 1, 1), alpha);
        Vector2d V_alpha = Vector2d::Zero();
        if (flow) V_alpha = InviscidVelocityDAlpha(state, xi_w, dgam_dalpha);
        Vector2d dwxj = dwake_x_dalpha.col(j);
        Vector2d dV = V_alpha + vd_xi.V_xi * dwxj;
        Vector2d dwtj = dwake_t_dalpha.col(j);
        dueinv(N + j) = t_w.dot(dV) + vd_xi.V.dot(dwtj);
    }

    /* Step 2: ∂(B·phi)/∂α at foil eval points. B has α-dependence ONLY in its
     * wake-source columns (k = N-1..Npanel-1), which depend on wake.x via
     * panel endpoints. Build rhs[i] = sum_{wake panel} ∂(a/b coefficients
     * contribution to B[i,col])/∂α · phi[col], then ∂(Bp·phi)/∂α =
     * -AIC^{-1} · rhs (since Bp = -AIC^{-1}·B and AIC has no α dependence). */
    int Npanel = N + Nw - 2;
    /* Reconstruct phi = sigma_m · mass from state. */
    VectorXd mass(Nsys);
    for (int j = 0; j < Nsys; ++j) {
        mass(j) = state.glob.U(1, j) * state.glob.U(3, j);
    }
    VectorXd phi = state.vsol.sigma_m * mass; /* [Npanel] */

    /* For each wake panel j, the contribution-to-B uses PanelLinSourceStream
     * on left half [wake.x[j], Xm] and right half [Xm, X_right]. The α-deriv
     * of (a, b) flows through dXj endpoints (no xi α-dep since xi = foil.x). */
    auto dwx_col = [&](int j) { return dwake_x_dalpha.col(j); };

    VectorXd rhs_alpha = VectorXd::Zero(N + 1); /* top N rows; Kutta row = 0 */
    for (int j = 0; j < Nw - 1; ++j) {
        Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
        Vector2d dXm = 0.5 * (dwx_col(j) + dwx_col(j + 1));
        Vector2d X_right;
        Vector2d dXr;
        if (j == Nw - 2) {
            X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
            dXr = 2.0 * dwx_col(j + 1) - dXm;
        } else {
            X_right = state.wake.x.col(j + 1);
            dXr = dwx_col(j + 1);
        }
        /* phi-weighting per CalcUeM:
         *   left half: if j>0  → a contributes to phi[N-1+j] (w=0.5) AND
         *   phi[N-1+j-1] (w=0.5)
         *                       b contributes to phi[N-1+j] (w=1.0)
         *              if j==0 → a no contribution (skipped); b contributes to
         *              phi[N-1] (w=1.0)
         *   right half: if j<Nw-2 → a contributes to phi[N-1+j] (w=1.0)
         *                          b contributes to phi[N-1+j] (w=0.5) AND
         *                          phi[N-1+j+1] (w=0.5)
         *               if j==Nw-2 → a contributes to phi[N-1+j] (w=1.0)
         *                          b contributes to phi[N-1+j] (w=1.0) */
        double wL_a, wL_b, wR_a, wR_b;
        double phi_jm1 = (j > 0) ? phi(N - 1 + j - 1) : 0.0;
        double phi_j = phi(N - 1 + j);
        double phi_jp1 = (j < Nw - 2) ? phi(N - 1 + j + 1) : 0.0;
        if (j > 0) {
            wL_a = 0.5 * phi_j + 0.5 * phi_jm1;
            wL_b = phi_j;
        } else {
            wL_a = 0.0;
            wL_b = phi_j;
        }
        if (j < Nw - 2) {
            wR_a = phi_j;
            wR_b = 0.5 * phi_j + 0.5 * phi_jp1;
        } else {
            wR_a = phi_j;
            wR_b = phi_j;
        }
        Matrix2d Xj_L;
        Xj_L.col(0) = state.wake.x.col(j);
        Xj_L.col(1) = Xm;
        Matrix2d Xj_R;
        Xj_R.col(0) = Xm;
        Xj_R.col(1) = X_right;
        /* Endpoint α-derivs: */
        Vector2d dXjL_0 = dwx_col(j);
        const Vector2d& dXjL_1 = dXm;
        const Vector2d& dXjR_0 = dXm;
        Vector2d dXjR_1 = dXr;
        for (int i = 0; i < N; ++i) {
            Vector2d xi = state.foil.x.col(i);
            LinSourceStreamDeriv dL = PanelLinSourceStreamDeriv(Xj_L, xi);
            LinSourceStreamDeriv dR_p = PanelLinSourceStreamDeriv(Xj_R, xi);
            /* ∂a/∂α = a_Xj[0:2]·dXjL_0 + a_Xj[2:4]·dXjL_1 (no xi α-dep) */
            auto chain = [&](const Vector4d& gXj, const Vector2d& dX0,
                             const Vector2d& dX1) {
                return gXj(0) * dX0(0) + gXj(1) * dX0(1) + gXj(2) * dX1(0) +
                       gXj(3) * dX1(1);
            };
            double daL = chain(dL.a_Xj, dXjL_0, dXjL_1);
            double dbL = chain(dL.b_Xj, dXjL_0, dXjL_1);
            double daR = chain(dR_p.a_Xj, dXjR_0, dXjR_1);
            double dbR = chain(dR_p.b_Xj, dXjR_0, dXjR_1);
            rhs_alpha(i) += wL_a * daL + wL_b * dbL + wR_a * daR + wR_b * dbR;
        }
    }
    /* ∂(Bp · phi)/∂α = -AIC^{-1} · rhs_alpha. AIC is [N+1 × N+1]. */
    VectorXd dBphi_dalpha(N + 1);
    {
        const MatrixXd& AIC = state.isol.AIC;
        dBphi_dalpha = -AIC.partialPivLu().solve(rhs_alpha);
    }

    /* Step 3: ue_sigma · phi α-derivative.
     *   foil row i (i < N): (ue_sigma · phi)[i] = (Bp · phi)[i] =
     *   dBphi_dalpha[i]. wake row N+0: (Dw · phi)[0] = (Bp[N-1] · phi) by
     *   continuity override. wake row N+j (j > 0): (Dw · phi)[j] = (Cgam[j,:] ·
     *   Bp · phi) +
     *                                          (Csig[j,:] · phi).
     *
     * Build v = Bp · phi (length N) from existing source_term: source_term =
     * ue_m · mass = sgue · ue_sigma · phi. So ue_sigma · phi =
     * source_term/sgue. */
    VectorXd source_term = state.vsol.ue_m * mass;
    VectorXd v(N);
    for (int k = 0; k < N; ++k) {
        v(k) = source_term(k) / state.isol.sgnue(k);
    }
    VectorXd v_pad = VectorXd::Zero(N + 1);
    v_pad.head(N) = v;

    /* For the wake rows j > 0, compute:
     *   d(Cgam[j,:] · v)/dα  +  Cgam[j,:] · dBphi_dalpha.head(N)
     *                        +  d(Csig[j,:] · phi)/dα
     *
     * d(Cgam[j,:] · v)/dα — Cgam[j, k] = V_G_k(wake.x[j]) · wake.t[j].
     * V_G is pure-geometry (no α-dep at fixed xi), so:
     *   d(Cgam[j,:] · v) = (V_G · v)_xi · dwxj  · t_w   + (V_G · v) · dwtj
     * where (V_G · v) at xi=wake.x[j] is obtained via
     * InviscidVelocityWithDerivExt with vortex=v (treating v as gammas), then
     * subtracting freestream.
     *
     * d(Csig[j,:] · phi)/dα — Csig involves wake-source panel velocities at
     * wake.x[j] projected onto wake.t[j]. Chain through both eval-point
     * dwxj (via V_xi) AND source-panel endpoints (chain through panel-velocity
     * derivative). */
    VectorXd dCgam_v_dalpha(Nw);
    VectorXd dCsig_phi_dalpha(Nw);
    dCgam_v_dalpha.setZero();
    dCsig_phi_dalpha.setZero();
    for (int jw = 1; jw < Nw; ++jw) {
        Vector2d xi_w = state.wake.x.col(jw);
        Vector2d t_w = state.wake.t.col(jw);
        Vector2d dwxj = dwake_x_dalpha.col(jw);
        Vector2d dwtj = dwake_t_dalpha.col(jw);

        /* (V_G · v) at wake.x[jw]: use InviscidVelocityWithDerivExt with
         * vortex strengths = v_pad, dvortex = 0. Returns V_full = V_G·v + V_∞;
         * subtract V_∞ to get V_G·v. */
        MatrixXd dvortex_zero = MatrixXd::Zero(N + 1, 1);
        VelDeriv vd_v = InviscidVelocityWithDerivExt(state, xi_w, v_pad,
                                                     dvortex_zero, alpha);
        Vector2d Vinf(state.oper.Vinf * Cosd(alpha),
                      state.oper.Vinf * Sind(alpha));
        Vector2d VGv = vd_v.V - Vinf; /* V_G · v at wake.x[jw] */
        Matrix2d VGv_xi = vd_v.V_xi;  /* (V_G · v)_xi (Vinf is xi-indep) */
        dCgam_v_dalpha(jw) = t_w.dot(VGv_xi * dwxj) + VGv.dot(dwtj);

        /* dCsig[j,:]·phi/dα: chain through wake-source panel velocities
         * dotted with t_w. We replicate the Csig assembly (CalcUeM lines
         * 824-911) but in forward mode w.r.t. α. Each wake source panel
         * contributes via (a, b) Vector2d evaluated at eval point xi=wake.x[jw]
         * and projected on t_w; the α-derivative comes from xi-deriv chain
         * (dwxj), t_w chain (dwtj), and source-endpoint chain (dwx at the
         * panel's endpoints). */
        for (int wp = 0; wp < Nw; ++wp) {
            if (wp == jw) {
                /* Self-influence (eval point on the same wake panel as source).
                 * Three cases per CalcUeM (mfoil.py:1543-1551):
                 *   interior (1 ≤ wp ≤ Nw-2): aa = (0.25/kPi)·log(d1/d2) where
                 *     d1 = |Xm - Xl|, d2 = |Xr - Xm|, with Xm = wake.x[wp],
                 *     Xl = 0.5(wake.x[wp-1]+wake.x[wp]), Xr = 0.5(wake.x[wp]+
                 *     wake.x[wp+1]). Csig contributes aa to two phi cols.
                 *   wp == 0 (TE three-panel system): unreachable here — the
                 *     outer loop starts at jw=1, and wake row 0 of ue_sigma is
                 *     overridden by Bp[N-1,:] for continuity (handled at
                 *     line ~682 via dgam_dalpha), so Csig[0,:] never enters
                 *     the chain.
                 *   wp == Nw-1 (ghost-extended last panel): identically zero
                 *     (mfoil.py:1550-1551, "no self effect"). */
                if (wp >= 1 && wp <= Nw - 2) {
                    Vector2d Xl =
                        0.5 * (state.wake.x.col(wp - 1) + state.wake.x.col(wp));
                    Vector2d Xm = state.wake.x.col(wp);
                    Vector2d Xr =
                        0.5 * (state.wake.x.col(wp) + state.wake.x.col(wp + 1));
                    Vector2d dXl = 0.5 * (dwx_col(wp - 1) + dwx_col(wp));
                    Vector2d dXm = dwx_col(wp);
                    Vector2d dXr = 0.5 * (dwx_col(wp) + dwx_col(wp + 1));
                    Vector2d e1 = Xm - Xl;
                    Vector2d e2 = Xr - Xm;
                    double d1 = e1.norm();
                    double d2 = e2.norm();
                    if (d1 > 1e-14 && d2 > 1e-14) {
                        double dd1 = (e1 / d1).dot(dXm - dXl);
                        double dd2 = (e2 / d2).dot(dXr - dXm);
                        double daa = (0.25 / kPi) * (dd1 / d1 - dd2 / d2);
                        int col_prev = N - 1 + wp - 1;
                        int col_curr = N - 1 + wp;
                        dCsig_phi_dalpha(jw) +=
                            daa * (phi(col_prev) + phi(col_curr));
                    }
                }
                /* wp == 0 unreachable (see above); wp == Nw-1 contributes 0. */
                continue;
            }
            /* Non-self contributions: build panel velocity derivatives.
             * Per CalcUeM lines 877-908:
             *   wp == 0: half panel Xj = [wake.x[0],
             *   0.5*(wake.x[0]+wake.x[1])].
             *            v.dot(ti) → cols (N-1) [from b]; (0), (N-2) [from a].
             *   wp == Nw-1: const panel Xj = [Xl, Xr] (with Xr = 2*Xm - Xl,
             *                Xm = wake.x[wp], Xl =
             *                0.5(wake.x[wp-1]+wake.x[wp])). v.dot(ti) → col
             *                (Npanel-1).
             *   interior wp: left half [Xl, Xm] + right half [Xm, Xr].
             *                left a → col (N-1+wp-1) (w=1), b → col (N-1+wp-1)
             *                (w=0.5) right a → col (N-1+wp) (w=0.5), b → col
             *                (N-1+wp) (w=1) */
            if (wp == 0) {
                Vector2d Xm0 =
                    0.5 * (state.wake.x.col(0) + state.wake.x.col(1));
                Vector2d dXm0 = 0.5 * (dwx_col(0) + dwx_col(1));
                Matrix2d Xj;
                Xj.col(0) = state.wake.x.col(0);
                Xj.col(1) = Xm0;
                Vector2d dXj0 = dwx_col(0);
                Vector2d dXj1 = dXm0;
                LinSourceVelocityDeriv lvd =
                    PanelLinSourceVelocityDeriv(Xj, xi_w);
                /* Endpoint chain: ∂a/∂α (and ∂b/∂α) via Xj endpoints + xi. */
                auto da_alpha = [&]() {
                    return lvd.a_Xj.col(0) * dXj0(0) +
                           lvd.a_Xj.col(1) * dXj0(1) +
                           lvd.a_Xj.col(2) * dXj1(0) +
                           lvd.a_Xj.col(3) * dXj1(1) + lvd.a_xi * dwxj;
                };
                auto db_alpha = [&]() {
                    return lvd.b_Xj.col(0) * dXj0(0) +
                           lvd.b_Xj.col(1) * dXj0(1) +
                           lvd.b_Xj.col(2) * dXj1(0) +
                           lvd.b_Xj.col(3) * dXj1(1) + lvd.b_xi * dwxj;
                };
                Vector2d da = da_alpha();
                Vector2d db = db_alpha();
                /* dCsig·phi contributions (per CalcUeM lines 877-885):
                 *   col(N-1) += b.dot(ti)*phi[N-1]
                 *   col(0) += a.dot(ti)*phi[0]
                 *   col(N-2) += a.dot(ti)*phi[N-2] */
                double w_a = phi(0) + phi(N - 2);
                double w_b = phi(N - 1);
                dCsig_phi_dalpha(jw) += w_a * (da.dot(t_w) + lvd.a.dot(dwtj)) +
                                        w_b * (db.dot(t_w) + lvd.b.dot(dwtj));
            } else if (wp == Nw - 1) {
                Vector2d Xl =
                    0.5 * (state.wake.x.col(wp - 1) + state.wake.x.col(wp));
                Vector2d Xm = state.wake.x.col(wp);
                Vector2d Xr = 2.0 * Xm - Xl;
                Vector2d dXl = 0.5 * (dwx_col(wp - 1) + dwx_col(wp));
                Vector2d dXm = dwx_col(wp);
                Vector2d dXr = 2.0 * dXm - dXl;
                Matrix2d Xj;
                Xj.col(0) = Xl;
                Xj.col(1) = Xr;
                ConstSourceVelocityDeriv csv =
                    PanelConstSourceVelocityDeriv(Xj, xi_w);
                Vector2d dV = csv.val_Xj.col(0) * dXl(0) +
                              csv.val_Xj.col(1) * dXl(1) +
                              csv.val_Xj.col(2) * dXr(0) +
                              csv.val_Xj.col(3) * dXr(1) + csv.val_xi * dwxj;
                double w = phi(Npanel - 1);
                dCsig_phi_dalpha(jw) += w * (dV.dot(t_w) + csv.val.dot(dwtj));
            } else {
                /* Interior wp ∈ [1, Nw-2], wp != jw */
                Vector2d Xl =
                    0.5 * (state.wake.x.col(wp - 1) + state.wake.x.col(wp));
                Vector2d Xm = state.wake.x.col(wp);
                Vector2d Xr =
                    0.5 * (state.wake.x.col(wp) + state.wake.x.col(wp + 1));
                Vector2d dXl = 0.5 * (dwx_col(wp - 1) + dwx_col(wp));
                Vector2d dXm = dwx_col(wp);
                Vector2d dXr = 0.5 * (dwx_col(wp) + dwx_col(wp + 1));
                /* Left half [Xl, Xm] */
                {
                    Matrix2d Xj_L;
                    Xj_L.col(0) = Xl;
                    Xj_L.col(1) = Xm;
                    LinSourceVelocityDeriv lvd =
                        PanelLinSourceVelocityDeriv(Xj_L, xi_w);
                    Vector2d da = lvd.a_Xj.col(0) * dXl(0) +
                                  lvd.a_Xj.col(1) * dXl(1) +
                                  lvd.a_Xj.col(2) * dXm(0) +
                                  lvd.a_Xj.col(3) * dXm(1) + lvd.a_xi * dwxj;
                    Vector2d db = lvd.b_Xj.col(0) * dXl(0) +
                                  lvd.b_Xj.col(1) * dXl(1) +
                                  lvd.b_Xj.col(2) * dXm(0) +
                                  lvd.b_Xj.col(3) * dXm(1) + lvd.b_xi * dwxj;
                    double phi_prev = phi(N - 1 + wp - 1);
                    dCsig_phi_dalpha(jw) +=
                        phi_prev * (da.dot(t_w) + lvd.a.dot(dwtj)) +
                        0.5 * phi_prev * (db.dot(t_w) + lvd.b.dot(dwtj));
                }
                /* Right half [Xm, Xr] */
                {
                    Matrix2d Xj_R;
                    Xj_R.col(0) = Xm;
                    Xj_R.col(1) = Xr;
                    LinSourceVelocityDeriv lvd =
                        PanelLinSourceVelocityDeriv(Xj_R, xi_w);
                    Vector2d da = lvd.a_Xj.col(0) * dXm(0) +
                                  lvd.a_Xj.col(1) * dXm(1) +
                                  lvd.a_Xj.col(2) * dXr(0) +
                                  lvd.a_Xj.col(3) * dXr(1) + lvd.a_xi * dwxj;
                    Vector2d db = lvd.b_Xj.col(0) * dXm(0) +
                                  lvd.b_Xj.col(1) * dXm(1) +
                                  lvd.b_Xj.col(2) * dXr(0) +
                                  lvd.b_Xj.col(3) * dXr(1) + lvd.b_xi * dwxj;
                    double phi_curr = phi(N - 1 + wp);
                    dCsig_phi_dalpha(jw) +=
                        0.5 * phi_curr * (da.dot(t_w) + lvd.a.dot(dwtj)) +
                        phi_curr * (db.dot(t_w) + lvd.b.dot(dwtj));
                }
            }
        }
        /* Foil-source panel α-deriv (in Csig): wp = 0..N-2.
         * Csig[jw, wp] += v.dot(ti) where v = PanelConstSourceVelocity(Xj_foil,
         * xi=wake.x[jw], &ti). α-chain via dwxj (eval point) and dwtj
         * (projection); panel endpoints are foil.x (α-indep). */
        int jstart = (jw == 0) ? 1 : 0;
        int jend = (jw == 0) ? N - 2 : N - 1;
        for (int jf = jstart; jf < jend; ++jf) {
            Matrix2d Xj_foil;
            Xj_foil.col(0) = state.foil.x.col(jf);
            Xj_foil.col(1) = state.foil.x.col(jf + 1);
            ConstSourceVelocityDeriv csv =
                PanelConstSourceVelocityDeriv(Xj_foil, xi_w);
            Vector2d dV = csv.val_xi * dwxj;
            dCsig_phi_dalpha(jw) += phi(jf) * (dV.dot(t_w) + csv.val.dot(dwtj));
        }
    }

    /* Step 4: Assemble R_ue α-derivative.
     *   foil rows i < N:     dR_ue/dα = -dueinv - sgnue[i] · dBphi_dalpha[i]
     *   wake row N+0:        dR_ue/dα = -dueinv - 1 · dBphi_dalpha[N-1]
     *   wake rows N+j j>0:   dR_ue/dα = -dueinv - 1 · (dCgam·v + Cgam·dBphi +
     *   dCsig·phi) */
    for (int i = 0; i < N; ++i) {
        dR(3 * Nsys + i) = -dueinv(i) - state.isol.sgnue(i) * dBphi_dalpha(i);
    }
    dR(3 * Nsys + N) =
        -dueinv(N) - dBphi_dalpha(N - 1); /* wake[0] continuity */
    for (int j = 1; j < Nw; ++j) {
        /* Cgam[j,:] · dBphi.head(N): rebuild Cgam[j,:] from stored Cgam? It's
         * not cached. Recompute via InviscidVelocity returning V_G [2 x N]. */
        Vector2d xi_w = state.wake.x.col(j);
        Vector2d t_w = state.wake.t.col(j);
        MatrixXd V_G;
        InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi_w, true,
                         &V_G);
        double Cgam_dBphi = 0.0;
        for (int k = 0; k < N; ++k) {
            double Cgam_jk = V_G(0, k) * t_w(0) + V_G(1, k) * t_w(1);
            Cgam_dBphi += Cgam_jk * dBphi_dalpha(k);
        }
        double d_ue_sigma_phi =
            dCgam_v_dalpha(j) + Cgam_dBphi + dCsig_phi_dalpha(j);
        dR(3 * Nsys + N + j) = -dueinv(N + j) - d_ue_sigma_phi;
    }
    return dR;
}

VectorXd DRVDMa(FoilState& state, double h) {
    /* Analytic ∂R_V/∂Ma via the thermo-param chain (KTb, KTl, H0, mu0, rho0)
     * through every Ma-dependent closure in the BL residuals. BuildGlobSys
     * accumulates per-station ∂R/∂Ma·Ma into glob.dR_dMa_times_Ma; the
     * chain finishes with ∂R/∂Ma = dR_dMa_times_Ma / Ma. ue-equation rows
     * are Ma-independent (ueinv, sigma_m, ue_m, mass all Ma-indep at fixed U).
     */
    (void)h;
    int Nsys = state.glob.Nsys;
    state.glob.realloc = true;
    BuildGlobSys(state);
    VectorXd dR = VectorXd::Zero(4 * Nsys);
    double Ma = state.oper.Ma;
    if (Ma > 0) {
        dR.head(3 * Nsys) = state.glob.dR_dMa_times_Ma / Ma;
    }
    return dR;
}

VectorXd DRVDRe(FoilState& state, double h) {
    /* Analytic ∂R_V/∂Re via the Ret chain. Re only enters R_V through
     * Ret = ρ·ue·θ/μ_0 in BL closures, with ∂Ret/∂Re = Ret/Re. BuildGlobSys
     * accumulates per-station `lin_Ret` = sum_k ∂R_BL/∂Ret_k · Ret_k into
     * glob.dR_dRe_times_Re, so ∂R_BL/∂Re = dR_dRe_times_Re / Re. The
     * ue-equation rows have no Re dependence (ue, ueinv, ue_m all Re-indep). */
    (void)h;
    int Nsys = state.glob.Nsys;
    state.glob.realloc = true;
    BuildGlobSys(state);
    VectorXd dR = VectorXd::Zero(4 * Nsys);
    dR.head(3 * Nsys) = state.glob.dR_dRe_times_Re / state.oper.Re;
    return dR;
}

/**
 * @brief SolveViscousAdjoint — transpose solve through the converged viscous
 * R_V.
 *
 * R_V is assembled and LU-factored by StoreAdjointFactorization. Eigen's
 * SparseLU supports transpose solves via the `.transpose()` view, which is
 * what we use here.
 * @param state foil state with a valid adjoint LU factorization
 * @param rhs right-hand side vector of size 4*Nsys
 * @return adjoint solution vector λ of size 4*Nsys
 */
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
    VectorXd lambda = state.glob.lu_factorization->transpose().solve(rhs);
    VectorXd residual = state.glob.R_V.transpose() * lambda - rhs;
    // Check the original equations, never a regularized replacement.
    VectorXd scale = state.glob.R_V.cwiseAbs().transpose() * lambda.cwiseAbs();
    scale += rhs.cwiseAbs();
    if (!lambda.allFinite() || !residual.allFinite() ||
        residual.norm() > 1e-11 * std::max(scale.norm(), 1e-30)) {
        throw std::runtime_error(
            "SolveViscousAdjoint: inaccurate transpose solve");
    }
    return lambda;
}

}  // namespace scirocco
