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
*/
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

/* ============================================================================
 * Phase 1 geometric chain functions for the viscous adjoint.
 *
 * These provide d(geom_quantity)/d(foil.x) for the geometric quantities that
 * the viscous residual R_V depends on (xi, wake.x, wake.s, wgap, sstag).
 *
 * The geometric chain takes dfoil_s_dx [N x 2N], with column 2*k+dim
 * corresponding to node k, coordinate dim. DFoilSDFoilX selects the spline
 * derivative used by RefreshPanelGeometry, cumulative panel norms selected
 * by RecomputeArcLength, or zero for fixed constructor-time arc labels. The
 * wake-coordinate and tangent overloads accept the same Jacobian so the entire
 * chain uses one convention.
 *
 * Column layout for output matrices: column index = 2*k + dim, where
 *   k   ∈ [0, N)  is the foil node index
 *   dim ∈ {0, 1} is x (0) or y (1).
 * ============================================================================
 */

/**
 * @brief Compute d(foil.s)/d(foil.x) using cumulative panel norms (matches
 * scirocco's RecomputeArcLength). Returns matrix of shape [N x 2N].
 *
 * This is only the derivative of RecomputeArcLength. For spline arc lengths
 * from RefreshPanelGeometry, use DFoilSDFoilX or SplineArcLengthDeriv.
 *
 * @param foil Panel object containing foil geometry.
 * @return MatrixXd of shape [N x 2N].
 */
MatrixXd DFoilSDFoilXCumulative(const Panel& foil);

/** @brief Arc-length Jacobian for the panel's declared geometry convention. */
MatrixXd DFoilSDFoilX(const Panel& foil);

/**
 * @brief d(sstag)/d(foil.x) — stagnation point arc-length, gamma-based
 * interpolation. Returns VectorXd of length 2N.
 *
 * Chain pieces:
 *   1. d(sstag)/d(foil.s[Istag]) via interpolation weights (w1, w2)
 *   2. d(sstag)/d(gamma[Istag]) via sstag_g (already cached in state.isol)
 *      → reverse-mode adjoint of AIC: λ_i = (AIC^T)^{-1} e_{Istag[i]}
 *      → d(gamma_Istag[i])/d(foil.x) = -AICResidualSensVJP(λ_i, alpha)
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return VectorXd of length 2N.
 */
VectorXd DSstagDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief d(wake.s)/d(foil.x) — cumulative wake arc-length. Returns [Nw x 2N].
 *
 * wake.s[i] = foil.s[N-1] + sv[i]
 *   sv[0] = 0
 *   sv[i] depends on ds1 = 0.5*((foil.s[1]-foil.s[0]) +
 *   (foil.s[N-1]-foil.s[N-2])) and on chord (constant in foil.x → contributes
 *   0).
 *
 * Chain:
 *   d(wake.s[i])/d(foil.x) = d(foil.s[N-1])/d(foil.x) + d(sv[i])/d(ds1) *
 *   d(ds1)/d(foil.s) * d(foil.s)/d(foil.x)
 *
 * This is the chain the previous attempt missed (per memory:
 * "Missing d(wake.s)/d(foil.x) chain causes ~0.216 adjoint error at TE nodes").
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nw x 2N].
 */
MatrixXd DWakeSDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief d(xi)/d(foil.x) — combined for ALL Nsys stations (foil + wake).
 * Returns [Nsys x 2N] where Nsys = N + Nw.
 *
 * Foil-side (i ∈ [0, N)):
 *   xi[i] = |foil.s[i] - sstag|
 *   d(xi[i])/d(foil.x) = sign(foil.s[i] - sstag) *
 *                        (d(foil.s[i])/d(foil.x) - d(sstag)/d(foil.x))
 *
 * Wake-side (i ∈ [N, N+Nw)):
 *   xi[N+j] = wake.s[j] - sstag
 *   d(xi[N+j])/d(foil.x) = d(wake.s[j])/d(foil.x) - d(sstag)/d(foil.x)
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nsys x 2N].
 */
MatrixXd DXiDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief d(wgap)/d(foil.x) — wake gap from TE blunt-base extrapolation.
 * Returns [Nw x 2N].
 *
 * wgap[i] depends on:
 *   - hTE = |TE gap|  → depends on foil.x[0, N-1] only
 *   - dtdx = TE thickness slope → depends on foil.x[{0,1,N-2,N-1}]
 *   - tcp, tdp factors → depends on same TE nodes
 *   - xi[N+i] - xi[N] (wake-only arc-length from TE) → depends on
 *     wake.s[i] - foil.s[N-1] via cancellation with sstag
 *
 * The cubic formula is:
 *   wgap[i] = hTE * (1 + (2 + flen*dtdx)*xib) * (1-xib)^2,  xib = (xi-xi0)/Lw
 *   Lw = flen*hTE
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nw x 2N].
 */
MatrixXd DWgapDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief d(wake.x)/d(foil.x) — wake node coordinates from BuildWake's
 * predictor- corrector streamline tracing. Returns [2*Nw x 2N], where row d*Nw
 * + j = d(wake.x(d, j))/d(foil.x).
 *
 * Implementation: linearized predictor-corrector sweep matching BuildWake.
 * The initial wake.x[:,0] = x_te + 1e-5·chord·t_te depends on foil.x[:,0]
 * and foil.x[:,N-1] via the TE midpoint and the (n_te-rotated) unit tangent.
 * Each subsequent xyw[:,i+1] depends on xyw[:,i] and on
 * InviscidVelocityWithDeriv at xyw[:,i] and at the predictor point, where
 * the velocity sensitivity has both a parametric piece (V_xfoil via dgamma)
 * and a spatial piece (V_xi chained through dxyw_i). Step size ds_i =
 * sv(i+1)−sv(i) is chained through SpaceGeomWithDeriv → ds1 → foil.s.
 * Validated against geom_chain_reference.h5::dwake_x_dx (CS at 1e-30 step)
 * to 1e-10 rel / 2e-11 abs across all 8484 entries in viscous_chain_test.cpp.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2*Nw x 2N].
 */
MatrixXd DWakeXDFoilX(const FoilState& state);
MatrixXd DWakeXDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief DWakeXDFoilX with a precomputed d(gamma)/d(foil.x).
 *
 * Building d(gamma)/d(foil.x) dominates the cost of this function. Callers
 * that need several wake Jacobians for the same state compute it once with
 * DGammaDFoilXForWake and pass it here.
 *
 * @param dgam_dx_wake DGammaDFoilXForWake(state) [N+1 x 2N].
 */
MatrixXd DWakeXDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx,
                      const MatrixXd& dgam_dx_wake);

/**
 * @brief d(gamma)/d(foil.x) at the alpha the wake was traced with.
 *
 * Equal to DGammaDFoilX(state, state.oper.alpha) unless the wake was traced
 * at another alpha (state.wake.tracing_alpha), as in lift-target solves that
 * keep the wake frozen. The wake Jacobians need this version.
 *
 * @return MatrixXd of shape [N+1 x 2N].
 */
MatrixXd DGammaDFoilXForWake(const FoilState& state);

/**
 * @brief True when the wake was traced at a different alpha than
 * state.oper.alpha (lift-target solves with a frozen wake, or a warm restart
 * that kept the old wake).
 */
bool WakeTracedAtOtherAlpha(const FoilState& state);

/**
 * @brief d(wake.t)/d(foil.x) — wake unit-tangent vectors from BuildWake's
 * predictor- corrector sweep. Returns [2*Nw x 2N], same row layout as
 * DWakeXDFoilX (row d*Nw + j = d(wake.t(d, j))/d(foil.x)).
 *
 * For j < Nw-1, wake.t[j] = V/|V| at wake.x[j] (the v1 in iteration i=j).
 * For j == Nw-1, wake.t[Nw-1] = V/|V| at the PREDICTOR point x_pred from
 * iteration i=Nw-2 (i.e. not at wake.x[Nw-1] itself).
 *
 * The simple "(I - t t^T)/|V| * dV at wake.x[j]" formula is exact for
 * j < Nw-1 but wrong for j = Nw-1 because t is set from a different point
 * than wake.x[Nw-1]. This function does the iterative chain to get it right.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2*Nw x 2N].
 */
MatrixXd DWakeTDFoilX(const FoilState& state);
MatrixXd DWakeTDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief DWakeTDFoilX with a precomputed d(gamma)/d(foil.x).
 *
 * @param dgam_dx_wake DGammaDFoilXForWake(state) [N+1 x 2N].
 */
MatrixXd DWakeTDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx,
                      const MatrixXd& dgam_dx_wake);

/* ----- Alpha-direction (scalar) versions of the above. Used by the analytic
 * DRVDAlpha chain. Each returns a SINGLE-column derivative since α is scalar.
 */

/**
 * @brief d(gamma)/d(alpha) [N+1] — closed-form from gamref[:,0]·cos α +
 * gamref[:,1]·sin α.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length N+1.
 */
VectorXd DGammaDAlpha(const FoilState& state);

/**
 * @brief d(wake.x)/d(alpha) [2 x Nw] — mirrors DWakeXDFoilX's
 * predictor-corrector sweep with dgam_dα (length N+1) and direct freestream
 * alpha contribution in the velocity.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2 x Nw].
 */
MatrixXd DWakeXDAlpha(const FoilState& state);

/**
 * @brief d(wake.t)/d(alpha) [2 x Nw] — mirrors DWakeTDFoilX.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2 x Nw].
 */
MatrixXd DWakeTDAlpha(const FoilState& state);

/**
 * @brief d(wake.s)/d(alpha) [Nw] — cumulative panel norms along wake.x
 * perturbation.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length Nw.
 */
VectorXd DWakeSDAlpha(const FoilState& state);

/**
 * @brief d(sstag)/d(alpha) — scalar. Chain through dgamma at the stagnation
 * panel.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return Scalar derivative of sstag w.r.t. alpha.
 */
double DSstagDAlpha(const FoilState& state);

/**
 * @brief d(xi)/d(alpha) [Nsys] — chain through sstag (foil-side) and wake.s
 * (wake-side).
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length Nsys.
 */
VectorXd DXiDAlpha(const FoilState& state);

/**
 * @brief d(wgap)/d(alpha) [Nw] — chain through wake.s and xi (wgap's xib).
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length Nw.
 */
VectorXd DWgapDAlpha(const FoilState& state);

/**
 * @brief InviscidVelocity α-derivative: ∂V/∂α [2] at xi, with the chain through
 * V_G · dgam/dα (gamma part) plus freestream V_∞ ·(-sin α, cos α)·deg2rad.
 *
 * @param state FoilState containing the current solution and geometry.
 * @param xi Evaluation point in 2D space.
 * @param dgam_dalpha d(gamma)/d(alpha) vector of length N+1.
 * @return Vector2d of ∂V/∂α.
 */
Vector2d InviscidVelocityDAlpha(const FoilState& state, const Vector2d& xi,
                                const VectorXd& dgam_dalpha);

/**
 * @brief d(TE info)/d(foil.x[k]_dim): expose the static helper from
 * adjoint.cpp. Only nodes k ∈ {0, 1, N-2, N-1} produce nonzero output.
 *
 * @param foil_x Foil node coordinate matrix [2 x N].
 * @param k Node index (only k ∈ {0, 1, N-2, N-1} produce nonzero output).
 * @param dim Coordinate dimension: 0 for x, 1 for y.
 * @param dtcp Output: d(tcp)/d(foil.x[k]_dim).
 * @param dtdp Output: d(tdp)/d(foil.x[k]_dim).
 * @param dhTE Output: d(hTE)/d(foil.x[k]_dim).
 * @param ddtdx Output: d(dtdx)/d(foil.x[k]_dim).
 */
void GetTEInfoDeriv(const MatrixXd& foil_x, int k, int dim, double& dtcp,
                    double& dtdp, double& dhTE, double& ddtdx);

/**
 * @brief d(gamma)/d(foil.x): full [N+1 x 2N] Jacobian via reverse-mode AIC.
 *
 * @param state FoilState containing the current solution and geometry.
 * @param alpha Angle of attack in degrees.
 * @return MatrixXd of shape [N+1 x 2N].
 */
MatrixXd DGammaDFoilX(const FoilState& state, double alpha);

/**
 * @brief InviscidVelocity at xi with full forward sensitivity matrix.
 * Returns:
 *   V        : the velocity at xi
 *   V_xi     : d(V)/d(xi) [2x2]
 *   V_xfoil  : combined d(V)/d(foil.x_panels) at fixed xi PLUS V_G * dgam_dx
 *              [2 x 2N]
 */
struct VelDeriv {
    Vector2d V;
    Matrix2d V_xi;
    MatrixXd V_xfoil;
};

/**
 * @brief InviscidVelocity at xi with full forward sensitivity matrix.
 *
 * @param state FoilState containing the current solution and geometry.
 * @param xi Evaluation point in 2D space.
 * @param dgam_dx d(gamma)/d(foil.x) matrix of shape [N+1 x 2N].
 * @param alpha Angle of attack in degrees.
 * @return VelDeriv containing V, V_xi, and V_xfoil.
 */
VelDeriv InviscidVelocityWithDeriv(const FoilState& state, const Vector2d& xi,
                                   const MatrixXd& dgam_dx, double alpha);

/**
 * @brief Overload: use arbitrary vortex strengths `vortex` (length N+1, with
 * the streamfunction slot in entry N — but the velocity computation only uses
 * entries 0..N-1) and a corresponding dvortex_dx [N+1 x 2N] chain. Used by
 * Phase 2b to chain through V_for_v at wake stations (Cgam direct-derivative
 * chain). Pass dvortex_dx = MatrixXd::Zero(...) to hold the vortex strengths
 * fixed (no gamma chain).
 *
 * @param state FoilState containing the current solution and geometry.
 * @param xi Evaluation point in 2D space.
 * @param vortex Arbitrary vortex strengths vector of length N+1.
 * @param dvortex_dx Corresponding chain Jacobian of shape [N+1 x 2N].
 * @param alpha Angle of attack in degrees.
 * @return VelDeriv containing V, V_xi, and V_xfoil.
 */
VelDeriv InviscidVelocityWithDerivExt(const FoilState& state,
                                      const Vector2d& xi,
                                      const VectorXd& vortex,
                                      const MatrixXd& dvortex_dx, double alpha);

} /* namespace scirocco */
