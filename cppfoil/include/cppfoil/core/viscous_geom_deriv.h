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
*/
#ifndef CPPFOIL_CORE_VISCOUS_GEOM_DERIV_H
#define CPPFOIL_CORE_VISCOUS_GEOM_DERIV_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// ============================================================================
// Phase 1 geometric chain functions for the viscous adjoint.
//
// These provide d(geom_quantity)/d(foil.x) for the geometric quantities that
// the viscous residual R_V depends on (xi, wake.x, wake.s, wgap, sstag).
//
// All functions take dfoil_s_dx [N x 2N] (derivative of airfoil arc-length
// foil.s w.r.t. foil.x in row-major (k, dim)=(k*2+dim) column ordering) as
// an input parameter so the test can supply mfoil's spline-based version
// (matching the complex-step reference) while production can use the
// cumulative-panel-norm version from DFoilSDFoilXCumulative.
//
// Column layout for output matrices: column index = 2*k + dim, where
//   k   ∈ [0, N)  is the foil node index
//   dim ∈ {0, 1} is x (0) or y (1).
// ============================================================================

// Compute d(foil.s)/d(foil.x) using cumulative panel norms (matches cppfoil's
// RecomputeArcLength). Returns matrix of shape [N x 2N].
//
// NOTE: this does NOT match mfoil's spline-based foil.s (~4e-5 relative
// discrepancy near stagnation). Tests against mfoil_cs should load
// d(foil.s)/d(foil.x) from the HDF5 reference instead.
MatrixXd DFoilSDFoilXCumulative(const Panel& foil);

// d(sstag)/d(foil.x) — stagnation point arc-length, gamma-based interpolation.
// Returns VectorXd of length 2N.
//
// Chain pieces:
//   1. d(sstag)/d(foil.s[Istag]) via interpolation weights (w1, w2)
//   2. d(sstag)/d(gamma[Istag]) via sstag_g (already cached in state.isol)
//      → reverse-mode adjoint of AIC: λ_i = (AIC^T)^{-1} e_{Istag[i]}
//      → d(gamma_Istag[i])/d(foil.x) = -AICResidualSensVJP(λ_i, alpha)
VectorXd DSstagDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

// d(wake.s)/d(foil.x) — cumulative wake arc-length. Returns [Nw x 2N].
//
// wake.s[i] = foil.s[N-1] + sv[i]
//   sv[0] = 0
//   sv[i] depends on ds1 = 0.5*((foil.s[1]-foil.s[0]) + (foil.s[N-1]-foil.s[N-2]))
//   and on chord (constant in foil.x → contributes 0).
//
// Chain:
//   d(wake.s[i])/d(foil.x) = d(foil.s[N-1])/d(foil.x) + d(sv[i])/d(ds1) * d(ds1)/d(foil.s) * d(foil.s)/d(foil.x)
//
// This is the chain the previous attempt missed (per memory:
// "Missing d(wake.s)/d(foil.x) chain causes ~0.216 adjoint error at TE nodes").
MatrixXd DWakeSDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

// d(xi)/d(foil.x) — combined for ALL Nsys stations (foil + wake).
// Returns [Nsys x 2N] where Nsys = N + Nw.
//
// Foil-side (i ∈ [0, N)):
//   xi[i] = |foil.s[i] - sstag|
//   d(xi[i])/d(foil.x) = sign(foil.s[i] - sstag) *
//                        (d(foil.s[i])/d(foil.x) - d(sstag)/d(foil.x))
//
// Wake-side (i ∈ [N, N+Nw)):
//   xi[N+j] = wake.s[j] - sstag
//   d(xi[N+j])/d(foil.x) = d(wake.s[j])/d(foil.x) - d(sstag)/d(foil.x)
MatrixXd DXiDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

// d(wgap)/d(foil.x) — wake gap from TE blunt-base extrapolation.
// Returns [Nw x 2N].
//
// wgap[i] depends on:
//   - hTE = |TE gap|  → depends on foil.x[0, N-1] only
//   - dtdx = TE thickness slope → depends on foil.x[{0,1,N-2,N-1}]
//   - tcp, tdp factors → depends on same TE nodes
//   - xi[N+i] - xi[N] (wake-only arc-length from TE) → depends on
//     wake.s[i] - foil.s[N-1] via cancellation with sstag
//
// The cubic formula is:
//   wgap[i] = hTE * (1 + (2 + flen*dtdx)*xib) * (1-xib)^2,  xib = (xi-xi0)/Lw
//   Lw = flen*hTE
MatrixXd DWgapDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

// d(wake.x)/d(foil.x) — wake node coordinates from predictor-corrector
// streamline tracing. Returns [2*Nw x 2N], where row 2*j+d = wake.x(d, j).
//
// PHASE 1 STATUS: STUB (returns zeros).
//
// Implementation requires a forward sweep through BuildWake's predictor-
// corrector loop. Each xyw[:,i+1] depends on:
//   - xyw[:,i] (recursive)
//   - inviscid_velocity at xyw[:,i] and at predictor point xyw[:,i] + ds*v1
//   - foil.x, foil.gam, alpha
//
// This is the chain that historically failed at TE nodes (~0.216 error).
// It must be implemented carefully with a per-iteration forward sensitivity
// matrix dXw_i = d(wake.x[:,i])/d(foil.x) ∈ R^{2 x 2N}, then propagated.
//
// DEFERRED to follow-on work after the closed-form chains validate.
MatrixXd DWakeXDFoilX(const FoilState& state);

// d(TE info)/d(foil.x[k]_dim): expose the static helper from adjoint.cpp.
// Only nodes k ∈ {0, 1, N-2, N-1} produce nonzero output.
void GetTEInfoDeriv(const MatrixXd& foil_x, int k, int dim, double& dtcp,
                    double& dtdp, double& dhTE, double& ddtdx);

// d(gamma)/d(foil.x): full [N+1 x 2N] Jacobian via reverse-mode AIC.
MatrixXd DGammaDFoilX(const FoilState& state, double alpha);

// InviscidVelocity at xi with full forward sensitivity matrix.
// Returns:
//   V        : the velocity at xi
//   V_xi     : d(V)/d(xi) [2x2]
//   V_xfoil  : combined d(V)/d(foil.x_panels) at fixed xi PLUS V_G * dgam_dx
//              [2 x 2N]
struct VelDeriv {
    Vector2d V;
    Matrix2d V_xi;
    MatrixXd V_xfoil;
};
VelDeriv InviscidVelocityWithDeriv(const FoilState& state, const Vector2d& xi,
                                   const MatrixXd& dgam_dx, double alpha);

// Overload: use arbitrary vortex strengths `vortex` (length N+1, with the
// streamfunction slot in entry N — but the velocity computation only uses
// entries 0..N-1) and a corresponding dvortex_dx [N+1 x 2N] chain. Used by
// Phase 2b to chain through V_for_v at wake stations (Cgam direct-derivative
// chain). Pass dvortex_dx = MatrixXd::Zero(...) to hold the vortex strengths
// fixed (no gamma chain).
VelDeriv InviscidVelocityWithDerivExt(const FoilState& state,
                                      const Vector2d& xi,
                                      const VectorXd& vortex,
                                      const MatrixXd& dvortex_dx,
                                      double alpha);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_VISCOUS_GEOM_DERIV_H
