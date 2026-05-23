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
#ifndef CPPFOIL_CORE_BOUNDARY_LAYER_H
#define CPPFOIL_CORE_BOUNDARY_LAYER_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"
#include "cppfoil/util/linearization.h"

namespace cppfoil {

// Initialize boundary layer from inviscid edge velocity
void InitBoundaryLayer(FoilState& state);

// Compute residual at a station (non-transition)
ResidualLin ResidualStation(const Param& param, const Vector2d& x,
                            const Matrix42d& U, const Vector2d& Aux);

// Compute residual at a transition station (combined laminar + turbulent)
// U1 should be laminar, U2 should be turbulent
// Returns combined residual and Jacobian
ResidualLin ResidualTransition(FoilState& state, const Param& param_base,
                               const Vector2d& x, const Matrix42d& U,
                               const Vector2d& Aux);

// Thwaites correlation for stagnation initialization
// Returns (theta, delta*) at stagnation point
std::pair<double, double> ThwaitesInit(double K, double nu);

// Build parameter structure for a surface
Param BuildParam(const FoilState& state, int si);

// Update parameter for a specific station
void StationParam(const FoilState& state, Param& param, int i);

// Result for stagnation state extrapolation
struct StagnationResult {
    Vector4d Ust;     // Stagnation state [th, ds, sa, ue]
    Matrix48d Ust_U;  // Linearization w.r.t. U1 and U2 (4x8)
    Matrix42d Ust_x;  // Linearization w.r.t. x1 and x2 (4x2)
    double xst;       // Stagnation point location (small)
};

// Extrapolate two states to stagnation point
StagnationResult StagnationState(const Matrix42d& U, const Vector2d& x);

// Result for wake system
struct WakeResult {
    Vector3d R;      // Residual [th, ds, sa]
    Matrix312d R_U;  // Linearization w.r.t. lower, upper, wake states (3x12)
    std::array<int, 3> J;  // Indices: lower TE, upper TE, first wake
    // ∂R/∂Ma · Ma, non-zero only when at least one TE side is laminar
    // (because ctw is then built from GetCttr which is Ma-dependent).
    Vector3d lin_Ma_times_Ma{Vector3d::Zero()};
};

// Construct residual system for wake initialization
WakeResult WakeSys(const FoilState& state, Param& param);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_BOUNDARY_LAYER_H
