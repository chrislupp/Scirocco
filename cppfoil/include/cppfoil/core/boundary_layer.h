// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

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

// Thwaites correlation for stagnation initialization
// Returns (theta, delta*) at stagnation point
std::pair<double, double> ThwaitesInit(double K, double nu);

// Build parameter structure for a surface
Param BuildParam(const FoilState& state, int si);

// Update parameter for a specific station
void StationParam(const FoilState& state, Param& param, int i);

// Result for stagnation state extrapolation
struct StagnationResult {
    Vector4d Ust;        // Stagnation state [th, ds, sa, ue]
    Matrix48d Ust_U;     // Linearization w.r.t. U1 and U2 (4x8)
    Matrix42d Ust_x;     // Linearization w.r.t. x1 and x2 (4x2)
    double xst;          // Stagnation point location (small)
};

// Extrapolate two states to stagnation point
StagnationResult StagnationState(const Matrix42d& U, const Vector2d& x);

// Result for wake system
struct WakeResult {
    Vector3d R;          // Residual [th, ds, sa]
    Matrix312d R_U;      // Linearization w.r.t. lower, upper, wake states (3x12)
    std::array<int, 3> J; // Indices: lower TE, upper TE, first wake
};

// Construct residual system for wake initialization
WakeResult WakeSys(const FoilState& state, Param& param);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_BOUNDARY_LAYER_H
