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

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_BOUNDARY_LAYER_H
