// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_UTIL_LINEARIZATION_H
#define CPPFOIL_UTIL_LINEARIZATION_H

#include "cppfoil/config.h"

namespace cppfoil {

// Result containing value and its linearization (derivative)
// V = value type, L = linearization type
template<typename V, typename L>
struct LinResult {
    V value;
    L lin;  // Linearization w.r.t. state U
};

// Common type aliases for boundary layer property functions

// Scalar value with derivative w.r.t. 4-component state U
using ScalarLin4 = LinResult<double, Vector4d>;

// Scalar value with derivative w.r.t. 8-component state (U1, U2)
using ScalarLin8 = LinResult<double, Vector8d>;

// 3-vector residual with 3x8 Jacobian w.r.t. (U1, U2)
using Vec3Lin8 = LinResult<Vector3d, Matrix38d>;

// For properties that also depend on x (arc length)
template<typename V, typename LinU, typename LinX>
struct LinResultWithX {
    V value;
    LinU lin_U;  // Linearization w.r.t. U
    LinX lin_x;  // Linearization w.r.t. x
};

// Residual at a station: 3 residuals with derivatives w.r.t. U[8] and x[2]
using ResidualLin = LinResultWithX<Vector3d, Matrix38d, Matrix32d>;

}  // namespace cppfoil

#endif  // CPPFOIL_UTIL_LINEARIZATION_H
