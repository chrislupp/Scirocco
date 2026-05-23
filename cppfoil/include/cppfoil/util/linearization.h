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
#ifndef CPPFOIL_UTIL_LINEARIZATION_H
#define CPPFOIL_UTIL_LINEARIZATION_H

#include "cppfoil/config.h"

namespace cppfoil {

// Result containing value and its linearization (derivative)
// V = value type, L = linearization type
template <typename V, typename L>
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
template <typename V, typename LinU, typename LinX>
struct LinResultWithX {
    V value;
    LinU lin_U;  // Linearization w.r.t. U
    LinX lin_x;  // Linearization w.r.t. x
};

// Residual at a station: 3 residuals with derivatives w.r.t. U[8] and x[2],
// plus optional ∂R/∂Re·Re for the analytic Re chain. lin_Ret stores
// sum_k (∂R/∂Ret_k · Ret_k) where k ranges over every distinct Ret value
// used in the residual (Ret_1 at U1, Ret_2 at U2, Ret_m at Um=0.5(U1+U2),
// plus any others). Since Re only enters Ret as Ret(U, mu0) ∝ Re/mu0 with
// ∂Ret/∂Re = Ret/Re, the chain ∂R/∂Re = lin_Ret / Re. Pre-multiplying by
// Ret here keeps callers from needing to track per-evaluation Ret values.
//
// Default zero — ResidualStation populates when populate_lin_Ret is true.
struct ResidualLin {
    Vector3d value;
    Matrix38d lin_U;          // ∂R/∂[U1, U2]
    Matrix32d lin_x;          // ∂R/∂[xi1, xi2]
    Vector3d lin_Ret = Vector3d::Zero();  // sum_k ∂R/∂Ret_k · Ret_k
};

}  // namespace cppfoil

#endif  // CPPFOIL_UTIL_LINEARIZATION_H
