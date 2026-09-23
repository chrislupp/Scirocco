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
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief Compute the vector-Jacobian product: lambda^T * d(A*gamma -
 * rhs)/dx_foil
 *
 * This avoids forming the full [N+1 x N+1 x 2N] tensor.
 *
 * @param state FoilState with solved inviscid problem (AIC, gam, gamref)
 * @param lambda Adjoint vector [N+1] (solution of A^T * lambda = dJ/dgamma)
 * @param alpha Angle of attack [degrees]
 * @return VJP vector [2*N]: the total sensitivity contribution from the AIC
 * system
 */
VectorXd AICResidualSensVJP(const FoilState& state, const VectorXd& lambda,
                            double alpha);

/**
 * @brief Compute lambda^T * d(A * v)/dx_foil
 *
 * Same panel-chain machinery as AICResidualSensVJP but with the gamma vector
 * replaced by an arbitrary v (length N+1, with the streamfunction slot in
 * entry N) and WITHOUT the rhs derivative term. Used by Phase 2b's
 * ComputeRueVJP for the dBp chain where v = [Bp * phi; 0].
 *
 * @param state FoilState with solved inviscid problem
 * @param lambda Adjoint vector [N+1]
 * @param v Arbitrary vector [N+1] replacing gamma
 * @return VJP vector [2*N]
 */
VectorXd AICMatVecVJP(const FoilState& state, const VectorXd& lambda,
                      const VectorXd& v);

}  // namespace scirocco
