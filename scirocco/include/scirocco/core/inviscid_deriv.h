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
 * @brief d(A*v)/d(foil.x) for a fixed vector v: [N+1 x 2N].
 *
 * Row i is the derivative of the i-th row of the AIC system applied to v,
 * with v held fixed. Evaluates each panel-influence derivative once, so the
 * cost is O(N^2). This is the kernel behind AICMatVecVJP,
 * AICResidualSensVJP and DGammaDFoilX.
 *
 * @param state FoilState with solved inviscid problem
 * @param v Vector [N+1] (gamma-like, streamfunction slot in entry N)
 * @return Jacobian [N+1 x 2N]; the Kutta row N is zero.
 */
MatrixXd AICMatVecJacobianX(const FoilState& state, const VectorXd& v);

/**
 * @brief d(A*gamma - rhs)/d(foil.x) at the current gamma: [N+1 x 2N].
 *
 * The inviscid residual's partial derivative with respect to the node
 * coordinates, at fixed gamma. d(gamma)/d(foil.x) = -A^{-1} times this.
 *
 * @param state FoilState with solved inviscid problem (AIC, gam)
 * @param alpha Angle of attack [degrees]
 * @return Jacobian [N+1 x 2N]
 */
MatrixXd AICResidualJacobianX(const FoilState& state, double alpha);

/**
 * @brief Compute the vector-Jacobian product: lambda^T * d(A*gamma -
 * rhs)/dx_foil
 *
 * Equal to AICResidualJacobianX(state, alpha)^T * lambda.
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
