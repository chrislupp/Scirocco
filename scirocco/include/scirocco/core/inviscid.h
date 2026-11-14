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
 * @brief Build and solve the inviscid linear system for gamma
 * @param state FoilState with panels set up; filled with AIC matrix and
 * vorticity solutions
 * @param alpha Angle of attack [degrees]
 */
void BuildGamma(FoilState& state, double alpha);

/**
 * @brief Compute inviscid tangential velocity at every node
 * @param state FoilState with inviscid solution
 * @return Inviscid edge velocity at airfoil and wake points
 */
VectorXd GetUeInv(const FoilState& state);

/**
 * @brief Solve AIC x = rhs, reusing the factorization stored by BuildGamma
 * when it belongs to the current AIC.
 * @param state FoilState with an influence matrix
 * @param rhs Right-hand sides [N+1 x k]
 * @return AIC^-1 rhs
 */
MatrixXd SolveAIC(const FoilState& state, const MatrixXd& rhs);

/**
 * @brief Compute 0,90-deg reference inviscid tangential velocities
 * @param state FoilState with inviscid solution
 * @return 0,90-deg inviscid tangential velocity at all points [N+Nw x 2]
 */
MatrixXd GetUeInvRef(const FoilState& state);

/**
 * @brief Compute inviscid velocity at an arbitrary point
 * @param X Panel node coordinates [2 x N]
 * @param G Vortex strengths at each node [N]
 * @param Vinf Freestream speed magnitude
 * @param alpha Angle of attack [degrees]
 * @param xi Point at which velocity is desired
 * @param compute_lin True to also compute linearization
 * @param V_G (optional) Linearization d(V)/d(G) [2 x N]
 * @return Velocity vector [2]
 */
Vector2d InviscidVelocity(const MatrixXd& X, const VectorXd& G, double Vinf,
                          double alpha, const Vector2d& xi,
                          bool compute_lin = false, MatrixXd* V_G = nullptr);

/**
 * @brief Solve inviscid problem (top-level)
 * @param state FoilState with geometry and operating conditions; updated with
 * inviscid solution
 */
void SolveInviscid(FoilState& state);

/**
 * @brief Initialize thermodynamic parameters
 * @param state FoilState to initialize
 */
void InitThermo(FoilState& state);

/**
 * @brief Clear all solution data
 * @param state FoilState to clear
 */
void ClearSolution(FoilState& state);

}  // namespace scirocco
