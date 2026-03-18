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
#ifndef CPPFOIL_CORE_INVISCID_H
#define CPPFOIL_CORE_INVISCID_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// Build and solve the inviscid linear system for gamma
// INPUT:
//   state : FoilState with panels set up
//   alpha : Angle of attack [degrees]
// OUTPUT:
//   state.isol : Filled with AIC matrix and vorticity solutions
void BuildGamma(FoilState& state, double alpha);

// Compute inviscid tangential velocity at every node
// INPUT:
//   state : FoilState with inviscid solution
// OUTPUT:
//   Inviscid edge velocity at airfoil and wake points
VectorXd GetUeInv(const FoilState& state);

// Compute 0,90-deg reference inviscid tangential velocities
// INPUT:
//   state : FoilState with inviscid solution
// OUTPUT:
//   0,90-deg inviscid tangential velocity at all points [N+Nw x 2]
MatrixXd GetUeInvRef(const FoilState& state);

// Compute inviscid velocity at an arbitrary point
// INPUT:
//   X : Panel node coordinates [2 x N]
//   G : Vortex strengths at each node [N]
//   Vinf : Freestream speed magnitude
//   alpha : Angle of attack [degrees]
//   xi : Point at which velocity is desired
//   compute_lin : True to also compute linearization
// OUTPUT:
//   V : Velocity vector [2]
//   V_G : (optional) Linearization d(V)/d(G) [2 x N]
Vector2d InviscidVelocity(const MatrixXd& X, const VectorXd& G,
                           double Vinf, double alpha, const Vector2d& xi,
                           bool compute_lin = false,
                           MatrixXd* V_G = nullptr);

// Solve inviscid problem (top-level)
// INPUT:
//   state : FoilState with geometry and operating conditions
// OUTPUT:
//   state : Updated with inviscid solution
void SolveInviscid(FoilState& state);

// Initialize thermodynamic parameters
// INPUT/OUTPUT:
//   state : FoilState to initialize
void InitThermo(FoilState& state);

// Clear all solution data
// INPUT/OUTPUT:
//   state : FoilState to clear
void ClearSolution(FoilState& state);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_INVISCID_H
