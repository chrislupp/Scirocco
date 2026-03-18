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
#ifndef CPPFOIL_DATA_GLOB_H
#define CPPFOIL_DATA_GLOB_H

#include "cppfoil/config.h"

namespace cppfoil {

// Global system variables
struct Glob {
    int Nsys = 0;          // Number of equations and states
    MatrixXd U;            // Primary states [th, ds, sa, ue] [4 x Nsys]
    MatrixXd dU;           // Primary state update [4 x Nsys]
    double dalpha = 0.0;   // Angle of attack update
    bool conv = true;      // Converged flag
    VectorXd R;            // Residuals [3*Nsys]
    SparseMatrixd R_U;     // Residual Jacobian w.r.t. primary states
    SparseMatrixd R_x;     // Residual Jacobian w.r.t. xi [3*Nsys x Nsys]
    SparseMatrixd R_V;     // Global Jacobian [4*Nsys x 4*Nsys]
    bool realloc = false;  // If true, system Jacobians will be reallocated
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_GLOB_H
