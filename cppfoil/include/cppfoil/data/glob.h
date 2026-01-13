// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_GLOB_H
#define CPPFOIL_DATA_GLOB_H

#include "cppfoil/config.h"

namespace cppfoil {

// Global system variables
struct Glob {
    int Nsys = 0;             // Number of equations and states
    MatrixXd U;               // Primary states [th, ds, sa, ue] [4 x Nsys]
    MatrixXd dU;              // Primary state update [4 x Nsys]
    double dalpha = 0.0;      // Angle of attack update
    bool conv = true;         // Converged flag
    VectorXd R;               // Residuals [3*Nsys]
    SparseMatrixd R_U;        // Residual Jacobian w.r.t. primary states
    SparseMatrixd R_x;        // Residual Jacobian w.r.t. xi [3*Nsys x Nsys]
    SparseMatrixd R_V;        // Global Jacobian [4*Nsys x 4*Nsys]
    bool realloc = false;     // If true, system Jacobians will be reallocated
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_GLOB_H
