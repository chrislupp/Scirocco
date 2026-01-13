// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_VSOL_H
#define CPPFOIL_DATA_VSOL_H

#include <array>
#include <vector>
#include "cppfoil/config.h"

namespace cppfoil {

// Viscous solution variables
struct Vsol {
    VectorXd th;              // Momentum thickness [Nsys]
    VectorXd ds;              // Displacement thickness [Nsys]
    std::array<std::vector<int>, 3> Is;  // Surface indices [lower, upper, wake]
    VectorXd wgap;            // Wake gap over wake points
    MatrixXd ue_m;            // Linearization of ue w.r.t. mass (all nodes)
    SparseMatrixd sigma_m;    // d(source)/d(mass) matrix
    MatrixXd ue_sigma;        // d(ue)/d(source) matrix
    std::vector<bool> turb;   // Flag indicating turbulent (true) or laminar (false)
    double xt = 0.0;          // Transition location (xi) on current surface
    Matrix2d Xt = Matrix2d::Zero();  // Transition xi/x for lower and upper surfaces [2x2]
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_VSOL_H
