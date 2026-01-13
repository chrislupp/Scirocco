// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_ISOL_H
#define CPPFOIL_DATA_ISOL_H

#include <array>
#include "cppfoil/config.h"

namespace cppfoil {

// Inviscid solution variables
struct Isol {
    MatrixXd AIC;              // Aerodynamic influence coefficient matrix
    MatrixXd gamref;           // 0,90-deg alpha vortex strengths at airfoil nodes [N x 2]
    VectorXd gam;              // Vortex strengths at airfoil nodes (current alpha)
    double sstag = 0.0;        // s location of stagnation point
    Vector2d sstag_g{0.0, 0.0};   // Linearization of sstag w.r.t. adjacent gammas
    Vector2d sstag_ue{0.0, 0.0};  // Linearization of sstag w.r.t. adjacent ue values
    std::array<int, 2> Istag{0, 0};  // Node indices before/after stagnation
    VectorXd sgnue;            // +/- 1 on upper/lower surface nodes
    VectorXd xi;               // Distance from stagnation point at all nodes
    VectorXd uewi;             // Inviscid edge velocity in wake
    MatrixXd uewiref;          // 0,90-deg alpha inviscid ue solutions on wake [Nw x 2]
    Vector2d xstag{0.0, 0.0};  // x,y location of stagnation point
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_ISOL_H
