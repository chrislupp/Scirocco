// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_PANEL_H
#define CPPFOIL_DATA_PANEL_H

#include "cppfoil/config.h"

namespace cppfoil {

// Panel data (used for both airfoil and wake)
struct Panel {
    int N = 0;      // Number of nodes
    MatrixXd x;     // Node coordinates [2 x N]
    VectorXd s;     // Arclength values at nodes
    MatrixXd t;     // Tangent vectors dx/ds, dy/ds [2 x N]
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_PANEL_H
