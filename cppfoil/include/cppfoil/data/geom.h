// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_GEOM_H
#define CPPFOIL_DATA_GEOM_H

#include <string>
#include "cppfoil/config.h"

namespace cppfoil {

// Geometry data
struct Geom {
    double chord = 1.0;           // Chord length
    double wakelen = 1.0;         // Wake extent length, in chords
    int npoint = 0;               // Number of geometry representation points
    std::string name = "noname";  // Airfoil name, e.g. "NACA 0012"
    MatrixXd xpoint;              // Point coordinates [2 x npoint]
    Vector2d xref{0.25, 0.0};     // Moment reference point
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_GEOM_H
