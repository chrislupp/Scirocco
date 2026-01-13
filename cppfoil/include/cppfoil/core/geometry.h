// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_GEOMETRY_H
#define CPPFOIL_CORE_GEOMETRY_H

#include <string>
#include "cppfoil/config.h"
#include "cppfoil/data/geom.h"

namespace cppfoil {

// Trailing edge information result
struct TEInfo {
    Vector2d t;      // Unit tangent at TE (average of upper/lower)
    double hTE;      // TE gap height
    double dtdx;     // d(theta)/dx at TE
    double tcp;      // TE constant source panel strength factor
    double tdp;      // TE vortex panel strength factor
};

// Generate NACA 4 or 5-digit airfoil points
// INPUT:
//   digits : NACA digit string (e.g., "0012", "23012")
// OUTPUT:
//   geom : Geom structure with point coordinates filled in
void NacaPoints(Geom& geom, const std::string& digits);

// Set coordinates from user-provided points
// INPUT:
//   geom : Geom structure to fill
//   coords : Point coordinates [2 x N], ordered from TE lower to TE upper
void SetCoords(Geom& geom, const MatrixXd& coords);

// Get trailing edge information from panel coordinates
// INPUT:
//   x : Panel node coordinates [2 x N]
// OUTPUT:
//   TEInfo structure with gap, tangent, and strength factors
TEInfo GetTEInfo(const MatrixXd& x);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_GEOMETRY_H
