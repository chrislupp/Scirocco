// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_PANEL_INFLUENCE_H
#define CPPFOIL_CORE_PANEL_INFLUENCE_H

#include "cppfoil/config.h"

namespace cppfoil {

// Panel geometry information
struct PanelInfo {
    double r1, r2;       // Distances from xi to panel endpoints
    double theta1, theta2;  // Angles from xi to panel endpoints
    double x1, x2;       // Local x coordinates (along panel)
    double y;            // Local y coordinate (perpendicular to panel)
    double len;          // Panel length
};

// Compute panel geometry information
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2] (j, j+1)
//   xi : Point at which to evaluate
// OUTPUT:
//   PanelInfo structure with geometric quantities
PanelInfo GetPanelInfo(const Matrix2d& Xj, const Vector2d& xi);

// Linear vortex panel velocity influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
//   vdir : Optional velocity direction for projection (nullptr for both components)
//   onmid : True if xi is at panel midpoint (use limit values)
// OUTPUT:
//   a, b : Influence coefficients for left and right vortex strengths
void PanelLinVortexVelocity(const Matrix2d& Xj, const Vector2d& xi,
                             const Vector2d* vdir, bool onmid,
                             Vector2d& a, Vector2d& b);

// Linear vortex panel streamfunction influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   a, b : Influence coefficients for left and right vortex strengths
void PanelLinVortexStream(const Matrix2d& Xj, const Vector2d& xi,
                           double& a, double& b);

// Constant source panel velocity influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
//   vdir : Optional velocity direction for projection (nullptr for both components)
// OUTPUT:
//   Velocity influence (or scalar if vdir provided)
Vector2d PanelConstSourceVelocity(const Matrix2d& Xj, const Vector2d& xi,
                                   const Vector2d* vdir);

// Constant source panel streamfunction influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   Streamfunction influence coefficient
double PanelConstSourceStream(const Matrix2d& Xj, const Vector2d& xi);

// Linear source panel streamfunction influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   a, b : Influence coefficients for left and right source strengths
void PanelLinSourceStream(const Matrix2d& Xj, const Vector2d& xi,
                           double& a, double& b);

// Linear source panel velocity influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   a, b : Velocity influence vectors for left and right source strengths
void PanelLinSourceVelocity(const Matrix2d& Xj, const Vector2d& xi,
                             Vector2d& a, Vector2d& b);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_PANEL_INFLUENCE_H
