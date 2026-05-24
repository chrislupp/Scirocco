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
#ifndef SCIROCCO_CORE_PANEL_INFLUENCE_H
#define SCIROCCO_CORE_PANEL_INFLUENCE_H

#include "scirocco/config.h"

namespace scirocco {

// Panel geometry information
struct PanelInfo {
    double r1, r2;          // Distances from xi to panel endpoints
    double theta1, theta2;  // Angles from xi to panel endpoints
    double x1, x2;          // Local x coordinates (along panel)
    double y;               // Local y coordinate (perpendicular to panel)
    double len;             // Panel length
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
//   vdir : Optional velocity direction for projection (nullptr for both
//   components) onmid : True if xi is at panel midpoint (use limit values)
// OUTPUT:
//   a, b : Influence coefficients for left and right vortex strengths
void PanelLinVortexVelocity(const Matrix2d& Xj, const Vector2d& xi,
                            const Vector2d* vdir, bool onmid, Vector2d& a,
                            Vector2d& b);

// Linear vortex panel streamfunction influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   a, b : Influence coefficients for left and right vortex strengths
void PanelLinVortexStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b);

// Constant source panel velocity influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
//   vdir : Optional velocity direction for projection (nullptr for both
//   components)
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
void PanelLinSourceStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b);

// Linear source panel velocity influence
// INPUT:
//   Xj : Panel endpoint coordinates [2 x 2]
//   xi : Evaluation point
// OUTPUT:
//   a, b : Velocity influence vectors for left and right source strengths
void PanelLinSourceVelocity(const Matrix2d& Xj, const Vector2d& xi, Vector2d& a,
                            Vector2d& b);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_PANEL_INFLUENCE_H
