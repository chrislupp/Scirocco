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
#ifndef CPPFOIL_CORE_PANELING_H
#define CPPFOIL_CORE_PANELING_H

#include "cppfoil/config.h"
#include "cppfoil/data/geom.h"
#include "cppfoil/data/panel.h"
#include "cppfoil/core/spline.h"

namespace cppfoil {

// Result from curvature-based spacing
struct SpacingResult {
    MatrixXd X;       // New point coordinates [2 x N]
    VectorXd S;       // Arc-length values at points
    MatrixXd XS;      // Tangent vectors [2 x N]
};

// Create panels from geometry using curvature-based spacing
// INPUT:
//   geom : Geometry data with point coordinates
//   npanel : Number of panels (N-1 panels means N nodes)
//   stgt : Optional target s-values (empty for automatic spacing)
// OUTPUT:
//   panel : Panel structure with nodes, arc-lengths, and tangents
void MakePanels(const Geom& geom, int npanel, Panel& panel,
                const VectorXd& stgt = VectorXd());

// Spline points with curvature-based spacing
// INPUT:
//   points : Input points to spline [2 x N]
//   N : Number of output points
//   Ufac : Uniformity factor (1 = normal; >1 = more uniform)
//   TEfac : TE resolution factor (1 = normal; >1 = higher resolution at TE)
//   stgt : Optional target s-values (empty for automatic)
// OUTPUT:
//   SpacingResult with new points, arc-lengths, and tangents
SpacingResult SplineCurvature(const MatrixXd& points, int N,
                               double Ufac, double TEfac,
                               const VectorXd& stgt = VectorXd());

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_PANELING_H
