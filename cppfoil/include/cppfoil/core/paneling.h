// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

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
