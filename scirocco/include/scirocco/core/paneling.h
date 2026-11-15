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
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#pragma once

#include "scirocco/config.h"
#include "scirocco/core/spline.h"
#include "scirocco/data/geom.h"
#include "scirocco/data/panel.h"

namespace scirocco {

/**
 * @brief Result from curvature-based spacing
 */
struct SpacingResult {
    MatrixXd X;  /* New point coordinates [2 x N] */
    VectorXd S;  /* Arc-length values at points */
    MatrixXd XS; /* Tangent vectors [2 x N] */
};

/**
 * @brief Create panels from geometry using curvature-based spacing
 * @param geom Geometry data with point coordinates
 * @param npanel Number of panels (N-1 panels means N nodes)
 * @param panel Panel structure with nodes, arc-lengths, and tangents
 * @param stgt Optional target s-values (empty for automatic spacing)
 */
void MakePanels(const Geom& geom, int npanel, Panel& panel,
                const VectorXd& stgt = VectorXd());

/**
 * @brief Refresh panel.s (arc length) and panel.t (tangent) from panel.x,
 * leaving panel.x untouched. Use after directly setting panel.x to new node
 * coordinates (e.g. from an upstream geometry discipline) so derived
 * geometric state stays consistent with panel.x — without re-spacing,
 * which would break upstream gradient chains.
 * @param panel Panel structure to refresh
 */
void RefreshPanelGeometry(Panel& panel);

/**
 * @brief Spline points with curvature-based spacing
 * @param points Input points to spline [2 x N]
 * @param N Number of output points
 * @param Ufac Uniformity factor (1 = normal; >1 = more uniform)
 * @param TEfac TE resolution factor (1 = normal; >1 = higher resolution at TE)
 * @param stgt Optional target s-values (empty for automatic)
 * @return SpacingResult with new points, arc-lengths, and tangents
 */
SpacingResult SplineCurvature(const MatrixXd& points, int N, double Ufac,
                              double TEfac, const VectorXd& stgt = VectorXd());

}  // namespace scirocco
