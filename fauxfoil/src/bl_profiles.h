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

/// @file bl_profiles.h
/// Boundary-layer velocity profiles drawn on the airfoil (BL, BLC).

#include <vector>

#include "plots.h"

namespace fauxfoil {

/// A velocity profile at one point on the airfoil surface.
struct BlProfile {
    Point wall;        ///< Surface point (airfoil coordinates)
    double ds = 0.0;   ///< Displacement thickness
    double ret = 0.0;  ///< Re_theta
    double hk = 0.0;   ///< Kinematic shape parameter
    bool turbulent = false;
    /// Profile in airfoil coordinates: the wall normal up to the boundary-
    /// layer edge, then u/Qinf * weight across it.
    std::vector<Point> normal;
    std::vector<Point> curve;
};

/// u/Ue at height eta = y / delta through a boundary layer of shape
/// factor @p hk: the Pohlhausen quartic family when laminar, a power law
/// when turbulent. Also returns delta / delta*.
double ProfileU(double eta, double hk, bool turbulent,
                double* delta_over_ds = nullptr);

/// Profile at arc length @p s along the airfoil, measured from the upper
/// trailing edge round the leading edge to the lower trailing edge.
/// @param weight Plot length of u = Qinf, in chords.
/// @return False where there is no viscous solution.
bool ProfileAtArc(const Solution& solution, double s, double weight,
                  BlProfile& profile);

/// Profile at the surface point nearest @p point (airfoil coordinates).
bool ProfileNear(const Solution& solution, const Point& point, double weight,
                 BlProfile& profile);

/// Total arc length of the airfoil contour.
double ArcLength(const Solution& solution);

/// Draw @p profile on @p scene through the plot's airfoil mapping.
void DrawProfile(Scene& scene, const BlProfile& profile);

}  // namespace fauxfoil
