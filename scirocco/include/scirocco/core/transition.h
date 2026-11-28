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
#include "scirocco/data/state.h"

namespace scirocco {

/** Where the trip (Param::xtrip_lower, xtrip_upper) sits on one surface. */
struct Trip {
    int station = -1;     /**< Index i in Is[si]: the trip is between stations
                               i-1 and i. Negative if the surface has no trip */
    double frac = 1.0;    /**< Fraction of that interval:
                               xi_trip = xi(i-1) + frac * (xi(i) - xi(i-1)) */
    bool at_node = false; /**< The trip was moved to a node (frac 0 or 1), so
                               it does not move with the geometry */
};

/**
 * @brief Locate the trip on a surface, as XFOIL's XIFSET does.
 *
 * The trip is at x = xtrip * chord on the geometric upper or lower surface,
 * interpolated linearly in x between panel nodes. It is absent if xtrip is 1
 * or more, if it lies at or behind the trailing-edge node, or if the
 * stagnation point has moved past it (XFOIL then also leaves that side
 * free). A trip ahead of the leading edge acts at the leading-edge node, and
 * one between the first two stations after the stagnation point, which stay
 * laminar, acts at the second.
 * @param state FoilState with surfaces (vsol.Is) and xi set
 * @param si Surface index (0=lower, 1=upper)
 */
Trip FindTrip(const FoilState& state, int si);

/**
 * @brief March amplification factor along a surface
 * @param state FoilState with current solution
 * @param si Surface index
 * @return Index of last laminar station in Is array (before transition):
 * the station before e^N transition or before the trip, whichever is first
 */
int MarchAmplification(FoilState& state, int si);

/**
 * @brief Amplification at node @p ii from node @p ip by the laminar
 * amplification equation, holding theta, delta* and ue at their current
 * values. This is the transition test of MarchAmplification and of the
 * boundary-layer march, so both place transition at the same interval.
 * @param state Foil state
 * @param ip Upstream node
 * @param ii Downstream node
 * @param converged Set to false if the scalar Newton did not converge
 * @return Amplification factor at @p ii
 */
double AmplifyInterval(const FoilState& state, int ip, int ii,
                       bool* converged = nullptr);

/**
 * @brief Update transition locations based on current solution
 * @param state FoilState with current solution
 */
void UpdateTransition(FoilState& state);

/**
 * @brief Store transition location for a surface
 * @param state FoilState with current solution
 * @param si Surface index (0=lower, 1=upper)
 * @param i Station index in Is array (transition is between i-1 and i)
 */
void StoreTransition(FoilState& state, int si, int i);

}  // namespace scirocco
