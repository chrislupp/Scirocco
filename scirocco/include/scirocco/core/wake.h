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

/** @brief Recompute foil.s (arc-length parametrization) from foil.x as
 * cumulative panel lengths. Call after perturbing foil.x so that foil.s is
 * consistent.
 *  @param state FoilState whose foil.x coordinates have been updated
 */
void RecomputeArcLength(FoilState& state);

/** @brief Build wake grid from inviscid solution
 *  @param state FoilState with inviscid solution
 *
 *  On output:
 *  - state.wake is filled with wake panel data
 *  - state.isol.uewi, uewiref contain wake edge velocities
 */
void BuildWake(FoilState& state);

/** @brief Find stagnation point on airfoil surface
 *  @param state FoilState with inviscid solution
 *
 *  On output:
 *  - state.isol.sstag, Istag, xstag contain stagnation point information
 *  - state.isol.sgnue contains the sign of edge velocity
 *  - state.isol.xi contains the distance from stagnation point (foil only)
 */
void StagPointFind(FoilState& state);

/** @brief Extend xi array to include wake points
 *
 *  Must be called after BuildWake.
 *  @param state FoilState with wake already built
 */
void ExtendXiToWake(FoilState& state);

/** @brief Move stagnation point during viscous iterations
 *  @param state FoilState with viscous solution (modified in place)
 */
void StagPointMove(FoilState& state);

/** @brief Identify surfaces (lower, upper, wake)
 *  @param state FoilState to update with surface indices
 */
void IdentifySurfaces(FoilState& state);

/** @brief Set wake gap for blunt trailing edge
 *  @param state FoilState to update
 */
void SetWakeGap(FoilState& state);

} /* namespace scirocco */
