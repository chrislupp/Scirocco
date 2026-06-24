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
#ifndef SCIROCCO_CORE_TRANSITION_H
#define SCIROCCO_CORE_TRANSITION_H

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief March amplification factor along a surface
 * @param state FoilState with current solution
 * @param si Surface index
 * @return Index of last laminar station in Is array (before transition)
 */
int MarchAmplification(FoilState& state, int si);

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

#endif  // SCIROCCO_CORE_TRANSITION_H
