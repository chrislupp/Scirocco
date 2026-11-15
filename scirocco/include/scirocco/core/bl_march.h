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
*/

/**
 * @file bl_march.h
 * @brief Mixed-mode boundary-layer march run before each coupled Newton
 * step (XFOIL's MRCHDU).
 */
#pragma once

#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief Re-solve every boundary-layer station by a downstream march at the
 * current edge velocity, locating transition as it goes.
 *
 * Each station solves its three boundary-layer equations, with the upstream
 * station fixed, together with one constraint that holds the station near
 * the current iterate: a combination of ue and Hk taken quasi-normal to the
 * boundary layer's own ue-Hk characteristic. In attached flow this pins ue
 * (direct mode); approaching separation, where ue stops responding to Hk,
 * it pins Hk (inverse mode), so the Goldstein singularity is never met.
 *
 * On laminar stretches the transition test is re-applied at every local
 * iteration with AmplifyInterval, the test UpdateTransition uses, so the
 * march and UpdateTransition agree on where transition sits at a converged
 * state. Stations whose local solve fails are extrapolated from upstream.
 * The stagnation pair and the first wake station are left to the Newton
 * step.
 *
 * The residuals are the global system's, so a converged state is a fixed
 * point.
 * @param state Foil state with a coupled solution in progress
 * @return Number of stations whose local solve failed
 */
int BoundaryLayerMarch(FoilState& state);

}  // namespace scirocco
