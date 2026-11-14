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
 * @file coupled_march.h
 * @brief Quasi-simultaneous boundary-layer march used to smooth coupled
 * Newton iterates.
 */
#pragma once

#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief One downstream sweep that re-solves every boundary-layer station
 * with its upstream neighbour fixed (nonlinear Gauss-Seidel on the coupled
 * system), locating transition as it goes.
 *
 * At each station the three boundary-layer equations are solved together
 * with that node's edge-velocity equation,
 *     ue_i = ueinv_i + sum_j ue_m(i, j) m_j,
 * holding the mass defect m_j of every other node at its current value
 * (upstream nodes already swept). Keeping the diagonal interaction term makes
 * the local problem regular through separation (quasi-simultaneous
 * interaction). On laminar stretches the amplification equation decides
 * where transition falls, so the transition point always agrees with the
 * boundary layer it sits in. This is the role of MRCHDU in XFOIL, which runs
 * before every global Newton step.
 *
 * The residuals are the global system's, so a converged state is a fixed
 * point. Stations whose local solve fails keep their previous values. The
 * stagnation pair and the first wake station are left to the global step.
 * @param state Foil state with a coupled solution in progress
 * @return Number of stations whose local solve failed
 */
int CoupledMarch(FoilState& state);

}  // namespace scirocco
