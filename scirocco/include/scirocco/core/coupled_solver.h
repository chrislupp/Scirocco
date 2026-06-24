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
#ifndef SCIROCCO_CORE_COUPLED_SOLVER_H
#define SCIROCCO_CORE_COUPLED_SOLVER_H

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief Solve viscous problem (top-level)
 * @param state Foil state to solve in place.
 */
void SolveViscous(FoilState& state);

/**
 * @brief Coupled Newton iteration loop
 * @param state Foil state iterated to convergence.
 */
void SolveCoupled(FoilState& state);

/**
 * @brief Build global residual system
 * @param state Foil state whose global residual is assembled.
 */
void BuildGlobSys(FoilState& state);

/**
 * @brief Solve the global linear system
 * @param state Foil state whose global linear system is solved.
 */
void SolveGlob(FoilState& state);

/**
 * @brief Update state with under-relaxation
 * @param state Foil state updated in place.
 */
void UpdateState(FoilState& state);

/**
 * @brief Calculate ue sensitivities w.r.t. mass sources
 * @param state Foil state in which ue sensitivities are computed.
 */
void CalcUeM(FoilState& state);

/**
 * @brief Rebuild ue sensitivities after stagnation point move
 * @param state Foil state in which ue sensitivities are rebuilt.
 */
void RebuildUeM(FoilState& state);

/**
 * @brief Rebuild inviscid solution after alpha change
 * @param state Foil state in which the inviscid solution is rebuilt.
 */
void RebuildIsol(FoilState& state);

/**
 * @brief Store SparseLU factorization for adjoint solves (post-convergence)
 * @param state Foil state whose R_V factorization is stored.
 */
void StoreAdjointFactorization(FoilState& state);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_COUPLED_SOLVER_H
