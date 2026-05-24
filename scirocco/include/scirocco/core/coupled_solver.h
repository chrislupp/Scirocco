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
#ifndef SCIROCCO_CORE_COUPLED_SOLVER_H
#define SCIROCCO_CORE_COUPLED_SOLVER_H

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

// Solve viscous problem (top-level)
void SolveViscous(FoilState& state);

// Coupled Newton iteration loop
void SolveCoupled(FoilState& state);

// Build global residual system
void BuildGlobSys(FoilState& state);

// Solve the global linear system
void SolveGlob(FoilState& state);

// Update state with under-relaxation
void UpdateState(FoilState& state);

// Calculate ue sensitivities w.r.t. mass sources
void CalcUeM(FoilState& state);

// Rebuild ue sensitivities after stagnation point move
void RebuildUeM(FoilState& state);

// Rebuild inviscid solution after alpha change
void RebuildIsol(FoilState& state);

// Store SparseLU factorization for adjoint solves (post-convergence)
void StoreAdjointFactorization(FoilState& state);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_COUPLED_SOLVER_H
