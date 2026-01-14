// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_COUPLED_SOLVER_H
#define CPPFOIL_CORE_COUPLED_SOLVER_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

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

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_COUPLED_SOLVER_H
