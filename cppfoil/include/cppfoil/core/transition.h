// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_TRANSITION_H
#define CPPFOIL_CORE_TRANSITION_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// March amplification factor along a surface
// Returns: index of last laminar station in Is array (before transition)
int MarchAmplification(FoilState& state, int si);

// Update transition locations based on current solution
void UpdateTransition(FoilState& state);

// Store transition location for a surface
void StoreTransition(FoilState& state, int si, int i);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_TRANSITION_H
