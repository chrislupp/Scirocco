// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_WAKE_H
#define CPPFOIL_CORE_WAKE_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// Build wake grid from inviscid solution
// INPUT:
//   state : FoilState with inviscid solution
// OUTPUT:
//   state.wake : Filled with wake panel data
//   state.isol.uewi, uewiref : Wake edge velocities
void BuildWake(FoilState& state);

// Find stagnation point on airfoil surface
// INPUT:
//   state : FoilState with inviscid solution
// OUTPUT:
//   state.isol.sstag, Istag, xstag : Stagnation point information
//   state.isol.sgnue : Sign of edge velocity
//   state.isol.xi : Distance from stagnation point (foil only)
void StagPointFind(FoilState& state);

// Extend xi array to include wake points
// Must be called after BuildWake
void ExtendXiToWake(FoilState& state);

// Move stagnation point during viscous iterations
// INPUT/OUTPUT:
//   state : FoilState with viscous solution
void StagPointMove(FoilState& state);

// Identify surfaces (lower, upper, wake)
// INPUT/OUTPUT:
//   state : FoilState to update with surface indices
void IdentifySurfaces(FoilState& state);

// Set wake gap for blunt trailing edge
// INPUT/OUTPUT:
//   state : FoilState to update
void SetWakeGap(FoilState& state);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_WAKE_H
