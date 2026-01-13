// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_FORCES_H
#define CPPFOIL_CORE_FORCES_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// Calculate force and moment coefficients
// INPUT:
//   state : FoilState with solution (inviscid or viscous)
// OUTPUT:
//   state.post : Filled with cl, cd, cm, cp, etc.
void CalcForce(FoilState& state);

// Compute distributions (theta, delta*, cf, etc.)
// INPUT:
//   state : FoilState with viscous solution
// OUTPUT:
//   state.post : Distribution quantities calculated
void GetDistributions(FoilState& state);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_FORCES_H
