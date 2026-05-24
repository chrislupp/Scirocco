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
#ifndef SCIROCCO_CORE_WAKE_H
#define SCIROCCO_CORE_WAKE_H

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

// Recompute foil.s (arc-length parametrization) from foil.x as cumulative
// panel lengths. Call after perturbing foil.x so that foil.s is consistent.
void RecomputeArcLength(FoilState& state);

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

}  // namespace scirocco

#endif  // SCIROCCO_CORE_WAKE_H
