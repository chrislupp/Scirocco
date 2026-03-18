/*
    cppfoil: C++ Airfoil Analysis Library

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
#ifndef CPPFOIL_DATA_ISOL_H
#define CPPFOIL_DATA_ISOL_H

#include <array>
#include "cppfoil/config.h"

namespace cppfoil {

// Inviscid solution variables
struct Isol {
    MatrixXd AIC;              // Aerodynamic influence coefficient matrix
    MatrixXd gamref;           // 0,90-deg alpha vortex strengths at airfoil nodes [N x 2]
    VectorXd gam;              // Vortex strengths at airfoil nodes (current alpha)
    double sstag = 0.0;        // s location of stagnation point
    Vector2d sstag_g{0.0, 0.0};   // Linearization of sstag w.r.t. adjacent gammas
    Vector2d sstag_ue{0.0, 0.0};  // Linearization of sstag w.r.t. adjacent ue values
    std::array<int, 2> Istag{0, 0};  // Node indices before/after stagnation
    VectorXd sgnue;            // +/- 1 on upper/lower surface nodes
    VectorXd xi;               // Distance from stagnation point at all nodes
    VectorXd uewi;             // Inviscid edge velocity in wake
    MatrixXd uewiref;          // 0,90-deg alpha inviscid ue solutions on wake [Nw x 2]
    Vector2d xstag{0.0, 0.0};  // x,y location of stagnation point
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_ISOL_H
