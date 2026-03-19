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
#ifndef CPPFOIL_DATA_SENSITIVITY_H
#define CPPFOIL_DATA_SENSITIVITY_H

#include "cppfoil/config.h"

namespace cppfoil {

// Sensitivity of aerodynamic outputs w.r.t. panel node coordinates
// Each matrix is [2 x N] matching the layout of foil.x
struct Sensitivity {
    MatrixXd dcl_dx;   // d(cl)/d(foil.x) [2 x N]
    MatrixXd dcd_dx;   // d(cd)/d(foil.x) [2 x N]
    MatrixXd dcm_dx;   // d(cm)/d(foil.x) [2 x N]
    MatrixXd dcdpi_dx; // d(cdpi)/d(foil.x) [2 x N] (inviscid pressure drag)
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_SENSITIVITY_H
