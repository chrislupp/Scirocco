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
#ifndef CPPFOIL_DATA_POST_H
#define CPPFOIL_DATA_POST_H

#include "cppfoil/config.h"

namespace cppfoil {

// Post-processing outputs and distributions
struct Post {
    VectorXd cp;            // cp distribution
    VectorXd cpi;           // Inviscid cp distribution
    double cl = 0.0;        // Lift coefficient
    VectorXd cl_ue;         // Linearization of cl w.r.t. ue [N, airfoil only]
    double cl_alpha = 0.0;  // Linearization of cl w.r.t. alpha
    double cm = 0.0;        // Moment coefficient
    double cdpi = 0.0;      // Near-field pressure drag coefficient
    double cd = 0.0;        // Total drag coefficient
    Vector4d cd_U;          // Linearization of cd w.r.t. last wake state [4]
    double cdf = 0.0;       // Skin friction drag coefficient
    double cdp = 0.0;       // Pressure drag coefficient

    // Distributions
    VectorXd th;   // Momentum thickness distribution
    VectorXd ds;   // Displacement thickness distribution
    VectorXd sa;   // Amplification factor/shear lag coeff distribution
    VectorXd ue;   // Edge velocity (compressible) distribution
    VectorXd uei;  // Inviscid edge velocity (compressible) distribution
    VectorXd cf;   // Skin friction distribution
    VectorXd Ret;  // Re_theta distribution
    VectorXd Hk;   // Kinematic shape parameter distribution
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_POST_H
