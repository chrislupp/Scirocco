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
#ifndef SCIROCCO_DATA_OPER_H
#define SCIROCCO_DATA_OPER_H

namespace scirocco {

// Operating conditions and flags
struct Oper {
    double Vinf = 1.0;      // Velocity magnitude
    double alpha = 0.0;     // Angle of attack [degrees]
    double rho = 1.0;       // Density
    double cltgt = 0.0;     // Lift coefficient target
    bool givencl = false;   // True if cl is given instead of alpha
    bool initbl = true;     // True to initialize the boundary layer
    bool viscous = false;   // True for viscous analysis
    bool redowake = false;  // True to rebuild wake after alpha changes
    double Re = 1e5;        // Reynolds number
    double Ma = 0.0;        // Mach number
};

}  // namespace scirocco

#endif  // SCIROCCO_DATA_OPER_H
