// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_OPER_H
#define CPPFOIL_DATA_OPER_H

namespace cppfoil {

// Operating conditions and flags
struct Oper {
    double Vinf = 1.0;       // Velocity magnitude
    double alpha = 0.0;      // Angle of attack [degrees]
    double rho = 1.0;        // Density
    double cltgt = 0.0;      // Lift coefficient target
    bool givencl = false;    // True if cl is given instead of alpha
    bool initbl = true;      // True to initialize the boundary layer
    bool viscous = false;    // True for viscous analysis
    bool redowake = false;   // True to rebuild wake after alpha changes
    double Re = 1e5;         // Reynolds number
    double Ma = 0.0;         // Mach number
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_OPER_H
