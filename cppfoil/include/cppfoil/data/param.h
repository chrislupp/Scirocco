// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_DATA_PARAM_H
#define CPPFOIL_DATA_PARAM_H

namespace cppfoil {

// Parameters
struct Param {
    // Solver parameters
    int verb = 1;             // Printing verbosity level (higher = more verbose)
    double rtol = 1e-10;      // Residual tolerance for Newton
    int niglob = 50;          // Maximum number of global iterations

    // Viscous parameters
    double ncrit = 9.0;       // Critical amplification factor
    double Cuq = 1.0;         // Scales the uq term in the shear lag equation
    double Dlr = 0.9;         // Wall/wake dissipation length ratio
    double SlagK = 5.6;       // Shear lag constant

    // Initial Ctau after transition
    double CtauC = 1.8;       // Ctau constant
    double CtauE = 3.3;       // Ctau exponent

    // G-beta locus: G = GA*sqrt(1+GB*beta) + GC/(H*Rt*sqrt(cf/2))
    double GA = 6.7;          // G-beta A constant
    double GB = 0.75;         // G-beta B constant
    double GC = 18.0;         // G-beta C constant

    // Operating conditions and thermodynamics
    double Minf = 0.0;        // Freestream Mach number
    double Vinf = 0.0;        // Freestream speed
    double muinf = 0.0;       // Freestream dynamic viscosity
    double mu0 = 0.0;         // Stagnation dynamic viscosity
    double rho0 = 1.0;        // Stagnation density
    double H0 = 0.0;          // Stagnation enthalpy
    double Tsrat = 0.35;      // Sutherland Ts/Tref
    double gam = 1.4;         // Ratio of specific heats
    double KTb = 1.0;         // Karman-Tsien beta = sqrt(1-Minf^2)
    double KTl = 0.0;         // Karman-Tsien lambda = Minf^2/(1+KTb)^2
    double cps = 0.0;         // Sonic cp

    // Station information
    bool simi = false;        // True at a similarity station
    bool turb = false;        // True at a turbulent station
    bool wake = false;        // True at a wake station
};

}  // namespace cppfoil

#endif  // CPPFOIL_DATA_PARAM_H
