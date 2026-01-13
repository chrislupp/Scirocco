// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/cppfoil.h"

#include <iostream>
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/paneling.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"

namespace cppfoil {

Cppfoil::Cppfoil(const std::string& naca, int npanel) {
    NacaPoints(state_.geom, naca);
    MakePanels(state_.geom, npanel, state_.foil);
}

Cppfoil::Cppfoil(const MatrixXd& coords, int npanel) {
    SetCoords(state_.geom, coords);
    MakePanels(state_.geom, npanel, state_.foil);
}

void Cppfoil::SetOper(double alpha, double cl, double Re, double Ma) {
    state_.oper.alpha = alpha;
    state_.oper.Ma = Ma;

    if (cl >= 0) {
        state_.oper.cltgt = cl;
        state_.oper.givencl = true;
    } else {
        state_.oper.givencl = false;
    }

    if (Re > 0) {
        state_.oper.Re = Re;
        state_.oper.viscous = true;
    } else {
        state_.oper.viscous = false;
    }

    state_.oper.Vinf = 1.0;  // Normalized
    state_.oper.rho = 1.0;   // Normalized density
}

bool Cppfoil::Solve() {
    std::cout << "Cppfoil::Solve viscous=" << state_.oper.viscous
              << ", Re=" << state_.oper.Re << std::endl;
    if (state_.oper.viscous) {
        // Initialize parameters for viscous solve
        state_.oper.initbl = true;
        state_.param.verb = 1;  // Normal verbose output
        state_.param.niglob = 25;  // Max iterations
        state_.param.rtol = 1e-6;  // Convergence tolerance

        // Initialize thermodynamics
        state_.param.Vinf = state_.oper.Vinf;
        state_.param.Minf = state_.oper.Ma;
        // muinf = rho * V * chord / Re
        state_.param.muinf = state_.oper.rho * state_.oper.Vinf *
                             state_.geom.chord / state_.oper.Re;
        state_.param.mu0 = state_.param.muinf;  // Same for incompressible
        state_.param.rho0 = state_.oper.rho;
        state_.param.KTb = 1.0;  // Incompressible
        state_.param.KTl = 0.0;

        // Viscous solver
        SolveViscous(state_);
        return state_.glob.conv;
    } else {
        // Inviscid solver
        SolveInviscid(state_);
        CalcForce(state_);
        return true;
    }
}

}  // namespace cppfoil
