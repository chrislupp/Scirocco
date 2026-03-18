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
#include "cppfoil/cppfoil.h"

#include <iostream>
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/paneling.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/data/state.h"

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
    if (state_.param.verb >= 1) {
        std::cout << "Cppfoil::Solve viscous=" << state_.oper.viscous
                  << ", Re=" << state_.oper.Re << std::endl;
    }
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

void Cppfoil::GeomFlap(const Vector2d& xzhinge, double eta) {
    cppfoil::GeomFlap(state_, xzhinge, eta);
}

void Cppfoil::GeomAddCamber(const MatrixXd& xzcamb) {
    cppfoil::GeomAddCamber(state_, xzcamb);
}

void Cppfoil::GeomDerotate() {
    cppfoil::GeomDerotate(state_);
}

}  // namespace cppfoil
