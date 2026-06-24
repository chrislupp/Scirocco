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
#include "scirocco/scirocco.h"

#include <cmath>
#include <iostream>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/paneling.h"
#include "scirocco/core/wake.h"
#include "scirocco/data/state.h"

namespace scirocco {

Scirocco::Scirocco(const std::string& naca, int npanel) {
    NacaPoints(state_.geom, naca);
    MakePanels(state_.geom, npanel, state_.foil);
}

Scirocco::Scirocco(const MatrixXd& coords, int npanel) {
    SetCoords(state_.geom, coords);
    MakePanels(state_.geom, npanel, state_.foil);
}

void Scirocco::SetOper(double alpha, double cl, double Re, double Ma) {
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

    state_.oper.Vinf = 1.0;  /* Normalized */
    state_.oper.rho = 1.0;   /* Normalized density */
}

bool Scirocco::Solve() {
    if (state_.param.verb >= 1) {
        std::cout << "Scirocco::Solve viscous=" << state_.oper.viscous
                  << ", Re=" << state_.oper.Re << '\n';
    }
    if (state_.oper.viscous) {
        /* Initialize parameters for viscous solve.
         * Do NOT override rtol/niglob/verb here — those are settable defaults
         * (Param: rtol=1e-11, niglob=50, verb=0). Earlier code clamped rtol
         * to 1e-6 which is far too loose for adjoint precision. */
        state_.oper.initbl = true;

        /* Initialize thermodynamics */
        state_.param.Vinf = state_.oper.Vinf;
        state_.param.Minf = state_.oper.Ma;
        /* muinf = rho * V * chord / Re */
        state_.param.muinf = state_.oper.rho * state_.oper.Vinf *
                             state_.geom.chord / state_.oper.Re;
        state_.param.mu0 = state_.param.muinf;
        state_.param.rho0 = state_.oper.rho;
        /* Karman-Tsien compressibility correction
         *   KTb = sqrt(1 - Ma^2),  KTl = Ma^2 / (1 + sqrt(1 - Ma^2))^2
         * At Ma=0: KTb=1, KTl=0 (no correction). Earlier code hardcoded
         * KTl=0 unconditionally, killing Ma sensitivity. */
        if (state_.oper.Ma > 0.0) {
            double b = std::sqrt(1.0 - state_.oper.Ma * state_.oper.Ma);
            state_.param.KTb = b;
            state_.param.KTl =
                (state_.oper.Ma * state_.oper.Ma) / ((1.0 + b) * (1.0 + b));
        } else {
            state_.param.KTb = 1.0;
            state_.param.KTl = 0.0;
        }

        /* Viscous solver */
        SolveViscous(state_);
        return state_.glob.conv;
    } else {
        /* Inviscid solver */
        SolveInviscid(state_);
        CalcForce(state_);
        return true;
    }
}

void Scirocco::GeomFlap(const Vector2d& xzhinge, double eta) {
    scirocco::GeomFlap(state_, xzhinge, eta);
}

void Scirocco::GeomAddCamber(const MatrixXd& xzcamb) {
    scirocco::GeomAddCamber(state_, xzcamb);
}

void Scirocco::GeomDerotate() { scirocco::GeomDerotate(state_); }

Sensitivity Scirocco::ComputeSensitivities() {
    if (state_.oper.viscous && state_.glob.U.cols() > 0) {
        sensitivity_ = ViscousSensitivities(state_);
        /* Alpha/Mach/Re scalar viscous chains: Phase 5b/5c — deferred.
         * For inviscid-derivable scalar sensitivities, leave the existing
         * inviscid helpers' values (zeros if not invoked). */
    } else {
        sensitivity_ = InviscidSensitivities(state_);
        AlphaSensitivities(state_, sensitivity_);
        MachSensitivities(state_, sensitivity_);
    }
    return sensitivity_;
}

}  // namespace scirocco
