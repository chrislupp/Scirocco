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

#include <algorithm>
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
namespace {
/// Newton iterations allowed to a warm-started viscous solve.
constexpr int kWarmIterations = 30;
void ValidateOperatingConditions(double alpha, double cl, double Re,
                                 double Ma) {
    if (!std::isfinite(alpha) || !std::isfinite(cl) || !std::isfinite(Re) ||
        !std::isfinite(Ma) || Ma < 0.0 || Ma >= 1.0)
        throw std::invalid_argument(
            "Operating conditions must be finite with 0 <= Ma < 1");
    // Finite negative cl/Re values retain their alpha-mode/inviscid meaning.
}

void ValidateResults(const FoilState& state) {
    if (!std::isfinite(state.post.cl) || !std::isfinite(state.post.cd) ||
        !std::isfinite(state.post.cm) || !std::isfinite(state.post.cdf) ||
        !std::isfinite(state.post.cdp) || !state.post.cp.allFinite() ||
        !state.isol.gam.allFinite() || !state.post.cf.allFinite() ||
        !state.post.th.allFinite() || !state.post.ds.allFinite() ||
        !state.wake.x.allFinite() ||
        (state.oper.viscous && !state.glob.U.allFinite()))
        throw std::runtime_error(
            "Scirocco solution contains nonfinite results");
}
}  // namespace

Scirocco::Scirocco(const std::string& naca, int npanel) {
    ValidatePanelCount(npanel);
    NacaPoints(state_.geom, naca);
    MakePanels(state_.geom, npanel, state_.foil);
}

Scirocco::Scirocco(const MatrixXd& coords, int npanel) {
    ValidatePanelCount(npanel);
    SetCoords(state_.geom, coords);
    MakePanels(state_.geom, npanel, state_.foil);
}

void Scirocco::SetOper(double alpha, double cl, double Re, double Ma) {
    sensitivity_valid_ = false;
    state_.glob.conv = false;
    state_.glob.coupled_factorization.reset();
    ValidateOperatingConditions(alpha, cl, Re, Ma);
    if (state_.oper.viscous && Re <= 0.0) ClearSolution(state_);
    state_.oper.Re = Re > 0.0 ? Re : 0.0;
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

    state_.oper.Vinf = 1.0; /* Normalized */
    state_.oper.rho = 1.0;  /* Normalized density */
}

bool Scirocco::Solve() {
    sensitivity_valid_ = false;
    state_.glob.conv = false;
    state_.glob.coupled_factorization.reset();
    try {
        ValidateOperatingConditions(state_.oper.alpha, state_.oper.cltgt,
                                    state_.oper.Re, state_.oper.Ma);
        ValidateCoordinates(state_.foil.x);
        if (state_.foil.N != state_.foil.x.cols())
            throw std::runtime_error("Panel count does not match coordinates");
        if (state_.param.verb >= 1) {
            std::cout << "Scirocco::Solve viscous=" << state_.oper.viscous
                      << ", Re=" << state_.oper.Re << '\n';
        }
        if (state_.oper.viscous) {
            /* Initialize parameters for viscous solve.
             * Do NOT override rtol/niglob/verb here — those are settable
             * defaults (Param: rtol=1e-11, niglob=50, verb=0). Earlier code
             * clamped rtol to 1e-6 which is far too loose for adjoint
             * precision. */
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

            /* Viscous solver: warm start from the previous converged
             * boundary layer, then cold start if that does not converge.
             * The warm attempt skips continuation; the cold one does not. */
            bool solved = false;
            if (state_.param.warm_start && warm_U_.cols() > 0) {
                FoilState cold = state_;
                state_.oper.initbl = false;
                state_.glob.U = warm_U_;
                state_.glob.Nsys = static_cast<int>(warm_U_.cols());
                state_.vsol.turb = warm_turb_;
                state_.isol.Istag = warm_istag_;
                const bool continuation = state_.param.use_continuation;
                const int niglob = state_.param.niglob;
                state_.param.use_continuation = false;
                /* A good initial guess converges in a few iterations; give
                 * up early on a bad one (e.g. the stagnation point crossed a
                 * node) rather than spend the full budget. */
                state_.param.niglob = std::min(niglob, kWarmIterations);
                try {
                    SolveViscous(state_);
                    solved = state_.glob.conv;
                } catch (const std::exception&) {
                    solved = false;
                }
                state_.param.use_continuation = continuation;
                state_.param.niglob = niglob;
                if (!solved) state_ = std::move(cold);
            }
            if (!solved) SolveViscous(state_);
            if (state_.glob.conv) {
                ValidateResults(state_);
                warm_U_ = state_.glob.U;
                warm_turb_ = state_.vsol.turb;
                warm_istag_ = state_.isol.Istag;
            }
            return state_.glob.conv;
        } else {
            /* Inviscid solver */
            SolveInviscid(state_);
            CalcForce(state_);
            ValidateResults(state_);
            if (state_.oper.givencl) {
                state_.glob.conv = false;
                for (int iteration = 0; iteration < state_.param.niglob;
                     ++iteration) {
                    const double residual = state_.post.cl - state_.oper.cltgt;
                    if (std::abs(residual) < state_.param.rtol) {
                        state_.glob.conv = true;
                        return true;
                    }
                    CalcForceLin(state_);
                    Sensitivity slope;
                    AlphaSensitivities(state_, slope);
                    if (!std::isfinite(slope.dcl_dalpha) ||
                        std::abs(slope.dcl_dalpha) < 1e-10)
                        return false;
                    state_.oper.alpha -=
                        std::clamp(residual / slope.dcl_dalpha, -2.0, 2.0);
                    BuildGamma(state_, state_.oper.alpha);
                    CalcForce(state_);
                    ValidateResults(state_);
                }
                return false;
            }
            return true;
        }
    } catch (...) {
        state_.glob.conv = false;
        state_.glob.coupled_factorization.reset();
        throw;
    }
}

void Scirocco::GeomFlap(const Vector2d& xzhinge, double eta) {
    sensitivity_valid_ = false;
    scirocco::GeomFlap(state_, xzhinge, eta);
}

void Scirocco::GeomAddCamber(const MatrixXd& xzcamb) {
    sensitivity_valid_ = false;
    scirocco::GeomAddCamber(state_, xzcamb);
}

void Scirocco::GeomDerotate() {
    sensitivity_valid_ = false;
    scirocco::GeomDerotate(state_);
}

Sensitivity Scirocco::ComputeSensitivities() {
    sensitivity_valid_ = false;
    if (!state_.glob.conv) {
        throw std::runtime_error(
            "ComputeSensitivities: solve must converge first");
    }
    if (state_.oper.viscous && state_.glob.U.cols() > 0) {
        sensitivity_ = ViscousSensitivities(state_);
    } else {
        sensitivity_ = InviscidSensitivities(state_);
    }
    sensitivity_valid_ = true;
    return sensitivity_;
}

void Scirocco::CheckSensitivityAvailable() const {
    if (!sensitivity_valid_ || !state_.glob.conv)
        throw std::runtime_error(
            "ComputeSensitivities must succeed at the current state before "
            "reading gradients");
}

}  // namespace scirocco
