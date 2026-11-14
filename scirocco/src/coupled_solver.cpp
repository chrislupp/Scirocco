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
#include "scirocco/core/coupled_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

#include "scirocco/core/bl_march.h"
#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/captured_transition.h"
#include "scirocco/core/coupled_linear.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/panel_influence.h"
#include "scirocco/core/transition.h"
#include "scirocco/core/wake.h"
#include "scirocco/properties/bl_properties.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

namespace {
VectorXd CoupledResidual(const FoilState& state) {
    const int n = state.glob.Nsys;
    VectorXd residual(4 * n + (state.oper.givencl ? 1 : 0));
    residual.head(3 * n) = state.glob.R;
    const VectorXd ue = state.glob.U.row(3).transpose();
    residual.segment(3 * n, n) =
        ue - GetUeInv(state) -
        state.vsol.ue_m * state.glob.U.row(1).transpose().cwiseProduct(ue);
    if (state.oper.givencl) residual(4 * n) = state.post.cl - state.oper.cltgt;
    return residual;
}

/**
 * @brief Re-index a boundary layer from a previous solution to a new
 * stagnation point: station k from the stagnation point on each surface
 * takes the values of station k of the old solution (XFOIL's STMOVE).
 * Stations past the end of the old surface repeat its last station.
 * @param old Stagnation node pair the boundary layer was solved with
 */
void ShiftToStagnation(FoilState& state, const std::array<int, 2>& old) {
    const auto& now = state.isol.Istag;
    if (now == old) return;
    const int N = state.foil.N;
    const MatrixXd U = state.glob.U;
    const auto turb = state.vsol.turb;
    auto copy = [&](int dst, int src) {
        state.glob.U.col(dst) = U.col(src);
        state.vsol.turb[dst] = turb[src];
    };
    for (int k = 0; k <= now[0]; ++k) copy(now[0] - k, std::max(old[0] - k, 0));
    for (int k = 0; now[1] + k < N; ++k)
        copy(now[1] + k, std::min(old[1] + k, N - 1));
}

/**
 * @brief Locate transition from the amplification march (UpdateTransition).
 * @return true if no station changed between laminar and turbulent
 */
bool TransitionSettled(FoilState& state) {
    const auto turb = state.vsol.turb;
    UpdateTransition(state);
    return state.vsol.turb == turb;
}

VectorXd PackedState(const FoilState& state) {
    const int n = state.glob.Nsys;
    VectorXd value(4 * n + (state.oper.givencl ? 1 : 0));
    value.head(4 * n) = Eigen::Map<const VectorXd>(state.glob.U.data(), 4 * n);
    if (state.oper.givencl) value(4 * n) = state.oper.alpha;
    return value;
}

bool ValidTrialState(const FoilState& state) {
    if (!state.glob.U.allFinite() || !std::isfinite(state.oper.alpha))
        return false;
    if ((state.glob.U.topRows(2).array() <= 0.0).any()) return false;
    const auto& bracket = state.isol.Istag;
    if (state.glob.U(3, bracket[0]) == 0.0 ||
        state.glob.U(3, bracket[1]) == 0.0)
        return false;
    int left = bracket[0], right = bracket[1];
    while (left >= 0 && state.glob.U(3, left) <= 0.0) --left;
    while (right < state.foil.N && state.glob.U(3, right) <= 0.0) ++right;
    return left > 0 && right < state.foil.N;
}
}  // namespace

/**
 * @brief Top-level viscous solution
 * @param state Foil state
 */
void SolveViscous(FoilState& state) {
    if (state.param.verb >= 1) {
        std::cout << "SolveViscous: starting" << '\n';
    }

    /* Stagnation nodes of a supplied initial boundary layer (warm start) */
    const std::array<int, 2> entry_istag = state.isol.Istag;

    /* Initialize thermodynamic parameters */
    state.oper.viscous = true;
    InitThermo(state);

    /* Build inviscid solution first */
    BuildGamma(state, state.oper.alpha);
    StagPointFind(state);

    if (state.param.verb >= 1) {
        std::cout << "  Inviscid done, N=" << state.foil.N
                  << ", Istag=" << state.isol.Istag[0] << ","
                  << state.isol.Istag[1] << '\n';
    }

    /* Build wake */
    BuildWake(state);

    /* Extend xi to include wake points */
    ExtendXiToWake(state);

    if (state.param.verb >= 1) {
        std::cout << "  Wake done, Nw=" << state.wake.N << '\n';
    }

    /* Identify surfaces */
    IdentifySurfaces(state);

    /* Set wake gap */
    SetWakeGap(state);

    if (state.param.verb >= 1) {
        std::cout << "  Surfaces: lower=" << state.vsol.Is[0].size()
                  << ", upper=" << state.vsol.Is[1].size()
                  << ", wake=" << state.vsol.Is[2].size() << '\n';
    }

    /* Calculate ue sensitivities */
    CalcUeM(state);

    /* Initialize boundary layer, unless the caller supplied an initial
     * state of the right size (warm start). */
    const int nsys = state.foil.N + state.wake.N;
    if (!state.oper.initbl && state.glob.U.cols() == nsys &&
        static_cast<int>(state.vsol.turb.size()) == nsys) {
        state.glob.Nsys = nsys;
        state.glob.realloc = true;
        ShiftToStagnation(state, entry_istag);
    } else {
        InitBoundaryLayer(state);
        if (state.param.verb >= 1) {
            std::cout << "  BL init done, Nsys=" << state.glob.Nsys << '\n';
        }
    }

    /* Retain the original initial state for a clean continuation restart. */
    state.glob.continuation_attempts = 0;
    state.glob.continuation_steps = 0;
    state.glob.continuation_iterations = 0;
    std::unique_ptr<FoilState> initial;
    if (state.param.use_continuation)
        initial = std::make_unique<FoilState>(state);

    /* Solve coupled system */
    SolveCoupled(state);
    if (!state.glob.conv && initial) SolveContinuation(*initial, state);

    if (state.glob.conv &&
        state.param.transition == TransitionModel::Captured) {
        /* Converge the captured system from the cut-cell solution. Its
         * adjoint is not implemented yet, so no factorization is stored. */
        SolveCaptured(state);
    } else if (state.glob.conv) {
        /* Stay at the ue-based xi from Newton convergence — DO NOT swap to
         * gamma-based xi here. The previous code did StagPointFind after
         * Newton to "match mfoil's CS reference", but mfoil_cs's actual
         * sensitivity computation (via complex-step on alpha/foil.x) uses
         * mfoil's natural ue-based xi from solve_viscous; only the
         * extract_comprehensive_reference.py script does an extra
         * stagpoint_find for Phase 2 chain-partial validation. Swapping
         * here makes the adjoint's R_V (ue-based) inconsistent with the
         * VJP-side R_x (was gamma-based) and with sstag_ue (still ue-based),
         * and that inconsistency manifested as a ~0.1 abs error in R_V at
         * the ue-Istag column → 1e-5 rel floor on all viscous gradients.
         *
         * Phase 2 (residual_deriv_test) now explicitly applies the gamma-xi
         * swap in its test fixture; production sensitivities use ue-based. */
        /* The adjoint factorization is built on demand by
         * ViscousSensitivities / SolveViscousAdjoint. */
    }

    /* Calculate forces */
    CalcForce(state);
    GetDistributions(state);

    if (state.param.verb >= 1) {
        std::cout << "SolveViscous: done, cl=" << state.post.cl
                  << ", cd=" << state.post.cd << '\n';
    }
}

/**
 * @brief Coupled Newton iteration with pseudo-transient continuation (PTC).
 * SolveGlob reads state.glob.tau and adds (1/tau)*I to the BL rows of
 * the augmented Jacobian; SER ramps tau as the residual decreases so
 * we recover full Newton once we're close enough to converge cleanly.
 * @param state Foil state
 */
void SolveCoupled(FoilState& state) {
    state.glob.conv = false;
    state.glob.coupled_factorization.reset();
    int niglob = state.param.niglob;
    double rtol = state.param.rtol;

    state.glob.tau = state.param.tau_init;
    state.glob.R_norm_prev = 0.0;
    state.glob.trust_radius = state.param.trust_radius_initial;
    state.glob.trust_accepted = 0;
    state.glob.trust_rejected = 0;
    state.glob.coupled_iterations = 0;

    /* As in XFOIL, march the boundary layer before every Newton step,
     * starting with the initial guess. */
    const bool march = state.param.bl_march;
    if (march) BoundaryLayerMarch(state);
    for (int iter = 0; iter < niglob; ++iter) {
        ++state.glob.coupled_iterations;
        BuildGlobSys(state);
        // Python mfoil updates force values and the lift constraint each
        // iteration.
        CalcForce(state);
        double rnorm = CoupledResidual(state).norm();
        if (state.param.verb >= 1) {
            std::cout << "  Iteration " << iter << ": |R| = " << rnorm
                      << ", tau = " << state.glob.tau << '\n';
        }

        /* Near the round-off floor (above rtol on some cases, e.g. a
         * stagnation point on a node) the residual stops contracting.
         * Quadratic convergence would cut it by orders of magnitude, so a
         * step that fails to halve it below rtol_floor means the floor. */
        const bool at_floor = iter > 0 && rnorm < state.param.rtol_floor &&
                              rnorm > 0.5 * state.glob.R_norm_prev;
        if (rnorm < rtol || at_floor) {
            /* With the boundary-layer march placing transition, accept
             * only if the amplification march agrees with where it sits. */
            if (!march || TransitionSettled(state)) {
                state.glob.conv = true;
                return;
            }
            continue;
        }

        if (!std::isfinite(rnorm)) {
            state.glob.conv = false;
            if (state.param.verb >= 1) {
                std::cout << "  Warning: NaN in residual, stopping" << '\n';
            }
            return;
        }

        /* SER tau update — skip on first iter (no previous norm yet).
         * ratio = ||R_{n-1}|| / ||R_n|| : > 1 when residual decreased
         * (grow tau, less damping); < 1 when residual increased (shrink
         * tau, more damping). Clamp the ratio per-iter to avoid wild
         * swings and clamp tau itself to [tau_min, tau_max]. */
        if (state.param.use_ptc && iter > 0 && state.glob.R_norm_prev > 0.0) {
            double ratio = state.glob.R_norm_prev / rnorm;
            ratio = std::clamp(ratio, state.param.tau_shrink_clamp,
                               state.param.tau_grow_clamp);
            state.glob.tau =
                std::clamp(state.glob.tau * ratio, state.param.tau_min,
                           state.param.tau_max);
        }
        state.glob.R_norm_prev = rnorm;

        /* Solve and update */
        try {
            SolveGlob(state);
        } catch (const std::runtime_error&) {
            state.glob.conv = false;
            return;
        }
        if (state.param.use_trust_region) {
            if (!TakeTrustRegionStep(state)) return;
        } else {
            UpdateState(state);
        }

        // xi and its Jacobian must describe the current U even near
        // convergence. Freezing stagnation data below a residual threshold
        // changes the root while the adjoint still differentiates the
        // moving-stagnation model.
        if (!state.param.use_trust_region ||
            !state.param.trust_update_topology) {
            StagPointMove(state);
            if (march)
                BoundaryLayerMarch(state);
            else
                UpdateTransition(state);
        }
    }

    state.glob.conv = false;
    if (state.param.verb >= 1) {
        std::cout << "  Warning: viscous solver did not converge" << '\n';
    }
}

bool SolveContinuation(const FoilState& target, FoilState& result) {
    // The target wake is deliberately retained throughout the homotopy.
    // Rebuilding it at intermediate angles would change the final model
    // relative to a cold solve at the requested operating conditions.
    if (target.oper.redowake || target.param.continuation_max_attempts <= 0 ||
        target.param.continuation_max_iterations <= 0)
        return false;
    int attempts = 0, steps = 0, iterations = 0;
    auto solve = [&](FoilState& trial) {
        const int remaining =
            target.param.continuation_max_iterations - iterations;
        if (attempts >= target.param.continuation_max_attempts ||
            remaining <= 0)
            return false;
        trial.param.niglob = std::min(target.param.niglob, remaining);
        ++attempts;
        SolveCoupled(trial);
        iterations += trial.glob.coupled_iterations;
        trial.param.niglob = target.param.niglob;
        return trial.glob.conv;
    };
    auto record = [&]() {
        result.glob.continuation_attempts = attempts;
        result.glob.continuation_steps = steps;
        result.glob.continuation_iterations = iterations;
    };
    try {
        FoilState seed = target;
        seed.glob.coupled_factorization.reset();
        seed.oper.givencl = false;
        seed.oper.alpha = std::abs(target.oper.alpha) < 0.5 ? 2.0 : 0.0;
        RebuildIsol(seed);
        StagPointFind(seed);
        ExtendXiToWake(seed);
        IdentifySurfaces(seed);
        SetWakeGap(seed);
        CalcUeM(seed);
        seed.vsol.turb.clear();
        InitBoundaryLayer(seed);
        const FoilState seed_initial = seed;
        if (!solve(seed)) {
            seed = seed_initial;
            seed.param.use_ptc = !seed.param.use_ptc;
            if (!solve(seed)) {
                record();
                return false;
            }
        }
        const double start =
            target.oper.givencl ? seed.post.cl : seed.oper.alpha;
        const double finish =
            target.oper.givencl ? target.oper.cltgt : target.oper.alpha;
        const double span = finish - start;
        const double increment = target.oper.givencl ? 0.2 : 2.0;
        double progress = 0.0;
        double step =
            std::min(1.0, increment / std::max(std::abs(span), increment));
        while (attempts < target.param.continuation_max_attempts &&
               iterations < target.param.continuation_max_iterations) {
            const double next = std::min(1.0, progress + step);
            FoilState trial = seed;
            trial.oper.givencl = target.oper.givencl;
            const double value = next == 1.0 ? finish : start + next * span;
            if (target.oper.givencl)
                trial.oper.cltgt = value;
            else
                trial.oper.alpha = value;
            RebuildIsol(trial);
            if (target.param.verb >= 1)
                std::cout << "  Continuation: fraction=" << next
                          << ", target=" << value << '\n';
            if (solve(trial)) {
                ++steps;
                if (next == 1.0) {
                    trial.param = target.param;
                    trial.oper.cltgt = target.oper.cltgt;
                    result = std::move(trial);
                    record();
                    return true;
                }
                seed = std::move(trial);
                progress = next;
                step = std::min(1.0 - progress, 1.5 * step);
            } else {
                step *= 0.5;
                if (step < 1e-4) break;
            }
        }
    } catch (const std::exception& error) {
        // The caller retains the failed direct solve at the requested point,
        // never an intermediate solution at a different operating condition.
        if (target.param.verb >= 1)
            std::cout << "  Continuation failed: " << error.what() << '\n';
        record();
        return false;
    }
    record();
    return false;
}

/**
 * @brief Build global residual and Jacobian
 * @param state Foil state
 */
void BuildGlobSys(FoilState& state) {
    int Nsys = state.glob.Nsys;

    /* Allocate residual and Jacobians if needed */
    bool alloc = state.glob.realloc || state.glob.R.size() != 3 * Nsys;
    if (alloc) {
        state.glob.R.resize(3 * Nsys);
        state.glob.R_U.resize(3 * Nsys, 4 * Nsys);
        state.glob.R_x.resize(3 * Nsys, Nsys);
        state.glob.realloc = false;
    }

    state.glob.R.setZero();
    /* Allocate/zero the ∂R/∂Re·Re and ∂R/∂Ma·Ma accumulators (BL rows only). */
    if (state.glob.dR_dRe_times_Re.size() != 3 * Nsys) {
        state.glob.dR_dRe_times_Re = VectorXd::Zero(3 * Nsys);
    } else {
        state.glob.dR_dRe_times_Re.setZero();
    }
    if (state.glob.dR_dMa_times_Ma.size() != 3 * Nsys) {
        state.glob.dR_dMa_times_Ma = VectorXd::Zero(3 * Nsys);
    } else {
        state.glob.dR_dMa_times_Ma.setZero();
    }

    /* Build using triplets for sparse matrix */
    std::vector<Triplet> triplets_U, triplets_x;
    triplets_U.reserve(Nsys * 24); /* Rough estimate */
    triplets_x.reserve(Nsys * 6);

    /* Loop over surfaces (lower, upper, wake) */
    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        if (Is.empty()) continue;

        int N = static_cast<int>(Is.size());

        /* Get parameter structure for this surface */
        Param param = BuildParam(state, si);

        /* Get xi values for this surface */
        VectorXd xi(N);
        for (int idx = 0; idx < N; ++idx) {
            xi(idx) = state.isol.xi(Is[idx]);
        }

        /* Get auxiliary data (wake gap) */
        VectorXd Aux = VectorXd::Zero(N);
        if (si == 2 && state.vsol.wgap.size() > 0) {
            for (int idx = 0;
                 idx < N && idx < static_cast<int>(state.vsol.wgap.size());
                 ++idx) {
                Aux(idx) = state.vsol.wgap(idx);
            }
        }

        /* Determine starting index (skip if first xi is tiny - stagnation) */
        int i0 = (si < 2 && xi(0) < 1e-8 * xi(N - 1)) ? 1 : 0;

        /* First point system */
        Vector3d R1;
        Matrix38d R1_U;
        Matrix32d R1_x; /* 3x2 for x linearization */
        std::vector<int> J;

        if (si < 2) {
            /* Stagnation state for lower/upper surfaces */
            int ip0 = i0, ip1 = i0 + 1;
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip0]);
            Uin.col(1) = state.glob.U.col(Is[ip1]);
            Vector2d xin(xi(ip0), xi(ip1));

            auto stag = StagnationState(Uin, xin);

            /* Similarity station residual */
            param.turb = false;
            param.simi = true;
            Matrix42d Ust_pair;
            Ust_pair.col(0) = stag.Ust;
            Ust_pair.col(1) = stag.Ust;
            Vector2d xst_pair(stag.xst, stag.xst);
            Vector2d Aux_pair(Aux(ip0), Aux(ip0));

            auto res = ResidualStation(param, xst_pair, Ust_pair, Aux_pair);
            param.simi = false;

            /* ∂R/∂Re·Re for the simi station: res.lin_Ret is computed for
             * (Ust, Ust) inputs where Ret_1=Ret_2=Ret_m=Ret(Ust), so it
             * already represents the total Ret-derivative · Ret_Ust. */
            state.glob.dR_dRe_times_Re.segment<3>(3 * Is[i0]) += res.lin_Ret;
            state.glob.dR_dMa_times_Ma.segment<3>(3 * Is[i0]) += res.lin_Ma;

            /* Chain rule: R1_U = R1_Ut * [Ust_U]
             * R1_Ust is 3x4: sum of both 3x4 blocks from res.lin_U */
            Matrix34d R1_Ust =
                res.lin_U.block<3, 4>(0, 0) + res.lin_U.block<3, 4>(0, 4);
            R1 = res.value;
            R1_U.setZero();
            R1_U.block<3, 4>(0, 0) = R1_Ust * stag.Ust_U.block<4, 4>(0, 0);
            R1_U.block<3, 4>(0, 4) = R1_Ust * stag.Ust_U.block<4, 4>(0, 4);
            R1_x = R1_Ust * stag.Ust_x;

            J = {Is[ip0], Is[ip1]};

            /* Handle i0=1 case (stagnation on node) */
            if (i0 == 1) {
                int Ig = 3 * Is[0];
                state.glob.R.segment<3>(Ig) =
                    state.glob.U.col(Is[0]).head<3>() - stag.Ust.head<3>();
                /* d/dU[Is[0]] */
                for (int k = 0; k < 3; ++k) {
                    triplets_U.emplace_back(Ig + k, 4 * Is[0] + k, 1.0);
                }
                /* d/dU[J] */
                for (int k = 0; k < 3; ++k) {
                    for (int c = 0; c < 4; ++c) {
                        if (std::abs(stag.Ust_U(k, c)) > 0.0)
                            triplets_U.emplace_back(Ig + k, 4 * J[0] + c,
                                                    -stag.Ust_U(k, c));
                        if (std::abs(stag.Ust_U(k, c + 4)) > 0.0)
                            triplets_U.emplace_back(Ig + k, 4 * J[1] + c,
                                                    -stag.Ust_U(k, c + 4));
                    }
                }
                /* d/dx[J] (matching mfoil: R_x[Ig, J] = -Ust_x[0:3,:]) */
                for (int k = 0; k < 3; ++k) {
                    if (std::abs(stag.Ust_x(k, 0)) > 0.0)
                        triplets_x.emplace_back(Ig + k, J[0],
                                                -stag.Ust_x(k, 0));
                    if (std::abs(stag.Ust_x(k, 1)) > 0.0)
                        triplets_x.emplace_back(Ig + k, J[1],
                                                -stag.Ust_x(k, 1));
                }
            }
        } else {
            /* Wake initialization */
            auto wake = WakeSys(state, param);
            R1 = wake.R;
            R1_U.setZero();
            /* Map 3x12 to 3x8 format (we only store the needed blocks) */
            R1_U.block<3, 4>(0, 0) = wake.R_U.block<3, 4>(0, 0); /* Lower */
            R1_U.block<3, 4>(0, 4) = wake.R_U.block<3, 4>(0, 4); /* Upper */
            /* Wake contribution (indices 8-11 in R_U) stored separately */
            R1_x.setZero();
            J = {wake.J[0], wake.J[1]};
            param.turb = true;
            param.wake = true;

            /* Store first point system */
            int Ig = 3 * Is[i0];
            state.glob.R.segment<3>(Ig) = R1;
            state.glob.dR_dMa_times_Ma.segment<3>(Ig) += wake.lin_Ma_times_Ma;
            state.glob.dR_dRe_times_Re.segment<3>(Ig) += wake.lin_Re_times_Re;
            /* Add wake state contribution */
            for (int k = 0; k < 3; ++k) {
                for (int c = 0; c < 4; ++c) {
                    if (std::abs(wake.R_U(k, c)) > 0.0)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[0] + c,
                                                wake.R_U(k, c));
                    if (std::abs(wake.R_U(k, c + 4)) > 0.0)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[1] + c,
                                                wake.R_U(k, c + 4));
                    if (std::abs(wake.R_U(k, c + 8)) > 0.0)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[2] + c,
                                                wake.R_U(k, c + 8));
                }
            }
        }

        /* Store first point system for non-wake surfaces */
        if (si < 2) {
            int Ig = 3 * Is[i0];
            state.glob.R.segment<3>(Ig) = R1;
            for (size_t j = 0; j < J.size(); ++j) {
                for (int k = 0; k < 3; ++k) {
                    for (int c = 0; c < 4; ++c) {
                        double val = R1_U(k, 4 * j + c);
                        if (std::abs(val) > 0.0)
                            triplets_U.emplace_back(Ig + k, 4 * J[j] + c, val);
                    }
                }
            }
            /* R1_x linearization w.r.t. xi is handled through stagnation point
             * sens */
            for (size_t j = 0; j < J.size(); ++j) {
                for (int k = 0; k < 3; ++k) {
                    double val = R1_x(k, j);
                    if (std::abs(val) > 0.0)
                        triplets_x.emplace_back(Ig + k, J[j], val);
                }
            }
        }

        /* March over rest of points */
        for (int i = i0 + 1; i < N; ++i) {
            int ip = i - 1;
            bool tran = (Is[ip] < static_cast<int>(state.vsol.turb.size()) &&
                         Is[i] < static_cast<int>(state.vsol.turb.size())) &&
                        (state.vsol.turb[Is[ip]] != state.vsol.turb[Is[i]]);

            /* Get states */
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip]);
            Uin.col(1) = state.glob.U.col(Is[i]);
            Vector2d xin(xi(ip), xi(i));
            Vector2d Aux_in(Aux(ip), Aux(i));

            /* Compute residual - use transition residual at transition stations
             */
            ResidualLin res;
            if (tran) {
                /* Transition station: combined laminar + turbulent residual */
                res = ResidualTransition(state, param, xin, Uin, Aux_in);
                /* After transition, subsequent stations are turbulent */
                param.turb = true;
            } else {
                /* Update param for this station */
                StationParam(state, param, Is[i]);
                /* Regular residual */
                res = ResidualStation(param, xin, Uin, Aux_in);
            }

            /* Store in global system */
            int Ig = 3 * Is[i];
            state.glob.R.segment<3>(Ig) += res.value;
            state.glob.dR_dRe_times_Re.segment<3>(Ig) += res.lin_Ret;
            state.glob.dR_dMa_times_Ma.segment<3>(Ig) += res.lin_Ma;

            /* Store Jacobian entries */
            for (int k = 0; k < 3; ++k) {
                for (int c = 0; c < 4; ++c) {
                    double val1 = res.lin_U(k, c);
                    double val2 = res.lin_U(k, c + 4);
                    if (std::abs(val1) > 0.0)
                        triplets_U.emplace_back(Ig + k, 4 * Is[ip] + c, val1);
                    if (std::abs(val2) > 0.0)
                        triplets_U.emplace_back(Ig + k, 4 * Is[i] + c, val2);
                }
                double val_x1 = res.lin_x(k, 0);
                double val_x2 = res.lin_x(k, 1);
                if (std::abs(val_x1) > 0.0)
                    triplets_x.emplace_back(Ig + k, Is[ip], val_x1);
                if (std::abs(val_x2) > 0.0)
                    triplets_x.emplace_back(Ig + k, Is[i], val_x2);
            }
        }
    }

    /* Include effects of R_x into R_U: R_ue += R_x * x_st * st_ue, where
     * x_st = -sgnue is the sensitivity of xi to the stagnation point
     * (the wake moves with the upper surface). Added as triplets, which
     * setFromTriplets sums with the existing entries. */
    const int N = state.foil.N;
    VectorXd R_st = VectorXd::Zero(3 * Nsys);
    for (const Triplet& t : triplets_x) {
        const double x_st = t.col() < N ? -state.isol.sgnue(t.col()) : -1.0;
        R_st(t.row()) += t.value() * x_st;
    }
    const auto& Istag = state.isol.Istag;
    const auto& st_ue = state.isol.sstag_ue;
    for (int row = 0; row < 3 * Nsys; ++row) {
        if (R_st(row) == 0.0) continue;
        triplets_U.emplace_back(row, 4 * Istag[0] + 3, R_st(row) * st_ue(0));
        triplets_U.emplace_back(row, 4 * Istag[1] + 3, R_st(row) * st_ue(1));
    }

    state.glob.R_U.setFromTriplets(triplets_U.begin(), triplets_U.end());
    state.glob.R_x.setFromTriplets(triplets_x.begin(), triplets_x.end());
}

namespace {

/// d(ueinv)/d(alpha) at every node [Nsys], for the Cl constraint.
VectorXd UeInvAlpha(const FoilState& state) {
    const int Nsys = state.glob.Nsys;
    const double alpha = state.oper.alpha;
    VectorXd out = VectorXd::Zero(Nsys);
    for (int i = 0; i < Nsys; ++i) {
        double duinv_dalpha = 0.0;
        if (i < state.foil.N) {
            duinv_dalpha =
                state.isol.sgnue(i) * (-state.isol.gamref(i, 0) * Sind(alpha) +
                                       state.isol.gamref(i, 1) * Cosd(alpha));
        } else if (i - state.foil.N <
                   static_cast<int>(state.isol.uewiref.rows())) {
            const int iw = i - state.foil.N;
            duinv_dalpha = -state.isol.uewiref(iw, 0) * Sind(alpha) +
                           state.isol.uewiref(iw, 1) * Cosd(alpha);
        }
        out(i) = duinv_dalpha * kDegToRad;
    }
    return out;
}

/// Edge velocity with the floor used in the coupled Jacobian.
VectorXd ClippedUe(const FoilState& state) {
    VectorXd ue = state.glob.U.row(3).transpose();
    const double uemax = ue.cwiseAbs().maxCoeff();
    for (int i = 0; i < ue.size(); ++i) ue(i) = std::max(ue(i), 1e-10 * uemax);
    return ue;
}

/**
 * @brief Coupled residual [R_BL; R_ue; (cl - cltgt)] at the current state.
 * @param ue Clipped edge velocity (see ClippedUe)
 */
VectorXd AugmentedResidual(const FoilState& state, const VectorXd& ue) {
    const int Nsys = state.glob.Nsys;
    const bool docl = state.oper.givencl;
    const VectorXd ds = state.glob.U.row(1).transpose();
    VectorXd R(4 * Nsys + (docl ? 1 : 0));
    R.head(3 * Nsys) = state.glob.R;
    R.segment(3 * Nsys, Nsys) =
        ue - (GetUeInv(state) + state.vsol.ue_m * ds.cwiseProduct(ue));
    if (docl) R(R.size() - 1) = state.post.cl - state.oper.cltgt;
    return R;
}

/**
 * @brief Assemble and factor the full augmented Jacobian R_V with an
 * equilibrated sparse LU. Used by the trust-region step, which needs R_V
 * and the scales.
 */
void SolveGlobAssembled(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    const bool docl = state.oper.givencl;
    const VectorXd ue = ClippedUe(state);
    const VectorXd ds = state.glob.U.row(1).transpose();
    VectorXd R = AugmentedResidual(state, ue);
    const int NN = static_cast<int>(R.size());
    state.glob.R_V.resize(NN, NN);

    /* Build Jacobian using triplets */
    std::vector<Triplet> triplets;
    triplets.reserve(state.glob.R_U.nonZeros() + 4 * Nsys);

    /* Copy existing R_U entries */
    for (int k = 0; k < state.glob.R_U.outerSize(); ++k) {
        for (SparseMatrixd::InnerIterator it(state.glob.R_U, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    /* Fourth equation Jacobian
     * d(R_ue)/d(ue) = I - ue_m * diag(ds)
     * d(R_ue)/d(ds) = -ue_m * diag(ue) */
    for (int i = 0; i < Nsys; ++i) {
        int row = 3 * Nsys + i;
        int ue_idx = 4 * i + 3; /* ue index in state vector */

        /* Identity part for ue */
        triplets.emplace_back(row, ue_idx, 1.0);

        /* ue_m contributions */
        for (int j = 0; j < Nsys; ++j) {
            double ue_m_ij = state.vsol.ue_m(i, j);
            if (std::abs(ue_m_ij) > 0.0) {
                /* d/d(ue_j): -ue_m(i,j) * ds(j) */
                triplets.emplace_back(row, 4 * j + 3, -ue_m_ij * ds(j));
                /* d/d(ds_j): -ue_m(i,j) * ue(j) */
                triplets.emplace_back(row, 4 * j + 1, -ue_m_ij * ue(j));
            }
        }
    }

    /* Cl-alpha constraint if requested */
    if (docl) {
        /* d(cl)/d(alpha) */
        triplets.emplace_back(NN - 1, NN - 1, state.post.cl_alpha);

        /* d(cl)/d(ue) at airfoil nodes */
        int N = state.foil.N;
        for (int i = 0; i < N; ++i) {
            if (i < static_cast<int>(state.post.cl_ue.size())) {
                triplets.emplace_back(NN - 1, 4 * i + 3, state.post.cl_ue(i));
            }
        }

        const VectorXd ueinv_alpha = UeInvAlpha(state);
        for (int i = 0; i < Nsys; ++i)
            triplets.emplace_back(3 * Nsys + i, NN - 1, -ueinv_alpha(i));
    }

    /* Build sparse matrix */
    state.glob.R_V.setFromTriplets(triplets.begin(), triplets.end());

    // Equilibrate only the linear solve. Keep the physical Jacobian and
    // residual intact for convergence checks and adjoint factorization.
    state.glob.coupled_residual = R;
    auto& uscale = state.glob.variable_scale;
    auto& rscale = state.glob.residual_scale;
    uscale = VectorXd::Ones(NN);
    rscale = VectorXd::Ones(NN);
    SparseMatrixd linear_system = state.glob.R_V;
    if (state.param.scale_coupled_system) {
        const double thickness = 1e-4 * state.geom.chord;
        for (int i = 0; i < Nsys; ++i) {
            uscale(4 * i) = std::max(std::abs(state.glob.U(0, i)), thickness);
            uscale(4 * i + 1) =
                std::max(std::abs(state.glob.U(1, i)), thickness);
            uscale(4 * i + 2) = state.vsol.turb[i] ? 0.03 : 1.0;
            uscale(4 * i + 3) = state.oper.Vinf;
        }
        VectorXd row_norm = VectorXd::Zero(NN);
        for (int col = 0; col < linear_system.outerSize(); ++col)
            for (SparseMatrixd::InnerIterator it(linear_system, col); it; ++it)
                row_norm(it.row()) += std::abs(it.value() * uscale(it.col()));
        for (int row = 0; row < NN; ++row)
            rscale(row) = 1.0 / std::max(row_norm(row), 1e-30);
        for (int col = 0; col < linear_system.outerSize(); ++col)
            for (SparseMatrixd::InnerIterator it(linear_system, col); it; ++it)
                it.valueRef() *= rscale(it.row()) * uscale(it.col());
    }

    /* Solve sparse system */
    Eigen::SparseLU<SparseMatrixd> solver;
    solver.compute(linear_system);
    if (solver.info() != Eigen::Success) {
        std::cerr << "  Warning: SparseLU decomposition failed" << '\n';
        throw std::runtime_error(
            "SolveGlob: linear factorization/solve failed");
    }

    VectorXd dV = -solver.solve(rscale.cwiseProduct(R));
    dV = uscale.cwiseProduct(dV);
    if (solver.info() != Eigen::Success || !dV.allFinite()) {
        std::cerr << "  Warning: SparseLU solve failed" << '\n';
        throw std::runtime_error(
            "SolveGlob: linear factorization/solve failed");
    }

    /* Store dU (reshaped from 4*Nsys vector to 4xNsys matrix) */
    if (state.glob.dU.cols() != Nsys) {
        state.glob.dU.resize(4, Nsys);
    }
    for (int i = 0; i < Nsys; ++i) {
        state.glob.dU.col(i) = dV.segment<4>(4 * i);
    }

    if (docl) {
        state.glob.dalpha = dV(NN - 1);
    } else {
        state.glob.dalpha = 0.0;
    }
}

}  // namespace

/**
 * @brief Solve the coupled Newton system for dU (and dalpha with a Cl target).
 *
 * Variables are [th, ds, sa, ue] at each node; the fourth equation is
 * ue = ueinv + ue_m (ds ue). The system is factored through the mass-defect
 * Schur complement (CoupledFactorization), which costs O(N^2) plus one dense
 * N x N LU per iteration. With a Cl target, the alpha row and column are a
 * bordering handled with a second solve.
 *
 * Pseudo-transient continuation adds (1/tau) to each equation's own variable:
 * boundary-layer equation k of node i on variable k of node i, and the
 * edge-velocity equation on ue. At convergence the term multiplies a zero
 * update, so the root and the adjoint are unchanged.
 * @param state Foil state
 */
void SolveGlob(FoilState& state) {
    if (state.param.use_trust_region) {
        SolveGlobAssembled(state);
        return;
    }
    const int Nsys = state.glob.Nsys;
    const bool docl = state.oper.givencl;
    const VectorXd ue = ClippedUe(state);
    const VectorXd ds = state.glob.U.row(1).transpose();
    const VectorXd R = AugmentedResidual(state, ue);
    state.glob.coupled_residual = R;

    const double lambda = (state.param.use_ptc && state.glob.tau > 0.0)
                              ? 1.0 / state.glob.tau
                              : 0.0;
    SparseMatrixd bl = state.glob.R_U;
    if (lambda > 0.0) {
        std::vector<Triplet> diag;
        diag.reserve(3 * Nsys);
        for (int i = 0; i < Nsys; ++i)
            for (int k = 0; k < 3; ++k)
                diag.emplace_back(3 * i + k, 4 * i + k, lambda);
        SparseMatrixd shift(3 * Nsys, 4 * Nsys);
        shift.setFromTriplets(diag.begin(), diag.end());
        bl += shift;
    }

    CoupledFactorization factorization;
    if (!factorization.Compute(bl, state.vsol.ue_m, ue, ds,
                               VectorXd::Constant(Nsys, 1.0 + lambda))) {
        throw std::runtime_error("SolveGlob: linear factorization failed");
    }
    const VectorXd x1 = factorization.Solve(-R.head(4 * Nsys));
    VectorXd dV = x1;
    state.glob.dalpha = 0.0;
    if (docl) {
        /* [J b; c^T d] [dV; da] = -[R; r_cl] */
        VectorXd b = VectorXd::Zero(4 * Nsys);
        b.tail(Nsys) = -UeInvAlpha(state);
        VectorXd c = VectorXd::Zero(4 * Nsys);
        for (int i = 0;
             i < state.foil.N && i < static_cast<int>(state.post.cl_ue.size());
             ++i)
            c(4 * i + 3) = state.post.cl_ue(i);
        const VectorXd x2 = factorization.Solve(b);
        const double denom = state.post.cl_alpha - c.dot(x2);
        if (!std::isfinite(denom) || denom == 0.0)
            throw std::runtime_error("SolveGlob: singular Cl constraint");
        state.glob.dalpha = (-R(R.size() - 1) - c.dot(x1)) / denom;
        dV -= state.glob.dalpha * x2;
    }
    if (!dV.allFinite())
        throw std::runtime_error("SolveGlob: linear solve failed");
    state.glob.dU = Eigen::Map<const MatrixXd>(dV.data(), 4, Nsys);
}

bool TakeTrustRegionStep(FoilState& state) {
    const VectorXd x = PackedState(state);
    const VectorXd uscale = state.glob.variable_scale;
    const VectorXd rscale = state.glob.residual_scale;
    const VectorXd f = rscale.cwiseProduct(CoupledResidual(state));
    SparseMatrixd a = state.glob.R_V;
    for (int col = 0; col < a.outerSize(); ++col)
        for (SparseMatrixd::InnerIterator it(a, col); it; ++it)
            it.valueRef() *= rscale(it.row()) * uscale(it.col());
    VectorXd newton(x.size());
    newton.head(state.glob.dU.size()) =
        Eigen::Map<const VectorXd>(state.glob.dU.data(), state.glob.dU.size());
    if (state.oper.givencl) newton(newton.size() - 1) = state.glob.dalpha;
    newton = newton.cwiseQuotient(uscale);
    const VectorXd gradient = a.transpose() * f;
    const double curvature = (a * gradient).squaredNorm();
    const double phi = 0.5 * f.squaredNorm();
    if (!std::isfinite(phi) || !(curvature > 0.0)) return false;
    const VectorXd cauchy = -(gradient.squaredNorm() / curvature) * gradient;
    double radius = state.glob.trust_radius;
    for (int attempt = 0; attempt < state.param.trust_max_trials; ++attempt) {
        VectorXd step;
        if (newton.norm() <= radius) {
            step = newton;
        } else if (cauchy.norm() >= radius) {
            step = (radius / cauchy.norm()) * cauchy;
        } else {
            const VectorXd edge = newton - cauchy;
            const double b = cauchy.dot(edge);
            const double t =
                (-b + std::sqrt(b * b +
                                edge.squaredNorm() *
                                    (radius * radius - cauchy.squaredNorm()))) /
                edge.squaredNorm();
            step = cauchy + t * edge;
        }
        FoilState trial = state;
        const VectorXd physical_step = uscale.cwiseProduct(step);
        trial.glob.dU = Eigen::Map<const MatrixXd>(physical_step.data(), 4,
                                                   state.glob.Nsys);
        trial.glob.dalpha =
            state.oper.givencl ? physical_step(physical_step.size() - 1) : 0.0;
        UpdateState(trial);
        // The local model uses the step before any representation changes
        // (ue sign flips, amplification -> turbulent shear). The actual
        // reduction includes those changes and all rebuilt dependent data.
        const VectorXd actual_step =
            (PackedState(trial) - x).cwiseQuotient(uscale);
        const double predicted =
            phi - 0.5 * (f + a * actual_step).squaredNorm();
        double actual = -std::numeric_limits<double>::infinity();
        if (ValidTrialState(trial)) {
            try {
                if (state.param.trust_update_topology) {
                    StagPointMove(trial);
                    UpdateTransition(trial);
                }
                BuildGlobSys(trial);
                CalcForce(trial);
                actual = phi - 0.5 * rscale.cwiseProduct(CoupledResidual(trial))
                                         .squaredNorm();
            } catch (const std::exception& error) {
                // Reject recoverable domain/linear-solve errors in a trial.
                // No mutation of the accepted state has taken place.
                actual = -std::numeric_limits<double>::infinity();
                if (state.param.verb >= 1)
                    std::cout << "  Trust trial rejected: " << error.what()
                              << '\n';
            }
        }
        const double ratio = predicted > 0.0 ? actual / predicted : -1.0;
        if (std::isfinite(ratio) && ratio > 1e-4 && actual > 0.0) {
            if (ratio < 0.25)
                radius *= 0.25;
            else if (ratio > 0.75 && step.norm() > 0.8 * radius)
                radius = std::min(2.0 * radius, 1e4);
            trial.glob.trust_radius = radius;
            trial.glob.trust_last_ratio = ratio;
            ++trial.glob.trust_accepted;
            state = std::move(trial);
            return true;
        }
        ++state.glob.trust_rejected;
        radius *= 0.25;
    }
    state.glob.trust_radius = radius;
    return false;
}

/**
 * @brief Update state with adaptive under-relaxation (matching mfoil approach)
 * @param state Foil state
 */
namespace {

/**
 * @brief mfoil under-relaxation: absolute caps on amplification, shear
 * stress and ue changes, and at most a 50% decrease in theta and delta*.
 */
double AbsoluteRelaxation(const FoilState& state) {
    const int Nsys = state.glob.Nsys;
    double omega = 1.0;

    /* Get max ctau for threshold calculations */
    double ctmax = 0.0;
    std::vector<int> turbIdx;
    for (int i = 0; i < Nsys; ++i) {
        if (i < static_cast<int>(state.vsol.turb.size()) &&
            state.vsol.turb[i]) {
            ctmax = std::max(ctmax, state.glob.U(2, i));
            turbIdx.push_back(i);
        }
    }

    /* 1. Limit theta and delta* decreases to 50% */
    for (int k = 0; k < 2; ++k) {
        double fmin = 0.0;
        for (int i = 0; i < Nsys; ++i) {
            double Uk = state.glob.U(k, i);
            double dUk = state.glob.dU(k, i);
            if (Uk > 1e-12) {
                fmin = std::min(fmin, dUk / Uk);
            }
        }
        if (fmin < -0.5) {
            double om = std::abs(0.5 / fmin);
            if (om < omega) {
                omega = om;
            }
        }
    }

    /* 2. Prevent negative amplification/ctau (mfoil approach) */
    for (int i = 0; i < Nsys; ++i) {
        if (i >= static_cast<int>(state.vsol.turb.size())) continue;
        double Uk = state.glob.U(2, i);
        double dUk = state.glob.dU(2, i);
        /* Skip small values (don't limit small amp/ctau) */
        if (!state.vsol.turb[i] && Uk < 0.2) continue;        /* small amp */
        if (state.vsol.turb[i] && Uk < 0.1 * ctmax) continue; /* small ctau */
        if (Uk == 0.0 || dUk == 0.0) continue;
        /* Check if update would make value negative */
        if (Uk + dUk < 0.0) {
            double om = 0.8 * std::abs(Uk / dUk);
            if (om < omega) {
                omega = om;
            }
        }
    }

    /* 3. Limit amplification changes (laminar regions) to 2.0
     * Include all laminar points (matching mfoil) */
    double dumax_amp = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        if (i < static_cast<int>(state.vsol.turb.size()) &&
            !state.vsol.turb[i]) {
            dumax_amp = std::max(dumax_amp, std::abs(state.glob.dU(2, i)));
        }
    }
    if (dumax_amp > 0) {
        double om = 2.0 / dumax_amp;
        if (om < omega) {
            omega = om;
        }
    }

    /* 4. Limit ctau changes (turbulent regions) to 0.05
     * Include all turbulent points (matching mfoil) */
    double dumax_ct = 0.0;
    for (int i : turbIdx) {
        dumax_ct = std::max(dumax_ct, std::abs(state.glob.dU(2, i)));
    }
    if (dumax_ct > 0) {
        double om = 0.05 / dumax_ct;
        if (om < omega) {
            omega = om;
        }
    }

    /* 5. Limit ue changes to 0.2 * Vinf */
    double fmax_ue = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        fmax_ue =
            std::max(fmax_ue, std::abs(state.glob.dU(3, i)) / state.oper.Vinf);
    }
    if (fmax_ue > 0.2) {
        double om = 0.2 / fmax_ue;
        if (om < omega) {
            omega = om;
        }
    }

    return omega;
}

/**
 * @brief XFOIL under-relaxation (UPDATE in xbl.f): normalize every change
 * (amplification by 10, shear stress, theta and delta* by their values, ue
 * by 0.25 Vinf) and scale the step so that none exceeds +1.5 or -0.5.
 */
double RelativeRelaxation(const FoilState& state) {
    constexpr double kHigh = 1.5;
    constexpr double kLow = -0.5;
    double omega = 1.0;
    auto limit = [&](double change) {
        if (omega * change > kHigh) omega = kHigh / change;
        if (omega * change < kLow) omega = kLow / change;
    };
    const auto& U = state.glob.U;
    const auto& dU = state.glob.dU;
    const auto& turb = state.vsol.turb;
    for (int i = 0; i < state.glob.Nsys; ++i) {
        const bool turbulent = i < static_cast<int>(turb.size()) && turb[i];
        if (turbulent) {
            if (U(2, i) > 0.0) limit(dU(2, i) / U(2, i));
        } else {
            limit(dU(2, i) / 10.0);
        }
        if (U(0, i) > 0.0) limit(dU(0, i) / U(0, i));
        if (U(1, i) > 0.0) limit(dU(1, i) / U(1, i));
        limit(std::abs(dU(3, i)) / (0.25 * state.oper.Vinf));
    }
    return omega;
}

}  // namespace

void UpdateState(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    double omega = state.param.step_limiter == StepLimiter::Relative
                       ? RelativeRelaxation(state)
                       : AbsoluteRelaxation(state);

    std::vector<int> turbIdx;
    for (int i = 0; i < Nsys; ++i)
        if (i < static_cast<int>(state.vsol.turb.size()) && state.vsol.turb[i])
            turbIdx.push_back(i);

    /* Limit alpha changes to 2 degrees */
    if (std::abs(state.glob.dalpha) > 2.0) {
        omega = std::min(omega, 2.0 / std::abs(state.glob.dalpha));
    }

    /* Apply update */
    double alpha_change = omega * state.glob.dalpha;
    state.glob.U += omega * state.glob.dU;
    state.oper.alpha += alpha_change;

    /* Fix bad Hk after update (using proper Hk calculation) */
    for (int si = 0; si < 3; ++si) {
        double Hkmin = (si == 2) ? 1.00005 : 1.02;
        const auto& Is = state.vsol.Is[si];
        Param param = BuildParam(state, si);
        for (int idx : Is) {
            if (idx >= Nsys) continue;
            StationParam(state, param, idx);
            Vector4d Uj = state.glob.U.col(idx);
            auto [Hk, Hk_U] = GetHk(Uj, param);
            if (Hk < Hkmin) {
                state.glob.U(1, idx) +=
                    2.0 * (Hkmin - Hk) * state.glob.U(1, idx);
            }
        }
    }

    /* Fix negative ctau after update */
    double ctmax = 0.0;
    for (int i : turbIdx) {
        ctmax = std::max(ctmax, state.glob.U(2, i));
    }
    for (int i : turbIdx) {
        if (state.glob.U(2, i) < 0) {
            state.glob.U(2, i) = 0.1 * ctmax;
        }
    }

    /* Ensure positivity constraints */
    for (int i = 0; i < Nsys; ++i) {
        state.glob.U(0, i) = std::max(state.glob.U(0, i), 1e-10);
        state.glob.U(1, i) = std::max(state.glob.U(1, i), 1e-10);
    }

    /* Rebuild inviscid solution if alpha changed */
    if (std::abs(alpha_change) > 1e-10) {
        RebuildIsol(state);
    }
}

/**
 * @brief Calculate d(ue)/d(mass) sensitivities
 * This implements the transpiration BC coupling
 * @param state Foil state
 */
void CalcUeM(FoilState& state) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    if (Nw == 0) {
        /* No wake - just use identity */
        state.vsol.ue_m.resize(Nsys, Nsys);
        state.vsol.ue_m.setIdentity();
        return;
    }

    /* 1. Cgam = d(wake uei)/d(gamma) [Nw x N]
     * Use full InviscidVelocity to get linearization w.r.t. gamma, which
     * includes TE source and TE vortex contributions (matching mfoil's
     * inviscid_velocity call) */
    MatrixXd Cgam = MatrixXd::Zero(Nw, N);
    for (int i = 0; i < Nw; ++i) {
        Vector2d xi = state.wake.x.col(i);
        Vector2d ti = state.wake.t.col(i);

        MatrixXd V_G;
        InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi, true,
                         &V_G);

        /* Cgam(i,:) = V_G[0,:]*ti[0] + V_G[1,:]*ti[1] (tangent component) */
        for (int j = 0; j < N; ++j) {
            Cgam(i, j) = V_G(0, j) * ti(0) + V_G(1, j) * ti(1);
        }
    }

    /* 2. B = d(airfoil streamfunction)/d(source) [(N+1) x (N+Nw-2)] */
    int Npanel = N + Nw - 2; /* N-1 airfoil panels + Nw-1 wake panels */
    MatrixXd B = MatrixXd::Zero(N + 1, Npanel);

    for (int i = 0; i < N; ++i) {
        Vector2d xi = state.foil.x.col(i);

        /* Constant source panels on airfoil */
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);
            B(i, j) = PanelConstSourceStream(Xj, xi);
        }

        /* Linear source panels on wake */
        for (int j = 0; j < Nw - 1; ++j) {
            Matrix2d Xj;
            Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
            Xj.col(0) = state.wake.x.col(j);
            Xj.col(1) = Xm;

            /* Ghost extension at last point */
            Vector2d X_right;
            if (j == Nw - 2) {
                X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
            } else {
                X_right = state.wake.x.col(j + 1);
            }

            double a_left, b_left, a_right, b_right;

            /* Left half panel */
            PanelLinSourceStream(Xj, xi, a_left, b_left);
            if (j > 0) {
                B(i, N - 1 + j) += 0.5 * a_left + b_left;
                B(i, N - 1 + j - 1) += 0.5 * a_left;
            } else {
                B(i, N - 1 + j) += b_left;
            }

            /* Right half panel */
            Xj.col(0) = Xm;
            Xj.col(1) = X_right;
            PanelLinSourceStream(Xj, xi, a_right, b_right);
            B(i, N - 1 + j) += a_right + 0.5 * b_right;
            if (j < Nw - 2) {
                B(i, N - 1 + j + 1) += 0.5 * b_right;
            } else {
                B(i, N - 1 + j) += 0.5 * b_right;
            }
        }
    }

    /* 3. Bp = -inv(AIC) * B [N x (N+Nw-2)]
     * Note: AIC is (N+1) x (N+1), but last row is Kutta condition */
    MatrixXd Bp = -SolveAIC(state, B);
    Bp.conservativeResize(N, Npanel); /* Remove last row */

    /* 4. Csig = d(wake uei)/d(source) [Nw x (N+Nw-2)] */
    MatrixXd Csig = MatrixXd::Zero(Nw, Npanel);

    for (int i = 0; i < Nw; ++i) {
        Vector2d xi = state.wake.x.col(i);
        Vector2d ti = state.wake.t.col(i);

        /* Constant source panels on airfoil */
        int jstart = (i == 0) ? 1 : 0;
        int jend = (i == 0) ? N - 2 : N - 1;
        for (int j = jstart; j < jend; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);
            Vector2d v = PanelConstSourceVelocity(Xj, xi, &ti);
            Csig(i, j) = v(0);
        }

        /* Piecewise linear source on wake */
        for (int j = 0; j < Nw; ++j) {
            int il = std::max(j - 1, 0);
            int im = j;
            int ir = std::min(j + 1, Nw - 1);

            Vector2d Xl = 0.5 * (state.wake.x.col(il) + state.wake.x.col(im));
            Vector2d Xm = state.wake.x.col(im);
            Vector2d Xr = 0.5 * (state.wake.x.col(im) + state.wake.x.col(ir));

            if (j == Nw - 1) {
                Xr = 2.0 * Xm - Xl; /* Ghost extension */
            }

            double d1 = (Xm - Xl).norm();
            double d2 = (Xr - Xm).norm();

            if (i == j) {
                /* Self-influence (singular) */
                if (j == 0) {
                    /* First point: TE special case */
                    double dl =
                        (state.foil.x.col(1) - state.foil.x.col(0)).norm();
                    double du =
                        (state.foil.x.col(N - 1) - state.foil.x.col(N - 2))
                            .norm();
                    Csig(i, 0) += (0.5 / kPi) * (SafeLog(dl / d2) + 1.0);
                    Csig(i, N - 2) += (0.5 / kPi) * (SafeLog(du / d2) + 1.0);
                    Csig(i, N - 1) += -0.5 / kPi;
                } else if (j == Nw - 1) {
                    /* Last point: no self effect */
                } else {
                    double aa = (0.25 / kPi) * SafeLog(d1 / d2);
                    Csig(i, N - 1 + j - 1) += aa + 0.5 / kPi;
                    Csig(i, N - 1 + j) += aa - 0.5 / kPi;
                }
            } else {
                /* Non-self influence */
                if (j == 0) {
                    Matrix2d Xj;
                    Xj.col(0) = Xm;
                    Xj.col(1) = Xr;
                    Vector2d a, b;
                    PanelLinSourceVelocity(Xj, xi, a, b);
                    Csig(i, N - 1) += b.dot(ti);
                    Csig(i, 0) += a.dot(ti);
                    Csig(i, N - 2) += a.dot(ti);
                } else if (j == Nw - 1) {
                    Matrix2d Xj;
                    Xj.col(0) = Xl;
                    Xj.col(1) = Xr;
                    Vector2d v = PanelConstSourceVelocity(Xj, xi, &ti);
                    Csig(i, N + Nw - 3) += v(0);
                } else {
                    Matrix2d Xj;
                    Vector2d a1, b1, a2, b2;

                    /* Left half-panel */
                    Xj.col(0) = Xl;
                    Xj.col(1) = Xm;
                    PanelLinSourceVelocity(Xj, xi, a1, b1);

                    /* Right half-panel */
                    Xj.col(0) = Xm;
                    Xj.col(1) = Xr;
                    PanelLinSourceVelocity(Xj, xi, a2, b2);

                    Csig(i, N - 1 + j - 1) += a1.dot(ti) + 0.5 * b1.dot(ti);
                    Csig(i, N - 1 + j) += 0.5 * a2.dot(ti) + b2.dot(ti);
                }
            }
        }
    }

    /* 5. ue_sigma = d(unsigned ue)/d(source) [(N+Nw) x (N+Nw-2)]
     * Df = Bp for foil, Dw = Cgam*Bp + Csig for wake */
    MatrixXd CgamBp = Cgam * Bp;
    MatrixXd Dw = CgamBp + Csig;
    Dw.row(0) = Bp.row(N - 1); /* First wake point same as TE */

    state.vsol.ue_sigma.resize(Nsys, Npanel);
    state.vsol.ue_sigma.topRows(N) = Bp;
    state.vsol.ue_sigma.bottomRows(Nw) = Dw;

    /* 6. sigma_m = d(source)/d(mass) [(N+Nw-2) x (N+Nw)] */
    std::vector<Triplet> triplets;
    for (int i = 0; i < N - 1; ++i) {
        double ds = state.foil.s(i + 1) - state.foil.s(i);
        if (ds > 1e-12) {
            triplets.emplace_back(i, i, -state.isol.sgnue(i) / ds);
            triplets.emplace_back(i, i + 1, state.isol.sgnue(i + 1) / ds);
        }
    }
    for (int i = 0; i < Nw - 1; ++i) {
        double ds = state.wake.s(i + 1) - state.wake.s(i);
        if (ds > 1e-12) {
            triplets.emplace_back(N - 1 + i, N + i, -1.0 / ds);
            triplets.emplace_back(N - 1 + i, N + i + 1, 1.0 / ds);
        }
    }
    state.vsol.sigma_m.resize(Npanel, Nsys);
    state.vsol.sigma_m.setFromTriplets(triplets.begin(), triplets.end());

    /* 7. ue_m = sgnue * ue_sigma * sigma_m */
    VectorXd sgue(Nsys);
    sgue.head(N) = state.isol.sgnue;
    sgue.tail(Nw).setOnes();

    state.vsol.ue_m =
        sgue.asDiagonal() * state.vsol.ue_sigma * state.vsol.sigma_m;
}

/**
 * @brief Rebuild sensitivities after stagnation point move
 * @param state Foil state
 */
void RebuildUeM(FoilState& state) { CalcUeM(state); }

/**
 * @brief Factor the coupled Jacobian at the converged state for adjoint
 * transpose solves (state.glob.coupled_factorization). The factorization is
 * of the unregularized Jacobian: no pseudo-time terms.
 * @param state Foil state
 */
void StoreAdjointFactorization(FoilState& state) {
    /* Rebuild global system at converged state to get current R_U */
    BuildGlobSys(state);
    const int Nsys = state.glob.Nsys;
    auto factorization = std::make_shared<CoupledFactorization>();
    // A ridge solves a different derivative problem. Fail instead of returning
    // a gradient of a regularized model that Python mfoil did not evaluate.
    const bool ok = factorization->Compute(
        state.glob.R_U, state.vsol.ue_m, ClippedUe(state),
        state.glob.U.row(1).transpose(), VectorXd::Ones(Nsys));
    state.glob.coupled_factorization = ok ? factorization : nullptr;
}

/**
 * @brief Rebuilds inviscid solution after angle of attack change
 * Recompute gamma from reference solutions
 * @param state Foil state
 */
void RebuildIsol(FoilState& state) {
    double alpha = state.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    /* Rebuild gamma on airfoil: gam = gamref[:,0]*cos(alpha) +
     * gamref[:,1]*sin(alpha) */
    int N = state.foil.N;
    for (int i = 0; i < N; ++i) {
        state.isol.gam(i) =
            state.isol.gamref(i, 0) * cosa + state.isol.gamref(i, 1) * sina;
    }

    /* For viscous flow, may need to rebuild wake and ue_m */
    if (state.oper.viscous && state.oper.redowake) {
        BuildWake(state);
        IdentifySurfaces(state);
        CalcUeM(state);
    }
}

} /* namespace scirocco */
