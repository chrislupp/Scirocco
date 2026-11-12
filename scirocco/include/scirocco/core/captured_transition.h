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
*/
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/param.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @file captured_transition.h
 * @brief Captured laminar-turbulent transition (TransitionModel::Captured).
 *
 * The cut cell switches a station from laminar to turbulent equations when
 * the transition point crosses a node, so outputs are continuous but their
 * derivatives jump. Here the amplification factor n and the root shear
 * stress ct live at every node, and each interval carries the same four
 * equations, discretized as ResidualStation:
 *
 *   - momentum, shape: closures Hs, cf, cDi blended per station,
 *     f = (1 - g) f_lam + g f_turb, g = S((n - nc) / wc)
 *   - amplification: dn/dx = (1 - gs) damp + gs smax(damp, C / th),
 *     gs = S((n - ncrit) / ws)
 *   - shear lag: the standard lag equation for ct, with equilibrium target
 *     (1 - r) cttr + r cteq and its pressure-gradient term scaled by r,
 *     r = S((n - ncrit + 2.5 ws) / ws)
 *
 * S(t) = (1 + tanh t) / 2; in the wake g = r = 1 and n2 = n1. Onset is at
 * n = ncrit, as with the cut cell; the closures then blend over about
 * capt_L momentum thicknesses (C = 2.2 wc / capt_L, nc = ncrit + b ws +
 * 2.3 wc). This mirrors mfoil/mfoil_captured.py, which documents the
 * formulation and its evidence (scirocco/validation/captured_transition).
 */

/** Interval residual of the captured model. */
struct CapturedResidual {
    Eigen::Vector4d R;                /**< [mom, shape, amp, lag] */
    Eigen::Matrix<double, 4, 10> R_U; /**< w.r.t. [U1, U2], 5 states each */
    Eigen::Matrix<double, 4, 2> R_x;  /**< w.r.t. [x1, x2] */
};

/**
 * @brief Captured residual of one interval.
 * @param param boundary-layer parameters (turb/wake/simi as for the station)
 * @param x xi at the two stations
 * @param U states [th, ds, ct, ue, n] at the two stations (5 x 2)
 * @param Aux wake gap at the two stations
 */
CapturedResidual CapturedIntervalResidual(const Param& param,
                                          const Eigen::Vector2d& x,
                                          const Eigen::Matrix<double, 5, 2>& U,
                                          const Eigen::Vector2d& Aux);

/**
 * @brief Converge the captured system, starting from a converged cut-cell
 * viscous solution in state. On success glob.U row 2 holds ct at every node,
 * glob.Un holds n, vsol.turb marks g > 1/2, and glob.conv is true.
 */
void SolveCaptured(FoilState& state);

}  // namespace scirocco
