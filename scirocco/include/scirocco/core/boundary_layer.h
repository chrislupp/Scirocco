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
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/state.h"
#include "scirocco/util/linearization.h"

namespace scirocco {

/**
 * @brief Initialize boundary layer from inviscid edge velocity
 * @param state Foil state to initialize
 */
void InitBoundaryLayer(FoilState& state);

/**
 * @brief Compute residual at a station (non-transition)
 * @param param Parameter structure for the surface
 * @param x Station arc-length coordinates
 * @param U State matrix at the station
 * @param Aux Auxiliary variables at the station
 * @return Residual and linearization at the station
 */
ResidualLin ResidualStation(const Param& param, const Vector2d& x,
                            const Matrix42d& U, const Vector2d& Aux);

/**
 * @brief Compute residual at a transition station (combined laminar + turbulent)
 * U1 should be laminar, U2 should be turbulent
 * Returns combined residual and Jacobian
 * @param state Foil state
 * @param param_base Base parameter structure for the surface
 * @param x Station arc-length coordinates
 * @param U State matrix at the station
 * @param Aux Auxiliary variables at the station
 * @return Combined residual and Jacobian at the transition station
 */
ResidualLin ResidualTransition(FoilState& state, const Param& param_base,
                               const Vector2d& x, const Matrix42d& U,
                               const Vector2d& Aux);

/**
 * @brief Thwaites correlation for stagnation initialization
 * Returns (theta, delta*) at stagnation point
 * @param K Thwaites parameter at stagnation
 * @param nu Kinematic viscosity
 * @return Pair of (momentum thickness theta, displacement thickness delta*) at stagnation point
 */
std::pair<double, double> ThwaitesInit(double K, double nu);

/**
 * @brief Build parameter structure for a surface
 * @param state Foil state
 * @param si Surface index
 * @return Parameter structure for the specified surface
 */
Param BuildParam(const FoilState& state, int si);

/**
 * @brief Update parameter for a specific station
 * @param state Foil state
 * @param param Parameter structure to update
 * @param i Station index
 */
void StationParam(const FoilState& state, Param& param, int i);

/**
 * @brief Result for stagnation state extrapolation
 */
struct StagnationResult {
    Vector4d Ust;     /* Stagnation state [th, ds, sa, ue] */
    Matrix48d Ust_U;  /* Linearization w.r.t. U1 and U2 (4x8) */
    Matrix42d Ust_x;  /* Linearization w.r.t. x1 and x2 (4x2) */
    double xst;       /* Stagnation point location (small) */
};

/**
 * @brief Extrapolate two states to stagnation point
 * @param U State matrix containing the two bracketing states
 * @param x Arc-length coordinates of the two bracketing stations
 * @return Stagnation state and its linearizations
 */
StagnationResult StagnationState(const Matrix42d& U, const Vector2d& x);

/**
 * @brief Result for wake system
 */
struct WakeResult {
    Vector3d R;      /* Residual [th, ds, sa] */
    Matrix312d R_U;  /* Linearization w.r.t. lower, upper, wake states (3x12) */
    std::array<int, 3> J;  /* Indices: lower TE, upper TE, first wake */
    /* ∂R/∂Ma · Ma, non-zero only when at least one TE side is laminar
     * (because ctw is then built from GetCttr which is Ma-dependent). */
    Vector3d lin_Ma_times_Ma{Vector3d::Zero()};
};

/**
 * @brief Construct residual system for wake initialization
 * @param state Foil state
 * @param param Parameter structure for the wake surface
 * @return Wake residual, linearization, and station indices
 */
WakeResult WakeSys(const FoilState& state, Param& param);

}  /* namespace scirocco */
