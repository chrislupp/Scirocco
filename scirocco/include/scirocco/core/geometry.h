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

#include <string>

#include "scirocco/config.h"
#include "scirocco/data/geom.h"

namespace scirocco {

/** @brief Forward declaration */
struct FoilState;

/** @brief Trailing edge information result */
struct TEInfo {
    Vector2d t;  /**< Unit tangent at TE (average of upper/lower) */
    double hTE;  /**< TE gap height */
    double dtdx; /**< d(theta)/dx at TE */
    double tcp;  /**< TE constant source panel strength factor */
    double tdp;  /**< TE vortex panel strength factor */
};

/** @brief Generate NACA 4 or 5-digit airfoil points
 *  @param geom   Geom structure with point coordinates filled in on output
 *  @param digits NACA digit string (e.g., "0012", "23012")
 */
void NacaPoints(Geom& geom, const std::string& digits);

/** @brief Set coordinates from user-provided points
 *  @param geom   Geom structure to fill
 *  @param coords Point coordinates [2 x N], ordered from TE lower to TE upper
 */
void SetCoords(Geom& geom, const MatrixXd& coords);

/** @brief Get trailing edge information from panel coordinates
 *  @param x Panel node coordinates [2 x N]
 *  @return TEInfo structure with gap, tangent, and strength factors
 */
TEInfo GetTEInfo(const MatrixXd& x);

/** @brief Deploy a flap at hinge location with given deflection angle
 *  @param state    FoilState containing airfoil geometry
 *  @param xzhinge  Hinge location [x, z]
 *  @param eta      Flap deflection angle in degrees (positive = down)
 *
 *  On output:
 *  - state.geom.xpoint contains modified airfoil coordinates
 *  - state.foil is repaneled if panels exist
 */
void GeomFlap(FoilState& state, const Vector2d& xzhinge, double eta);

/** @brief Add camber to airfoil from given camber line increment
 *  @param state   FoilState containing airfoil geometry
 *  @param xzcamb  Camber increment points [2 x Nc] (x, z coords)
 *
 *  On output:
 *  - state.geom.xpoint contains modified airfoil coordinates
 *  - state.foil is repaneled if panels exist
 */
void GeomAddCamber(FoilState& state, const MatrixXd& xzcamb);

/** @brief Derotate airfoil to make chord line horizontal
 *  @param state FoilState containing airfoil geometry
 *
 *  On output:
 *  - state.geom.xpoint contains modified airfoil coordinates
 *  - state.foil is repaneled if panels exist
 */
void GeomDerotate(FoilState& state);

} /* namespace scirocco */
