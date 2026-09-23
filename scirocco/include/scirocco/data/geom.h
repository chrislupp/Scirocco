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

namespace scirocco {

/** @brief Geometry data */
struct Geom {
    double chord = 1.0;   /**< Chord length */
    double wakelen = 1.0; /**< Wake extent length, in chords */
    int npoint = 0;       /**< Number of geometry representation points */
    std::string name = "noname"; /**< Airfoil name, e.g. "NACA 0012" */
    MatrixXd xpoint;             /**< Point coordinates [2 x npoint] */
    Vector2d xref{0.25, 0.0};    /**< Moment reference point */
    double hTE = 0.0;            /**< Trailing edge gap height */
};

}  // namespace scirocco
