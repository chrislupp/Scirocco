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

namespace scirocco {

/** @brief Panel geometry information */
struct PanelInfo {
    double r1, r2;          /* Distances from xi to panel endpoints */
    double theta1, theta2;  /* Angles from xi to panel endpoints */
    double x1, x2;          /* Local x coordinates (along panel) */
    double y;               /* Local y coordinate (perpendicular to panel) */
    double len;             /* Panel length */
};

/**
 * @brief Compute panel geometry information
 * @param Xj Panel endpoint coordinates [2 x 2] (j, j+1)
 * @param xi Point at which to evaluate
 * @return PanelInfo structure with geometric quantities
 */
PanelInfo GetPanelInfo(const Matrix2d& Xj, const Vector2d& xi);

/**
 * @brief Linear vortex panel velocity influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @param vdir Optional velocity direction for projection (nullptr for both
 *   components)
 * @param onmid True if xi is at panel midpoint (use limit values)
 * @param a Influence coefficient for left vortex strength
 * @param b Influence coefficient for right vortex strength
 */
void PanelLinVortexVelocity(const Matrix2d& Xj, const Vector2d& xi,
                            const Vector2d* vdir, bool onmid, Vector2d& a,
                            Vector2d& b);

/**
 * @brief Linear vortex panel streamfunction influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @param a Influence coefficient for left vortex strength
 * @param b Influence coefficient for right vortex strength
 */
void PanelLinVortexStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b);

/**
 * @brief Constant source panel velocity influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @param vdir Optional velocity direction for projection (nullptr for both
 *   components)
 * @return Velocity influence (or scalar if vdir provided)
 */
Vector2d PanelConstSourceVelocity(const Matrix2d& Xj, const Vector2d& xi,
                                  const Vector2d* vdir);

/**
 * @brief Constant source panel streamfunction influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @return Streamfunction influence coefficient
 */
double PanelConstSourceStream(const Matrix2d& Xj, const Vector2d& xi);

/**
 * @brief Linear source panel streamfunction influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @param a Influence coefficient for left source strength
 * @param b Influence coefficient for right source strength
 */
void PanelLinSourceStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b);

/**
 * @brief Linear source panel velocity influence
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @param a Velocity influence vector for left source strength
 * @param b Velocity influence vector for right source strength
 */
void PanelLinSourceVelocity(const Matrix2d& Xj, const Vector2d& xi, Vector2d& a,
                            Vector2d& b);

}  /* namespace scirocco */
