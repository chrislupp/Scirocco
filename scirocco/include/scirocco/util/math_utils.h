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
#ifndef SCIROCCO_UTIL_MATH_UTILS_H
#define SCIROCCO_UTIL_MATH_UTILS_H

#include <cmath>

#include "scirocco/config.h"

namespace scirocco {

/**
 * @brief Sine of angle in degrees
 * @param alpha Angle in degrees.
 * @return Sine of alpha.
 */
inline double Sind(double alpha) { return std::sin(alpha * kDegToRad); }

/**
 * @brief Cosine of angle in degrees
 * @param alpha Angle in degrees.
 * @return Cosine of alpha.
 */
inline double Cosd(double alpha) { return std::cos(alpha * kDegToRad); }

/**
 * @brief L2 norm of a vector
 * @param x Input vector.
 * @return Euclidean norm of x.
 */
inline double Norm2(const VectorXd& x) { return x.norm(); }

/**
 * @brief Distance between two scalars (Euclidean)
 * @param a First scalar.
 * @param b Second scalar.
 * @return sqrt(a*a + b*b).
 */
inline double Dist(double a, double b) { return std::sqrt(a * a + b * b); }

/**
 * @brief Two-argument arctangent
 * @param y Numerator (y-component).
 * @param x Denominator (x-component).
 * @return atan2(y, x).
 */
inline double Atan2(double y, double x) { return std::atan2(y, x); }

/**
 * @brief Safe log that prevents log(0)
 * @param x Input value.
 * @return log(max(x, 1e-20)).
 */
inline double SafeLog(double x) { return std::log(std::max(x, 1e-20)); }

/**
 * @brief Clamp value to range [lo, hi]
 * @param val Value to clamp.
 * @param lo  Lower bound.
 * @param hi  Upper bound.
 * @return val clamped to [lo, hi].
 */
inline double Clamp(double val, double lo, double hi) {
    return std::max(lo, std::min(hi, val));
}

/**
 * @brief Sign function: returns -1, 0, or 1
 * @param x Input value.
 * @return -1.0, 0.0, or 1.0 depending on the sign of x.
 */
inline double Sign(double x) {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}

}  // namespace scirocco

#endif  // SCIROCCO_UTIL_MATH_UTILS_H
