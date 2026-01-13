// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_UTIL_MATH_UTILS_H
#define CPPFOIL_UTIL_MATH_UTILS_H

#include <cmath>
#include "cppfoil/config.h"

namespace cppfoil {

// Sine of angle in degrees
inline double Sind(double alpha) {
    return std::sin(alpha * kDegToRad);
}

// Cosine of angle in degrees
inline double Cosd(double alpha) {
    return std::cos(alpha * kDegToRad);
}

// L2 norm of a vector
inline double Norm2(const VectorXd& x) {
    return x.norm();
}

// Distance between two scalars (Euclidean)
inline double Dist(double a, double b) {
    return std::sqrt(a * a + b * b);
}

// Two-argument arctangent
inline double Atan2(double y, double x) {
    return std::atan2(y, x);
}

// Safe log that prevents log(0)
inline double SafeLog(double x) {
    return std::log(std::max(x, 1e-20));
}

// Clamp value to range [lo, hi]
inline double Clamp(double val, double lo, double hi) {
    return std::max(lo, std::min(hi, val));
}

// Sign function: returns -1, 0, or 1
inline double Sign(double x) {
    if (x > 0) return 1.0;
    if (x < 0) return -1.0;
    return 0.0;
}

}  // namespace cppfoil

#endif  // CPPFOIL_UTIL_MATH_UTILS_H
