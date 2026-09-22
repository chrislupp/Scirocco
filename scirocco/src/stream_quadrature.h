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
// Stable evaluation of vortex/source streamfunctions away from a panel.
#pragma once
#include <array>
#include <cmath>

#include "scirocco/core/panel_influence_deriv.h"
#include "scirocco/util/math_utils.h"

namespace scirocco::detail {
// In the far field, closed-form first moments subtract large terms to
// recover small panel contributions. Integrate the log-distance (vortex)
// and angle (source) kernels directly. At r_mid > 4*length, Gauss integration
// resolves the smooth integrand and its coordinate derivatives to double
// precision: the midpoint expansion parameter is < 1/8 and Gauss is exact
// through degree
// 23. Endpoint/near-panel singular limits stay on the closed-form path.
inline bool FarStreamPanel(const Matrix2d& Xj, const Vector2d& xi) {
    const Vector2d dx = Xj.col(1) - Xj.col(0);
    const Vector2d r = xi - 0.5 * (Xj.col(0) + Xj.col(1));
    return dx.squaredNorm() > 1e-24 && r.squaredNorm() > 16 * dx.squaredNorm();
}
inline constexpr std::array<double, 12> kStreamPoints{
    0.0092196828766403782, 0.047941371814762601, 0.11504866290284765,
    0.20634102285669126,   0.31608425050090994,  0.43738329574426554,
    0.56261670425573451,   0.68391574949909006,  0.79365897714330869,
    0.88495133709715235,   0.95205862818523745,  0.99078031712335957};
inline constexpr std::array<double, 12> kStreamWeights{
    0.02358766819325564, 0.05346966299765938, 0.08003916427167318,
    0.10158371336153295, 0.11674626826917744, 0.12457352290670148,
    0.12457352290670148, 0.11674626826917744, 0.10158371336153295,
    0.08003916427167318, 0.05346966299765938, 0.02358766819325564};
inline LinVortexStreamDeriv FarVortexStream(const Matrix2d& Xj,
                                            const Vector2d& xi,
                                            bool derivatives = true) {
    LinVortexStreamDeriv out{};
    out.a_Xj.setZero();
    out.b_Xj.setZero();
    out.a_xi.setZero();
    out.b_xi.setZero();
    const Vector2d dx = Xj.col(1) - Xj.col(0);
    const double length = dx.norm();
    const Vector2d tangent = dx / length;
    const double scale = length / (2 * kPi);
    for (size_t k = 0; k < kStreamPoints.size(); ++k) {
        const double t = kStreamPoints[k];
        const Vector2d r = xi - Xj.col(0) - t * dx;
        const double r2 = r.squaredNorm();
        const double value = scale * kStreamWeights[k] * 0.5 * std::log(r2);
        out.a += (1 - t) * value;
        out.b += t * value;
        if (!derivatives) continue;
        const Vector2d g = scale * kStreamWeights[k] * r / r2;
        out.a_xi += (1 - t) * g;
        out.b_xi += t * g;
        out.a_Xj.head<2>() -= (1 - t) * (1 - t) * g;
        out.a_Xj.tail<2>() -= (1 - t) * t * g;
        out.b_Xj.head<2>() -= t * (1 - t) * g;
        out.b_Xj.tail<2>() -= t * t * g;
    }
    if (derivatives) {
        out.a_Xj.head<2>() -= tangent * (out.a / length);
        out.a_Xj.tail<2>() += tangent * (out.a / length);
        out.b_Xj.head<2>() -= tangent * (out.b / length);
        out.b_Xj.tail<2>() += tangent * (out.b / length);
    }
    return out;
}

// Integrate the source angle kernel. linear selects Python's [0,2*pi) branch;
// constant-source coefficients additionally carry its quarter-panel offset.
inline LinSourceStreamDeriv FarSourceStream(const Matrix2d& Xj,
                                            const Vector2d& xi, bool linear,
                                            bool derivatives = true) {
    LinSourceStreamDeriv out{};
    out.a_Xj.setZero();
    out.b_Xj.setZero();
    out.a_xi.setZero();
    out.b_xi.setZero();
    const Vector2d dx = Xj.col(1) - Xj.col(0);
    const double length = dx.norm();
    const Vector2d tangent = dx / length;
    const Vector2d normal(-tangent(1), tangent(0));
    const Vector2d origin = xi - Xj.col(0);
    const double x = origin.dot(tangent), z = origin.dot(normal);
    const double shift =
        linear ? 0.0
               : (std::atan2(z, x) + std::atan2(z, x - length) > kPi ? -0.25
                                                                     : 0.75);
    const double scale = length / (2 * kPi);
    for (size_t k = 0; k < kStreamPoints.size(); ++k) {
        const double t = kStreamPoints[k];
        double theta = std::atan2(z, x - t * length);
        if (linear && theta < 0) theta += 2 * kPi;
        const double value =
            length * kStreamWeights[k] * (theta / (2 * kPi) + shift);
        out.a += (1 - t) * value;
        out.b += t * value;
        if (!derivatives) continue;
        const Vector2d r = origin - t * dx;
        const Vector2d g =
            scale * kStreamWeights[k] * Vector2d(-r(1), r(0)) / r.squaredNorm();
        // The angle is measured in the moving panel frame: subtract its
        // rotation in addition to differentiating the evaluation-point ray.
        const Vector2d rotation = kStreamWeights[k] * normal / (2 * kPi);
        out.a_xi += (1 - t) * g;
        out.b_xi += t * g;
        out.a_Xj.head<2>() += (1 - t) * (-(1 - t) * g + rotation);
        out.a_Xj.tail<2>() += (1 - t) * (-t * g - rotation);
        out.b_Xj.head<2>() += t * (-(1 - t) * g + rotation);
        out.b_Xj.tail<2>() += t * (-t * g - rotation);
    }
    if (derivatives) {
        out.a_Xj.head<2>() -= tangent * (out.a / length);
        out.a_Xj.tail<2>() += tangent * (out.a / length);
        out.b_Xj.head<2>() -= tangent * (out.b / length);
        out.b_Xj.tail<2>() += tangent * (out.b / length);
    }
    return out;
}
inline bool FarConstantSourcePanel(const Matrix2d& Xj, const Vector2d& xi) {
    if (!FarStreamPanel(Xj, xi)) return false;
    const auto info = GetPanelInfo(Xj, xi);
    // Preserve the original evaluation on its arbitrary streamfunction cut.
    return std::abs(info.theta1 + info.theta2 - kPi) > 1e-12;
}
}  // namespace scirocco::detail
