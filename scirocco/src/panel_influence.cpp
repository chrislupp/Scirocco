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
#include "scirocco/core/panel_influence.h"

#include <cmath>

#include "scirocco/util/math_utils.h"
#include "stream_quadrature.h"

namespace scirocco {

namespace {
constexpr double kEps = 1e-12;
}

PanelInfo GetPanelInfo(const Matrix2d& Xj, const Vector2d& xi) {
    PanelInfo info;

    Vector2d x1 = Xj.col(0);
    Vector2d x2 = Xj.col(1);
    Vector2d dx = x2 - x1;
    info.len = dx.norm();

    /* Panel tangent and normal */
    Vector2d t = dx / info.len;
    Vector2d n(-t(1), t(0));

    /* Transform xi to panel coordinates */
    Vector2d dxi1 = xi - x1;
    Vector2d dxi2 = xi - x2;

    info.x1 = dxi1.dot(t);
    info.x2 = dxi2.dot(t);
    info.y = dxi1.dot(n);

    /* Distances */
    info.r1 = dxi1.norm();
    info.r2 = dxi2.norm();

    /* Angles (atan2 of local coordinates) */
    info.theta1 = Atan2(info.y, info.x1);
    info.theta2 = Atan2(info.y, info.x2);

    return info;
}

void PanelLinVortexVelocity(const Matrix2d& Xj, const Vector2d& xi,
                            const Vector2d* vdir, bool onmid, Vector2d& a,
                            Vector2d& b) {
    /* Linear vortex panel: gamma varies linearly from ga at x1 to gb at x2
     * Velocity at xi due to this distribution */

    Vector2d x1 = Xj.col(0);
    Vector2d x2 = Xj.col(1);
    Vector2d dx = x2 - x1;
    double len = dx.norm();

    if (len < kEps) {
        a.setZero();
        b.setZero();
        return;
    }

    /* Panel tangent and normal */
    Vector2d t = dx / len;
    Vector2d n(-t(1), t(0));

    PanelInfo info = GetPanelInfo(Xj, xi);

    double dtheta = info.theta2 - info.theta1;
    double logr = SafeLog(info.r2 / std::max(info.r1, kEps));

    /* Velocity in panel coordinates
     * u = tangential, v = normal */
    double u_a, v_a, u_b, v_b;

    if (onmid) {
        /* Limit values at panel midpoint (matching mfoil) */
        u_a = 0.25;
        v_a = -1.0 / (2.0 * kPi);
        u_b = 0.25;
        v_b = 1.0 / (2.0 * kPi);
    } else {
        double inv_2pi = 1.0 / (2.0 * kPi);
        double inv_len = 1.0 / len;

        /* Linear vortex velocity (matching mfoil formulas)
         * temp1 = dtheta/(2*pi)
         * temp2 = (2*z*log(r1/r2) - 2*x*dtheta)/(4*pi*d) = (-z*logr -
         * x*dtheta)/(2*pi*d) ug1 = temp1 + temp2 = (d*dtheta - z*logr -
         * x*dtheta)/(2*pi*d)
         *                     = ((d-x)*dtheta - z*logr)/(2*pi*d)
         * ug2 = -temp2 = (z*logr + x*dtheta)/(2*pi*d) */
        double temp1 = dtheta * inv_2pi;
        double temp2 = (-info.y * logr - info.x1 * dtheta) * inv_2pi * inv_len;
        u_a = temp1 + temp2;
        u_b = -temp2;

        /* temp1 = logr/(2*pi)
         * temp2 = (x*log(r1/r2) - d + z*dtheta)/(2*pi*d) = (-x*logr - d +
         * z*dtheta)/(2*pi*d) wg1 = temp1 + temp2 = (d*logr - x*logr - d +
         * z*dtheta)/(2*pi*d)
         *                     = ((d-x)*logr - d + z*dtheta)/(2*pi*d)
         * wg2 = -temp2 = (x*logr + d - z*dtheta)/(2*pi*d) */
        temp1 = logr * inv_2pi;
        temp2 = (-info.x1 * logr - len + info.y * dtheta) * inv_2pi * inv_len;
        v_a = temp1 + temp2;
        v_b = -temp2;
    }

    /* Transform back to global coordinates
     * u is along t, v is along n */
    a = u_a * t + v_a * n;
    b = u_b * t + v_b * n;

    /* Project onto direction if requested */
    if (vdir != nullptr) {
        a = Vector2d(a.dot(*vdir), 0.0);
        b = Vector2d(b.dot(*vdir), 0.0);
    }
}

void PanelLinVortexStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b) {
    if (detail::FarStreamPanel(Xj, xi)) {
        const auto result = detail::FarVortexStream(Xj, xi, false);
        a = result.a;
        b = result.b;
        return;
    }

    /* Streamfunction influence from linear vortex panel
     * Matches Python: P1, P2 formulas with a = P1 - P2/d, b = P2/d */

    PanelInfo info = GetPanelInfo(Xj, xi);

    if (info.len < kEps) {
        a = 0.0;
        b = 0.0;
        return;
    }

    /* Variables matching Python notation:
     * d = panel length, x = local x coord, z = local y coord */
    double d = info.len;
    double x = info.x1;
    double z = info.y;
    double dtheta = info.theta2 - info.theta1;
    double logr1 = (info.r1 < 1e-9) ? 0.0 : std::log(info.r1);
    double logr2 = (info.r2 < 1e-9) ? 0.0 : std::log(info.r2);

    double inv_2pi = 0.5 / kPi;

    /* P1 = (0.5/pi)*(z*(theta2-theta1) - d + x*logr1 - (x-d)*logr2) */
    double P1 = inv_2pi * (z * dtheta - d + x * logr1 - (x - d) * logr2);

    /* P2 = x*P1 + (0.5/pi)*(0.5*r2^2*logr2 - 0.5*r1^2*logr1 - r2^2/4 + r1^2/4)
     */
    double r1_sq = info.r1 * info.r1;
    double r2_sq = info.r2 * info.r2;
    double P2 = x * P1 + inv_2pi * (0.5 * r2_sq * logr2 - 0.5 * r1_sq * logr1 -
                                    r2_sq / 4.0 + r1_sq / 4.0);

    /* a = P1 - P2/d, b = P2/d */
    a = P1 - P2 / d;
    b = P2 / d;
}

Vector2d PanelConstSourceVelocity(const Matrix2d& Xj, const Vector2d& xi,
                                  const Vector2d* vdir) {
    /* Velocity from constant source panel */

    PanelInfo info = GetPanelInfo(Xj, xi);

    if (info.len < kEps) {
        return Vector2d::Zero();
    }

    Vector2d dx = Xj.col(1) - Xj.col(0);
    Vector2d t = dx / info.len;
    Vector2d n(-t(1), t(0));

    double dtheta = info.theta2 - info.theta1;
    double logr = SafeLog(info.r2 / std::max(info.r1, kEps));

    double inv_2pi = 1.0 / (2.0 * kPi);

    /* Velocity in panel coordinates
     * Note: u = (log(r1) - log(r2)) / (2*pi) = -log(r2/r1) / (2*pi) */
    double u = -inv_2pi * logr;
    double v = inv_2pi * dtheta;

    /* Transform to global */
    Vector2d V = u * t + v * n;

    if (vdir != nullptr) {
        return {V.dot(*vdir), 0.0};
    }
    return V;
}

double PanelConstSourceStream(const Matrix2d& Xj, const Vector2d& xi) {
    if (detail::FarConstantSourcePanel(Xj, xi)) {
        const auto q = detail::FarSourceStream(Xj, xi, false, false);
        return q.a + q.b;
    }

    /* Streamfunction from constant source panel
     * Matches Python implementation */

    PanelInfo info = GetPanelInfo(Xj, xi);

    if (info.len < kEps) {
        return 0.0;
    }

    /* Variables matching Python notation */
    double d = info.len;
    double x = info.x1;
    double z = info.y;
    double theta1 = info.theta1;
    double theta2 = info.theta2;
    double r1 = info.r1;
    double r2 = info.r2;

    /* Handle singular cases */
    double logr1, logr2;
    double ep = 1e-9;
    if (r1 < ep) {
        logr1 = 0.0;
        theta1 = kPi;
        theta2 = kPi;
    } else {
        logr1 = std::log(r1);
    }
    if (r2 < ep) {
        logr2 = 0.0;
        theta1 = 0.0;
        theta2 = 0.0;
    } else {
        logr2 = std::log(r2);
    }

    /* P = (x*(theta1-theta2) + d*theta2 + z*logr1 - z*logr2)/(2*pi) */
    double P = (x * (theta1 - theta2) + d * theta2 + z * logr1 - z * logr2) /
               (2.0 * kPi);

    /* dP = d; adjust P based on theta sum */
    double dP = d;
    if ((theta1 + theta2) > kPi) {
        P = P - 0.25 * dP;
    } else {
        P = P + 0.75 * dP;
    }

    return P;
}

void PanelLinSourceStream(const Matrix2d& Xj, const Vector2d& xi, double& a,
                          double& b) {
    if (detail::FarStreamPanel(Xj, xi)) {
        const auto q = detail::FarSourceStream(Xj, xi, true, false);
        a = q.a;
        b = q.b;
        return;
    }

    /* Streamfunction from linear source panel
     * Matches mfoil panel_linsource_stream */

    PanelInfo info = GetPanelInfo(Xj, xi);

    if (info.len < kEps) {
        a = 0.0;
        b = 0.0;
        return;
    }

    /* Variables matching Python notation */
    double x = info.x1;
    double z = info.y;
    double d = info.len;
    double r1 = info.r1;
    double r2 = info.r2;
    double theta1 = info.theta1;
    double theta2 = info.theta2;

    /* Make branch cut at theta = 0 (mfoil does this) */
    if (theta1 < 0) theta1 += 2.0 * kPi;
    if (theta2 < 0) theta2 += 2.0 * kPi;

    /* Check for r1, r2 zero */
    double ep = 1e-9;
    double logr1, logr2;
    if (r1 < ep) {
        logr1 = 0.0;
        theta1 = kPi;
        theta2 = kPi;
    } else {
        logr1 = std::log(r1);
    }
    if (r2 < ep) {
        logr2 = 0.0;
        theta1 = 0.0;
        theta2 = 0.0;
    } else {
        logr2 = std::log(r2);
    }

    double inv_2pi = 0.5 / kPi;

    /* Streamfunction components (mfoil formula) */
    double P1 =
        inv_2pi * (x * (theta1 - theta2) + theta2 * d + z * logr1 - z * logr2);
    double P2 = x * P1 + inv_2pi * (0.5 * r2 * r2 * theta2 -
                                    0.5 * r1 * r1 * theta1 - 0.5 * z * d);

    /* Influence coefficients */
    a = P1 - P2 / d;
    b = P2 / d;
}

void PanelLinSourceVelocity(const Matrix2d& Xj, const Vector2d& xi, Vector2d& a,
                            Vector2d& b) {
    /* Velocity from linear source panel
     * Matches mfoil panel_linsource_velocity */

    PanelInfo info = GetPanelInfo(Xj, xi);

    if (info.len < kEps) {
        a.setZero();
        b.setZero();
        return;
    }

    Vector2d dx = Xj.col(1) - Xj.col(0);
    Vector2d t = dx / info.len;
    Vector2d n(-t(1), t(0));

    /* Variables matching Python notation */
    double x = info.x1;
    double z = info.y;
    double d = info.len;
    double r1 = info.r1;
    double r2 = info.r2;
    double dtheta = info.theta2 - info.theta1;

    double inv_2pi = 0.5 / kPi;
    double inv_d = 1.0 / d;

    /* log(r1/r2) used in mfoil */
    double logr12 = SafeLog(r1) - SafeLog(r2);

    /* Velocity in panel coordinates (mfoil formulas) */
    double temp1 = logr12 * inv_2pi;
    double temp2 = (x * logr12 - d + z * dtheta) * inv_2pi * inv_d;
    double ug1 = temp1 - temp2;
    double ug2 = temp2;

    temp1 = dtheta * inv_2pi;
    temp2 = (-z * logr12 + x * dtheta) * inv_2pi * inv_d;
    double wg1 = temp1 - temp2;
    double wg2 = temp2;

    /* Transform to global coordinates */
    a = ug1 * t + wg1 * n;
    b = ug2 * t + wg2 * n;
}

} /* namespace scirocco */
