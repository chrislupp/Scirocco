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
#include "scirocco/core/panel_influence_deriv.h"

#include <cmath>

#include "scirocco/util/math_utils.h"
#include "stream_quadrature.h"

namespace scirocco {

namespace {
constexpr double kEps = 1e-12;
}

/**
 * @brief Differentiate GetPanelInfo w.r.t. Xj and xi
 *
 * Xj = [x1_x, x1_z; x2_x, x2_z] (2x2 matrix, columns are endpoints)
 * xi = [xi_x, xi_z]
 *
 * Derivative ordering for Xj: [x1_x, x1_z, x2_x, x2_z] (column-major of Xj)
 * @param Xj Panel endpoint coordinates [2 x 2]
 * @param xi Evaluation point
 * @return PanelInfoDeriv containing primal info and all first derivatives
 */
PanelInfoDeriv GetPanelInfoDeriv(const Matrix2d& Xj, const Vector2d& xi) {
    PanelInfoDeriv d;
    d.info = GetPanelInfo(Xj, xi);

    Vector2d x1p = Xj.col(0); /* panel endpoint 1 */
    Vector2d x2p = Xj.col(1); /* panel endpoint 2 */
    Vector2d dx = x2p - x1p;
    double len = d.info.len;

    if (len < kEps) {
        /* Zero-length panel: zero all derivatives */
        d.x1_Xj.setZero();
        d.x2_Xj.setZero();
        d.y_Xj.setZero();
        d.x1_xi.setZero();
        d.x2_xi.setZero();
        d.y_xi.setZero();
        d.r1_Xj.setZero();
        d.r2_Xj.setZero();
        d.r1_xi.setZero();
        d.r2_xi.setZero();
        d.theta1_Xj.setZero();
        d.theta2_Xj.setZero();
        d.theta1_xi.setZero();
        d.theta2_xi.setZero();
        d.len_Xj.setZero();
        return d;
    }

    double inv_len = 1.0 / len;

    /* Panel tangent t = dx / len and normal n = [-t(1), t(0)] */
    Vector2d t = dx * inv_len;
    Vector2d n(-t(1), t(0));

    /* ---- Derivatives of t and n w.r.t. Xj ----
     * t = (x2p - x1p) / len
     * dt/d(x1p) = (-I*len - (x2p-x1p) * d(len)/d(x1p)) / len^2
     * d(len)/d(x1p) = -(x2p-x1p)^T / len = -t^T
     * d(len)/d(x2p) = (x2p-x1p)^T / len = t^T */

    /* d(len)/d(Xj) = [-t(0), -t(1), t(0), t(1)] */
    d.len_Xj << -t(0), -t(1), t(0), t(1);

    /* dt/d(Xj): 2x4 matrix, but we only need the chain rule products below
     * dt/d(x1p_x) = (-e_x * len + dx * t(0)) / len^2 = (-e_x + t*t(0)) / len
     * dt/d(x1p_z) = (-e_z * len + dx * t(1)) / len^2 = (-e_z + t*t(1)) / len
     * dt/d(x2p_x) = (e_x * len - dx * t(0)) / len^2 = (e_x - t*t(0)) / len
     * dt/d(x2p_z) = (e_z * len - dx * t(1)) / len^2 = (e_z - t*t(1)) / len */

    /* Store dt/d(Xj) as 2x4 matrix for chain rule */
    Eigen::Matrix<double, 2, 4> dt_dXj;
    dt_dXj.col(0) = (-Vector2d::Unit(0) + t * t(0)) * inv_len; /* d/d(x1p_x) */
    dt_dXj.col(1) = (-Vector2d::Unit(1) + t * t(1)) * inv_len; /* d/d(x1p_z) */
    dt_dXj.col(2) = (Vector2d::Unit(0) - t * t(0)) * inv_len;  /* d/d(x2p_x) */
    dt_dXj.col(3) = (Vector2d::Unit(1) - t * t(1)) * inv_len;  /* d/d(x2p_z) */

    /* dn/d(Xj): n = [-t(1), t(0)] */
    Eigen::Matrix<double, 2, 4> dn_dXj;
    for (int k = 0; k < 4; ++k) {
        dn_dXj(0, k) = -dt_dXj(1, k);
        dn_dXj(1, k) = dt_dXj(0, k);
    }

    /* ---- Local coordinates ----
     * dxi1 = xi - x1p, dxi2 = xi - x2p */
    Vector2d dxi1 = xi - x1p;
    Vector2d dxi2 = xi - x2p;

    /* x1_local = dxi1.dot(t), x2_local = dxi2.dot(t), y_local = dxi1.dot(n) */

    /* d(x1_local)/d(Xj_k) = d(dxi1)/d(Xj_k) . t + dxi1 . d(t)/d(Xj_k)
     * d(dxi1)/d(x1p_x) = -e_x, d(dxi1)/d(x1p_z) = -e_z
     * d(dxi1)/d(x2p_x) = 0, d(dxi1)/d(x2p_z) = 0 */

    for (int k = 0; k < 4; ++k) {
        Vector2d ddxi1_dXjk = Vector2d::Zero();
        Vector2d ddxi2_dXjk = Vector2d::Zero();
        if (k == 0) {
            ddxi1_dXjk(0) = -1.0;
        }
        if (k == 1) {
            ddxi1_dXjk(1) = -1.0;
        }
        if (k == 2) {
            ddxi2_dXjk(0) = -1.0;
        }
        if (k == 3) {
            ddxi2_dXjk(1) = -1.0;
        }

        d.x1_Xj(k) = ddxi1_dXjk.dot(t) + dxi1.dot(dt_dXj.col(k));
        d.x2_Xj(k) = ddxi2_dXjk.dot(t) + dxi2.dot(dt_dXj.col(k));
        d.y_Xj(k) = ddxi1_dXjk.dot(n) + dxi1.dot(dn_dXj.col(k));
    }

    /* d/d(xi): d(dxi1)/d(xi) = I, d(dxi2)/d(xi) = I
     * d(t)/d(xi) = 0, d(n)/d(xi) = 0 */
    d.x1_xi = t;
    d.x2_xi = t;
    d.y_xi = n;

    /* ---- Distances r1, r2 ----
     * r1 = ||dxi1||, r2 = ||dxi2|| */
    double r1 = d.info.r1;
    double r2 = d.info.r2;

    /* d(r1)/d(Xj_k) = dxi1 . d(dxi1)/d(Xj_k) / r1 */
    if (r1 > kEps) {
        double inv_r1 = 1.0 / r1;
        d.r1_Xj(0) = -dxi1(0) * inv_r1;
        d.r1_Xj(1) = -dxi1(1) * inv_r1;
        d.r1_Xj(2) = 0.0;
        d.r1_Xj(3) = 0.0;
        d.r1_xi = dxi1 * inv_r1;
    } else {
        d.r1_Xj.setZero();
        d.r1_xi.setZero();
    }

    if (r2 > kEps) {
        double inv_r2 = 1.0 / r2;
        d.r2_Xj(0) = 0.0;
        d.r2_Xj(1) = 0.0;
        d.r2_Xj(2) = -dxi2(0) * inv_r2;
        d.r2_Xj(3) = -dxi2(1) * inv_r2;
        d.r2_xi = dxi2 * inv_r2;
    } else {
        d.r2_Xj.setZero();
        d.r2_xi.setZero();
    }

    /* ---- Angles theta1, theta2 ----
     * theta1 = atan2(y, x1), theta2 = atan2(y, x2)
     * d(atan2(y,x))/dy = x/(x^2+y^2), d/dx = -y/(x^2+y^2) */
    double y_val = d.info.y;
    double x1_val = d.info.x1;
    double x2_val = d.info.x2;

    double r1sq = r1 * r1;
    double r2sq = r2 * r2;

    if (r1sq > kEps * kEps) {
        double inv_r1sq = 1.0 / r1sq;
        /* dtheta1/d(Xj) = (x1_val / r1sq) * dy/d(Xj) + (-y_val / r1sq) *
         * dx1/d(Xj) */
        d.theta1_Xj =
            (x1_val * inv_r1sq) * d.y_Xj + (-y_val * inv_r1sq) * d.x1_Xj;
        d.theta1_xi =
            (x1_val * inv_r1sq) * d.y_xi + (-y_val * inv_r1sq) * d.x1_xi;
    } else {
        d.theta1_Xj.setZero();
        d.theta1_xi.setZero();
    }

    if (r2sq > kEps * kEps) {
        double inv_r2sq = 1.0 / r2sq;
        d.theta2_Xj =
            (x2_val * inv_r2sq) * d.y_Xj + (-y_val * inv_r2sq) * d.x2_Xj;
        d.theta2_xi =
            (x2_val * inv_r2sq) * d.y_xi + (-y_val * inv_r2sq) * d.x2_xi;
    } else {
        d.theta2_Xj.setZero();
        d.theta2_xi.setZero();
    }

    return d;
}

/* ============================================================================
 * PanelLinVortexStreamDeriv
 * ============================================================================
 */
LinVortexStreamDeriv PanelLinVortexStreamDeriv(const Matrix2d& Xj,
                                               const Vector2d& xi) {
    if (detail::FarStreamPanel(Xj, xi)) return detail::FarVortexStream(Xj, xi);

    LinVortexStreamDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.a = 0.0;
        result.b = 0.0;
        result.a_Xj.setZero();
        result.b_Xj.setZero();
        result.a_xi.setZero();
        result.b_xi.setZero();
        return result;
    }

    double d = info.len;
    double x = info.x1;
    double z = info.y;
    double dtheta = info.theta2 - info.theta1;
    double logr1 = (info.r1 < 1e-9) ? 0.0 : std::log(info.r1);
    double logr2 = (info.r2 < 1e-9) ? 0.0 : std::log(info.r2);

    double inv_2pi = 0.5 / kPi;

    /* P1 = inv_2pi * (z*dtheta - d + x*logr1 - (x-d)*logr2) */
    double P1 = inv_2pi * (z * dtheta - d + x * logr1 - (x - d) * logr2);

    /* P2 = x*P1 + inv_2pi * (0.5*r2^2*logr2 - 0.5*r1^2*logr1 - r2^2/4 + r1^2/4)
     */
    double r1_sq = info.r1 * info.r1;
    double r2_sq = info.r2 * info.r2;
    double P2 = x * P1 + inv_2pi * (0.5 * r2_sq * logr2 - 0.5 * r1_sq * logr1 -
                                    r2_sq / 4.0 + r1_sq / 4.0);

    result.a = P1 - P2 / d;
    result.b = P2 / d;

    /* Now differentiate P1 and P2
     * Need derivatives of: x, z, d, dtheta, logr1, logr2, r1^2, r2^2 */

    /* dlogr1/d(...) = (1/r1) * dr1/d(...) if r1 > eps
     * dlogr2/d(...) = (1/r2) * dr2/d(...) if r2 > eps */

    auto diff = [&](const Vector4d& x_d, const Vector4d& z_d,
                    const Vector4d& d_d, const Vector4d& dtheta_d,
                    const Vector4d& logr1_d, const Vector4d& logr2_d,
                    const Vector4d& r1sq_d, const Vector4d& r2sq_d,
                    Vector4d& a_d, Vector4d& b_d) {
        /* dP1 = inv_2pi * (dz*dtheta + z*d(dtheta) - dd + dx*logr1 + x*dlogr1
         *        - (dx - dd)*logr2 - (x-d)*dlogr2) */
        Vector4d dP1 =
            inv_2pi * (z_d * dtheta + z * dtheta_d - d_d + x_d * logr1 +
                       x * logr1_d - (x_d - d_d) * logr2 - (x - d) * logr2_d);

        /* dP2 = dx*P1 + x*dP1 + inv_2pi*(0.5*dr2sq*logr2 + 0.5*r2sq*dlogr2
         *        - 0.5*dr1sq*logr1 - 0.5*r1sq*dlogr1 - dr2sq/4 + dr1sq/4) */
        Vector4d dP2 = x_d * P1 + x * dP1 +
                       inv_2pi * (0.5 * r2sq_d * logr2 + 0.5 * r2_sq * logr2_d -
                                  0.5 * r1sq_d * logr1 - 0.5 * r1_sq * logr1_d -
                                  r2sq_d / 4.0 + r1sq_d / 4.0);

        /* a = P1 - P2/d, b = P2/d
         * da = dP1 - (dP2*d - P2*dd) / d^2
         * db = (dP2*d - P2*dd) / d^2 */
        double inv_d = 1.0 / d;
        double inv_d2 = inv_d * inv_d;
        a_d = dP1 - (dP2 * d - P2 * d_d) * inv_d2;
        b_d = (dP2 * d - P2 * d_d) * inv_d2;
    };

    /* Build derivative vectors for Xj */
    Vector4d dtheta_Xj = pid.theta2_Xj - pid.theta1_Xj;
    Vector4d logr1_Xj =
        (info.r1 > 1e-9) ? (pid.r1_Xj / info.r1).eval() : Vector4d::Zero();
    Vector4d logr2_Xj =
        (info.r2 > 1e-9) ? (pid.r2_Xj / info.r2).eval() : Vector4d::Zero();
    Vector4d r1sq_Xj = 2.0 * info.r1 * pid.r1_Xj;
    Vector4d r2sq_Xj = 2.0 * info.r2 * pid.r2_Xj;

    diff(pid.x1_Xj, pid.y_Xj, pid.len_Xj, dtheta_Xj, logr1_Xj, logr2_Xj,
         r1sq_Xj, r2sq_Xj, result.a_Xj, result.b_Xj);

    /* Build derivative vectors for xi (len doesn't depend on xi) */
    Vector2d dtheta_xi = pid.theta2_xi - pid.theta1_xi;
    Vector2d logr1_xi =
        (info.r1 > 1e-9) ? (pid.r1_xi / info.r1).eval() : Vector2d::Zero();
    Vector2d logr2_xi =
        (info.r2 > 1e-9) ? (pid.r2_xi / info.r2).eval() : Vector2d::Zero();
    Vector2d r1sq_xi = 2.0 * info.r1 * pid.r1_xi;
    Vector2d r2sq_xi = 2.0 * info.r2 * pid.r2_xi;
    Vector2d len_xi = Vector2d::Zero();

    /* Reuse the diff lambda with Vector2d */
    {
        Vector2d dP1 =
            inv_2pi *
            (pid.y_xi * dtheta + z * dtheta_xi - len_xi + pid.x1_xi * logr1 +
             x * logr1_xi - (pid.x1_xi - len_xi) * logr2 - (x - d) * logr2_xi);

        Vector2d dP2 =
            pid.x1_xi * P1 + x * dP1 +
            inv_2pi * (0.5 * r2sq_xi * logr2 + 0.5 * r2_sq * logr2_xi -
                       0.5 * r1sq_xi * logr1 - 0.5 * r1_sq * logr1_xi -
                       r2sq_xi / 4.0 + r1sq_xi / 4.0);

        double inv_d_xi = 1.0 / d;
        result.a_xi = dP1 - dP2 * inv_d_xi;
        result.b_xi = dP2 * inv_d_xi;
    }

    return result;
}

/* ============================================================================
 * PanelConstSourceStreamDeriv
 * ============================================================================
 */
ConstSourceStreamDeriv PanelConstSourceStreamDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi) {
    if (detail::FarConstantSourcePanel(Xj, xi)) {
        const auto q = detail::FarSourceStream(Xj, xi, false);
        return {q.a + q.b, q.a_Xj + q.b_Xj, q.a_xi + q.b_xi};
    }

    ConstSourceStreamDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.val = 0.0;
        result.val_Xj.setZero();
        result.val_xi.setZero();
        return result;
    }

    double d = info.len;
    double x = info.x1;
    double z = info.y;
    double theta1 = info.theta1;
    double theta2 = info.theta2;
    double r1 = info.r1;
    double r2 = info.r2;

    double logr1, logr2;
    double ep = 1e-9;
    bool r1_zero = r1 < ep;
    bool r2_zero = r2 < ep;

    /* For derivative purposes, when r1 or r2 is zero the thetas are
     * replaced with constants and logr is set to 0. We treat these as
     * non-differentiable points and set the corresponding derivatives to zero.
     */
    if (r1_zero) {
        logr1 = 0.0;
        theta1 = kPi;
        theta2 = kPi;
    } else {
        logr1 = std::log(r1);
    }
    if (r2_zero) {
        logr2 = 0.0;
        theta1 = 0.0;
        theta2 = 0.0;
    } else {
        logr2 = std::log(r2);
    }

    /* P = (x*(theta1-theta2) + d*theta2 + z*logr1 - z*logr2) / (2*pi) */
    double P = (x * (theta1 - theta2) + d * theta2 + z * logr1 - z * logr2) /
               (2.0 * kPi);

    double dP_val = d;
    if ((theta1 + theta2) > kPi) {
        P = P - 0.25 * dP_val;
    } else {
        P = P + 0.75 * dP_val;
    }
    result.val = P;

    /* Differentiate P (not the branch-based offset which only depends on d)
     * Note: we treat the branch as fixed (not differentiating through the
     * branch choice) */
    double offset_sign = ((theta1 + theta2) > kPi) ? -0.25 : 0.75;

    /* When r is near zero, theta and logr derivatives are zero */
    Vector4d th1_Xj = r1_zero ? Vector4d::Zero() : pid.theta1_Xj;
    Vector4d th2_Xj = r2_zero ? Vector4d::Zero() : pid.theta2_Xj;
    Vector4d logr1_Xj_d =
        (r1 > ep && !r1_zero) ? (pid.r1_Xj / r1).eval() : Vector4d::Zero();
    Vector4d logr2_Xj_d =
        (r2 > ep && !r2_zero) ? (pid.r2_Xj / r2).eval() : Vector4d::Zero();

    result.val_Xj = (1.0 / (2.0 * kPi)) *
                        (pid.x1_Xj * (theta1 - theta2) + x * (th1_Xj - th2_Xj) +
                         pid.len_Xj * theta2 + d * th2_Xj + pid.y_Xj * logr1 +
                         z * logr1_Xj_d - pid.y_Xj * logr2 - z * logr2_Xj_d) +
                    offset_sign * pid.len_Xj;

    /* xi derivatives */
    Vector2d th1_xi = r1_zero ? Vector2d::Zero() : pid.theta1_xi;
    Vector2d th2_xi = r2_zero ? Vector2d::Zero() : pid.theta2_xi;
    Vector2d logr1_xi_d =
        (r1 > ep && !r1_zero) ? (pid.r1_xi / r1).eval() : Vector2d::Zero();
    Vector2d logr2_xi_d =
        (r2 > ep && !r2_zero) ? (pid.r2_xi / r2).eval() : Vector2d::Zero();

    result.val_xi =
        (1.0 / (2.0 * kPi)) *
        (pid.x1_xi * (theta1 - theta2) + x * (th1_xi - th2_xi) + d * th2_xi +
         pid.y_xi * logr1 + z * logr1_xi_d - pid.y_xi * logr2 - z * logr2_xi_d);

    return result;
}

/* ============================================================================
 * PanelLinVortexVelocityDeriv
 * ============================================================================
 */
LinVortexVelocityDeriv PanelLinVortexVelocityDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi) {
    LinVortexVelocityDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.a.setZero();
        result.b.setZero();
        result.a_Xj.setZero();
        result.b_Xj.setZero();
        result.a_xi.setZero();
        result.b_xi.setZero();
        return result;
    }

    double len = info.len;
    double inv_len = 1.0 / len;
    double inv_2pi = 1.0 / (2.0 * kPi);
    double dtheta = info.theta2 - info.theta1;
    double logr = SafeLog(info.r2 / std::max(info.r1, kEps));

    double x = info.x1;
    double z = info.y;

    /* Panel tangent and normal in global coords */
    Vector2d dx_vec = Xj.col(1) - Xj.col(0);
    Vector2d t = dx_vec / len;
    Vector2d n(-t(1), t(0));

    /* Velocity coefficients in panel coords (non-onmid case) */
    double temp1_u = dtheta * inv_2pi;
    double temp2_u = (-z * logr - x * dtheta) * inv_2pi * inv_len;
    double u_a = temp1_u + temp2_u;
    double u_b = -temp2_u;

    double temp1_v = logr * inv_2pi;
    double temp2_v = (-x * logr - len + z * dtheta) * inv_2pi * inv_len;
    double v_a = temp1_v + temp2_v;
    double v_b = -temp2_v;

    /* Global velocity: a = u_a*t + v_a*n, b = u_b*t + v_b*n */
    result.a = u_a * t + v_a * n;
    result.b = u_b * t + v_b * n;

    /* dt/dXj (2x4), dn/dXj (2x4) */
    Eigen::Matrix<double, 2, 4> dt_dXj;
    dt_dXj.col(0) = (-Vector2d::Unit(0) + t * t(0)) * inv_len;
    dt_dXj.col(1) = (-Vector2d::Unit(1) + t * t(1)) * inv_len;
    dt_dXj.col(2) = (Vector2d::Unit(0) - t * t(0)) * inv_len;
    dt_dXj.col(3) = (Vector2d::Unit(1) - t * t(1)) * inv_len;
    Eigen::Matrix<double, 2, 4> dn_dXj;
    for (int k = 0; k < 4; ++k) {
        dn_dXj(0, k) = -dt_dXj(1, k);
        dn_dXj(1, k) = dt_dXj(0, k);
    }

    /* logr derivatives
     * logr = log(r2) - log(r1) (when r1 > eps)
     * dlogr/d(...) = (1/r2)*dr2/d(...) - (1/r1)*dr1/d(...) */
    Vector4d logr_Xj = Vector4d::Zero();
    Vector2d logr_xi = Vector2d::Zero();
    if (info.r1 > kEps) {
        logr_Xj -= pid.r1_Xj / info.r1;
        logr_xi -= pid.r1_xi / info.r1;
    }
    if (info.r2 > kEps) {
        logr_Xj += pid.r2_Xj / info.r2;
        logr_xi += pid.r2_xi / info.r2;
    }

    Vector4d dtheta_Xj = pid.theta2_Xj - pid.theta1_Xj;
    Vector2d dtheta_xi = pid.theta2_xi - pid.theta1_xi;

    /* Differentiate u_a, u_b, v_a, v_b w.r.t. Xj and xi
     * temp2_u = (-z*logr - x*dtheta) * inv_2pi * inv_len
     * d(temp2_u) = inv_2pi * ( d[(-z*logr - x*dtheta)/len] )
     *            = inv_2pi * ( (-dz*logr - z*dlogr - dx*dtheta -
     *            x*d(dtheta))*len
     *                         - (-z*logr - x*dtheta)*dlen ) / len^2 */

    auto diffCoeffs = [&](const auto& x_d, const auto& z_d, const auto& d_d,
                          const auto& dtheta_d, const auto& logr_d, auto& ua_d,
                          auto& ub_d, auto& va_d, auto& vb_d) {
        /* dtemp1_u = d(dtheta)*inv_2pi */
        auto dtemp1_u = dtheta_d * inv_2pi;

        /* dtemp2_u = inv_2pi * [(-dz*logr - z*dlogr - dx*dtheta -
         * x*d(dtheta))/len
         *            - (-z*logr - x*dtheta)*dlen/len^2] */
        double numer = -z * logr - x * dtheta;
        auto dnumer = -z_d * logr - z * logr_d - x_d * dtheta - x * dtheta_d;
        auto dtemp2_u =
            inv_2pi * (dnumer * inv_len - numer * d_d * inv_len * inv_len);

        ua_d = dtemp1_u + dtemp2_u;
        ub_d = -dtemp2_u;

        /* temp1_v = logr * inv_2pi */
        auto dtemp1_v = logr_d * inv_2pi;

        /* temp2_v = (-x*logr - len + z*dtheta) * inv_2pi * inv_len */
        double numer2 = -x * logr - len + z * dtheta;
        auto dnumer2 =
            -x_d * logr - x * logr_d - d_d + z_d * dtheta + z * dtheta_d;
        auto dtemp2_v =
            inv_2pi * (dnumer2 * inv_len - numer2 * d_d * inv_len * inv_len);

        va_d = dtemp1_v + dtemp2_v;
        vb_d = -dtemp2_v;
    };

    /* Xj derivatives of coefficients */
    Vector4d ua_Xj, ub_Xj, va_Xj, vb_Xj;
    diffCoeffs(pid.x1_Xj, pid.y_Xj, pid.len_Xj, dtheta_Xj, logr_Xj, ua_Xj,
               ub_Xj, va_Xj, vb_Xj);

    /* xi derivatives of coefficients */
    Vector2d len_xi = Vector2d::Zero();
    Vector2d ua_xi, ub_xi, va_xi, vb_xi;
    diffCoeffs(pid.x1_xi, pid.y_xi, len_xi, dtheta_xi, logr_xi, ua_xi, ub_xi,
               va_xi, vb_xi);

    /* a = u_a*t + v_a*n => da/d(Xj_k) = dua*t + ua*dt + dva*n + va*dn */
    for (int k = 0; k < 4; ++k) {
        result.a_Xj.col(k) = ua_Xj(k) * t + u_a * dt_dXj.col(k) + va_Xj(k) * n +
                             v_a * dn_dXj.col(k);
        result.b_Xj.col(k) = ub_Xj(k) * t + u_b * dt_dXj.col(k) + vb_Xj(k) * n +
                             v_b * dn_dXj.col(k);
    }

    /* For xi, t and n don't depend on xi */
    for (int k = 0; k < 2; ++k) {
        result.a_xi.col(k) = ua_xi(k) * t + va_xi(k) * n;
        result.b_xi.col(k) = ub_xi(k) * t + vb_xi(k) * n;
    }

    return result;
}

/* ============================================================================
 * PanelConstSourceVelocityDeriv
 * ============================================================================
 */
ConstSourceVelocityDeriv PanelConstSourceVelocityDeriv(const Matrix2d& Xj,
                                                       const Vector2d& xi) {
    ConstSourceVelocityDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.val.setZero();
        result.val_Xj.setZero();
        result.val_xi.setZero();
        return result;
    }

    double inv_2pi = 1.0 / (2.0 * kPi);
    double dtheta = info.theta2 - info.theta1;
    double logr = SafeLog(info.r2 / std::max(info.r1, kEps));

    /* Velocity in panel coords: u = -logr*inv_2pi, v = dtheta*inv_2pi */
    double u = -inv_2pi * logr;
    double v = inv_2pi * dtheta;

    Vector2d dx_vec = Xj.col(1) - Xj.col(0);
    double len = info.len;
    double inv_len = 1.0 / len;
    Vector2d t = dx_vec * inv_len;
    Vector2d n(-t(1), t(0));

    result.val = u * t + v * n;

    /* dt/dXj, dn/dXj */
    Eigen::Matrix<double, 2, 4> dt_dXj;
    dt_dXj.col(0) = (-Vector2d::Unit(0) + t * t(0)) * inv_len;
    dt_dXj.col(1) = (-Vector2d::Unit(1) + t * t(1)) * inv_len;
    dt_dXj.col(2) = (Vector2d::Unit(0) - t * t(0)) * inv_len;
    dt_dXj.col(3) = (Vector2d::Unit(1) - t * t(1)) * inv_len;
    Eigen::Matrix<double, 2, 4> dn_dXj;
    for (int k = 0; k < 4; ++k) {
        dn_dXj(0, k) = -dt_dXj(1, k);
        dn_dXj(1, k) = dt_dXj(0, k);
    }

    /* logr and dtheta derivatives */
    Vector4d logr_Xj = Vector4d::Zero();
    Vector2d logr_xi = Vector2d::Zero();
    if (info.r1 > kEps) {
        logr_Xj -= pid.r1_Xj / info.r1;
        logr_xi -= pid.r1_xi / info.r1;
    }
    if (info.r2 > kEps) {
        logr_Xj += pid.r2_Xj / info.r2;
        logr_xi += pid.r2_xi / info.r2;
    }
    Vector4d dtheta_Xj = pid.theta2_Xj - pid.theta1_Xj;
    Vector2d dtheta_xi = pid.theta2_xi - pid.theta1_xi;

    Vector4d du_Xj = -inv_2pi * logr_Xj;
    Vector4d dv_Xj = inv_2pi * dtheta_Xj;
    Vector2d du_xi = -inv_2pi * logr_xi;
    Vector2d dv_xi = inv_2pi * dtheta_xi;

    for (int k = 0; k < 4; ++k) {
        result.val_Xj.col(k) =
            du_Xj(k) * t + u * dt_dXj.col(k) + dv_Xj(k) * n + v * dn_dXj.col(k);
    }
    for (int k = 0; k < 2; ++k) {
        result.val_xi.col(k) = du_xi(k) * t + dv_xi(k) * n;
    }

    return result;
}

/* ============================================================================
 * PanelLinSourceStreamDeriv
 * ============================================================================
 */
LinSourceStreamDeriv PanelLinSourceStreamDeriv(const Matrix2d& Xj,
                                               const Vector2d& xi) {
    if (detail::FarStreamPanel(Xj, xi))
        return detail::FarSourceStream(Xj, xi, true);

    LinSourceStreamDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.a = 0.0;
        result.b = 0.0;
        result.a_Xj.setZero();
        result.b_Xj.setZero();
        result.a_xi.setZero();
        result.b_xi.setZero();
        return result;
    }

    double x = info.x1;
    double z = info.y;
    double d = info.len;
    double r1 = info.r1;
    double r2 = info.r2;
    double theta1 = info.theta1;
    double theta2 = info.theta2;

    /* Branch cut at theta = 0 */
    if (theta1 < 0) theta1 += 2.0 * kPi;
    if (theta2 < 0) theta2 += 2.0 * kPi;

    double ep = 1e-9;
    bool r1_zero = r1 < ep;
    bool r2_zero = r2 < ep;
    double logr1, logr2;

    if (r1_zero) {
        logr1 = 0.0;
        theta1 = kPi;
        theta2 = kPi;
    } else {
        logr1 = std::log(r1);
    }
    if (r2_zero) {
        logr2 = 0.0;
        theta1 = 0.0;
        theta2 = 0.0;
    } else {
        logr2 = std::log(r2);
    }

    double inv_2pi = 0.5 / kPi;

    double P1 =
        inv_2pi * (x * (theta1 - theta2) + theta2 * d + z * logr1 - z * logr2);
    double P2 = x * P1 + inv_2pi * (0.5 * r2 * r2 * theta2 -
                                    0.5 * r1 * r1 * theta1 - 0.5 * z * d);

    result.a = P1 - P2 / d;
    result.b = P2 / d;

    /* Differentiate following same pattern as LinVortexStream but with
     * theta instead of logr in key places */

    /* For the derivative, we treat the branch-cut adjustment and the r~0
     * special cases as fixed (non-differentiable). When r1 or r2 are not near
     * zero, the theta derivatives come from the PanelInfoDeriv normally.
     * The branch cut adjustment theta += 2*pi is a constant shift so
     * dtheta/d(...) is unchanged. */

    double inv_d = 1.0 / d;
    double inv_d2 = inv_d * inv_d;

    /* Xj derivatives */
    Vector4d logr1_Xj =
        (!r1_zero && r1 > ep) ? (pid.r1_Xj / r1).eval() : Vector4d::Zero();
    Vector4d logr2_Xj =
        (!r2_zero && r2 > ep) ? (pid.r2_Xj / r2).eval() : Vector4d::Zero();
    Vector4d r1sq_Xj =
        r1_zero ? Vector4d::Zero() : (2.0 * r1 * pid.r1_Xj).eval();
    Vector4d r2sq_Xj =
        r2_zero ? Vector4d::Zero() : (2.0 * r2 * pid.r2_Xj).eval();
    Vector4d th1_Xj = r1_zero ? Vector4d::Zero() : pid.theta1_Xj;
    Vector4d th2_Xj = r2_zero ? Vector4d::Zero() : pid.theta2_Xj;

    {
        Vector4d dP1 =
            inv_2pi * (pid.x1_Xj * (theta1 - theta2) + x * (th1_Xj - th2_Xj) +
                       th2_Xj * d + theta2 * pid.len_Xj + pid.y_Xj * logr1 +
                       z * logr1_Xj - pid.y_Xj * logr2 - z * logr2_Xj);

        Vector4d dP2 =
            pid.x1_Xj * P1 + x * dP1 +
            inv_2pi * (0.5 * r2sq_Xj * theta2 + 0.5 * r2 * r2 * th2_Xj -
                       0.5 * r1sq_Xj * theta1 - 0.5 * r1 * r1 * th1_Xj -
                       0.5 * pid.y_Xj * d - 0.5 * z * pid.len_Xj);

        result.a_Xj = dP1 - (dP2 * d - P2 * pid.len_Xj) * inv_d2;
        result.b_Xj = (dP2 * d - P2 * pid.len_Xj) * inv_d2;
    }

    /* xi derivatives */
    Vector2d logr1_xi =
        (!r1_zero && r1 > ep) ? (pid.r1_xi / r1).eval() : Vector2d::Zero();
    Vector2d logr2_xi =
        (!r2_zero && r2 > ep) ? (pid.r2_xi / r2).eval() : Vector2d::Zero();
    Vector2d r1sq_xi =
        r1_zero ? Vector2d::Zero() : (2.0 * r1 * pid.r1_xi).eval();
    Vector2d r2sq_xi =
        r2_zero ? Vector2d::Zero() : (2.0 * r2 * pid.r2_xi).eval();
    Vector2d th1_xi = r1_zero ? Vector2d::Zero() : pid.theta1_xi;
    Vector2d th2_xi = r2_zero ? Vector2d::Zero() : pid.theta2_xi;

    {
        Vector2d dP1 =
            inv_2pi * (pid.x1_xi * (theta1 - theta2) + x * (th1_xi - th2_xi) +
                       th2_xi * d + pid.y_xi * logr1 + z * logr1_xi -
                       pid.y_xi * logr2 - z * logr2_xi);

        Vector2d dP2 =
            pid.x1_xi * P1 + x * dP1 +
            inv_2pi * (0.5 * r2sq_xi * theta2 + 0.5 * r2 * r2 * th2_xi -
                       0.5 * r1sq_xi * theta1 - 0.5 * r1 * r1 * th1_xi -
                       0.5 * pid.y_xi * d);

        result.a_xi = dP1 - dP2 * inv_d;
        result.b_xi = dP2 * inv_d;
    }

    return result;
}

/* ============================================================================
 * PanelLinSourceVelocityDeriv
 * ============================================================================
 */
LinSourceVelocityDeriv PanelLinSourceVelocityDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi) {
    LinSourceVelocityDeriv result;

    PanelInfoDeriv pid = GetPanelInfoDeriv(Xj, xi);
    const PanelInfo& info = pid.info;

    if (info.len < kEps) {
        result.a.setZero();
        result.b.setZero();
        result.a_Xj.setZero();
        result.b_Xj.setZero();
        result.a_xi.setZero();
        result.b_xi.setZero();
        return result;
    }

    double x = info.x1;
    double z = info.y;
    double d = info.len;
    double inv_len = 1.0 / d;
    double inv_2pi = 0.5 / kPi;
    double dtheta = info.theta2 - info.theta1;

    /* log(r1/r2) = log(r1) - log(r2) */
    double logr12 = SafeLog(info.r1) - SafeLog(info.r2);

    Vector2d dx_vec = Xj.col(1) - Xj.col(0);
    Vector2d t = dx_vec * inv_len;
    Vector2d n(-t(1), t(0));

    /* Panel-coord velocity coefficients */
    double temp1_u = logr12 * inv_2pi;
    double temp2_u = (x * logr12 - d + z * dtheta) * inv_2pi * inv_len;
    double ug1 = temp1_u - temp2_u;
    double ug2 = temp2_u;

    double temp1_v = dtheta * inv_2pi;
    double temp2_v = (-z * logr12 + x * dtheta) * inv_2pi * inv_len;
    double wg1 = temp1_v - temp2_v;
    double wg2 = temp2_v;

    result.a = ug1 * t + wg1 * n;
    result.b = ug2 * t + wg2 * n;

    /* dt/dXj, dn/dXj */
    Eigen::Matrix<double, 2, 4> dt_dXj;
    dt_dXj.col(0) = (-Vector2d::Unit(0) + t * t(0)) * inv_len;
    dt_dXj.col(1) = (-Vector2d::Unit(1) + t * t(1)) * inv_len;
    dt_dXj.col(2) = (Vector2d::Unit(0) - t * t(0)) * inv_len;
    dt_dXj.col(3) = (Vector2d::Unit(1) - t * t(1)) * inv_len;
    Eigen::Matrix<double, 2, 4> dn_dXj;
    for (int k = 0; k < 4; ++k) {
        dn_dXj(0, k) = -dt_dXj(1, k);
        dn_dXj(1, k) = dt_dXj(0, k);
    }

    /* logr12 = log(r1) - log(r2) derivatives */
    Vector4d logr12_Xj = Vector4d::Zero();
    Vector2d logr12_xi = Vector2d::Zero();
    if (info.r1 > kEps) {
        logr12_Xj += pid.r1_Xj / info.r1;
        logr12_xi += pid.r1_xi / info.r1;
    }
    if (info.r2 > kEps) {
        logr12_Xj -= pid.r2_Xj / info.r2;
        logr12_xi -= pid.r2_xi / info.r2;
    }

    Vector4d dtheta_Xj = pid.theta2_Xj - pid.theta1_Xj;
    Vector2d dtheta_xi = pid.theta2_xi - pid.theta1_xi;

    /* Differentiate velocity coefficients */
    auto diffCoeffs = [&](const auto& x_d, const auto& z_d, const auto& d_d,
                          const auto& dtheta_d, const auto& logr12_d,
                          auto& ug1_d, auto& ug2_d, auto& wg1_d, auto& wg2_d) {
        auto dtemp1_u = logr12_d * inv_2pi;
        double numer_u = x * logr12 - d + z * dtheta;
        auto dnumer_u =
            x_d * logr12 + x * logr12_d - d_d + z_d * dtheta + z * dtheta_d;
        auto dtemp2_u =
            inv_2pi * (dnumer_u * inv_len - numer_u * d_d * inv_len * inv_len);
        ug1_d = dtemp1_u - dtemp2_u;
        ug2_d = dtemp2_u;

        auto dtemp1_v = dtheta_d * inv_2pi;
        double numer_v = -z * logr12 + x * dtheta;
        auto dnumer_v =
            -z_d * logr12 - z * logr12_d + x_d * dtheta + x * dtheta_d;
        auto dtemp2_v =
            inv_2pi * (dnumer_v * inv_len - numer_v * d_d * inv_len * inv_len);
        wg1_d = dtemp1_v - dtemp2_v;
        wg2_d = dtemp2_v;
    };

    Vector4d ug1_Xj, ug2_Xj, wg1_Xj, wg2_Xj;
    diffCoeffs(pid.x1_Xj, pid.y_Xj, pid.len_Xj, dtheta_Xj, logr12_Xj, ug1_Xj,
               ug2_Xj, wg1_Xj, wg2_Xj);

    Vector2d len_xi = Vector2d::Zero();
    Vector2d ug1_xi, ug2_xi, wg1_xi, wg2_xi;
    diffCoeffs(pid.x1_xi, pid.y_xi, len_xi, dtheta_xi, logr12_xi, ug1_xi,
               ug2_xi, wg1_xi, wg2_xi);

    for (int k = 0; k < 4; ++k) {
        result.a_Xj.col(k) = ug1_Xj(k) * t + ug1 * dt_dXj.col(k) +
                             wg1_Xj(k) * n + wg1 * dn_dXj.col(k);
        result.b_Xj.col(k) = ug2_Xj(k) * t + ug2 * dt_dXj.col(k) +
                             wg2_Xj(k) * n + wg2 * dn_dXj.col(k);
    }
    for (int k = 0; k < 2; ++k) {
        result.a_xi.col(k) = ug1_xi(k) * t + wg1_xi(k) * n;
        result.b_xi.col(k) = ug2_xi(k) * t + wg2_xi(k) * n;
    }

    return result;
}

} /* namespace scirocco */
