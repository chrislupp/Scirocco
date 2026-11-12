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
#include "scirocco/core/spline.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace scirocco {

namespace {

/* Find segment index for parameter value */
int FindSegment(const VectorXd& s, double sval) {
    int n = static_cast<int>(s.size());
    if (sval <= s(0)) return 0;
    if (sval >= s(n - 1)) return n - 2;

    /* Binary search */
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (s(mid) > sval) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return lo;
}

/* Build cubic spline with not-a-knot boundary conditions (like
 * scipy.CubicSpline) Returns coefficient matrix [4 x (N-1)] */
MatrixXd BuildCubicSpline(const VectorXd& s, const VectorXd& y) {
    int n = static_cast<int>(s.size());
    if (n < 2) {
        throw std::invalid_argument("Need at least 2 points for spline");
    }

    /* Compute intervals */
    VectorXd h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h(i) = s(i + 1) - s(i);
    }

    VectorXd d2y = VectorXd::Zero(n);

    if (n == 2) {
        /* Linear interpolation (no curvature) */
        d2y.setZero();
    } else if (n == 3) {
        /* Quadratic - single d2y value (constant curvature) */
        double dy0 = (y(1) - y(0)) / h(0);
        double dy1 = (y(2) - y(1)) / h(1);
        double d2 = 2.0 * (dy1 - dy0) / (h(0) + h(1));
        d2y(0) = d2;
        d2y(1) = d2;
        d2y(2) = d2;
    } else {
        /* Not-a-knot: the endpoint equations only couple three adjacent
         * second derivatives. Eliminate d2y[0] and d2y[n-1] from the first
         * and last interior equations, leaving a tridiagonal system for
         * d2y[1..n-2]. This is algebraically equivalent to the former dense
         * solve, but costs O(n) rather than O(n^3). */
        const int m = n - 2;
        VectorXd lower = VectorXd::Zero(m);
        VectorXd diagonal = VectorXd::Zero(m);
        VectorXd upper = VectorXd::Zero(m);
        VectorXd rhs = VectorXd::Zero(m);

        for (int i = 1; i < n - 1; ++i) {
            rhs(i - 1) =
                6.0 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i - 1));
        }

        /* First row after eliminating d2y[0]. */
        diagonal(0) = (h(0) + h(1)) * (h(0) + 2.0 * h(1)) / h(1);
        upper(0) = h(1) - h(0) * h(0) / h(1);

        /* Unmodified interior rows. */
        for (int i = 2; i < n - 2; ++i) {
            const int row = i - 1;
            lower(row) = h(i - 1);
            diagonal(row) = 2.0 * (h(i - 1) + h(i));
            upper(row) = h(i);
        }

        /* Last row after eliminating d2y[n-1]. */
        lower(m - 1) = h(n - 3) - h(n - 2) * h(n - 2) / h(n - 3);
        diagonal(m - 1) =
            (h(n - 3) + h(n - 2)) * (2.0 * h(n - 3) + h(n - 2)) / h(n - 3);

        /* Thomas algorithm. The system is nonsingular for strictly increasing
         * knots, as required by the spline interpolation problem. */
        for (int i = 1; i < m; ++i) {
            const double factor = lower(i) / diagonal(i - 1);
            diagonal(i) -= factor * upper(i - 1);
            rhs(i) -= factor * rhs(i - 1);
        }
        VectorXd interior(m);
        interior(m - 1) = rhs(m - 1) / diagonal(m - 1);
        for (int i = m - 2; i >= 0; --i) {
            interior(i) = (rhs(i) - upper(i) * interior(i + 1)) / diagonal(i);
        }

        d2y.segment(1, m) = interior;
        d2y(0) = ((h(0) + h(1)) * d2y(1) - h(0) * d2y(2)) / h(1);
        d2y(n - 1) =
            ((h(n - 3) + h(n - 2)) * d2y(n - 2) - h(n - 2) * d2y(n - 3)) /
            h(n - 3);
    }

    /* Compute coefficients for each segment
     * f(t) = c0*t^3 + c1*t^2 + c2*t + c3 where t = s - s[i] */
    MatrixXd c(4, n - 1);
    for (int i = 0; i < n - 1; ++i) {
        double hi = h(i);
        c(3, i) = y(i); /* constant term */
        c(2, i) = (y(i + 1) - y(i)) / hi -
                  hi * (2.0 * d2y(i) + d2y(i + 1)) / 6.0; /* linear term */
        c(1, i) = d2y(i) / 2.0;                           /* quadratic term */
        c(0, i) = (d2y(i + 1) - d2y(i)) / (6.0 * hi);     /* cubic term */
    }

    return c;
}

}  // namespace

double CubicSpline1D::Eval(double sval) const {
    int i = FindSegment(s, sval);
    double t = sval - s(i);
    return ((c(0, i) * t + c(1, i)) * t + c(2, i)) * t + c(3, i);
}

double CubicSpline1D::EvalDeriv(double sval) const {
    int i = FindSegment(s, sval);
    double t = sval - s(i);
    return (3.0 * c(0, i) * t + 2.0 * c(1, i)) * t + c(2, i);
}

double CubicSpline1D::EvalDeriv2(double sval) const {
    int i = FindSegment(s, sval);
    double t = sval - s(i);
    return 6.0 * c(0, i) * t + 2.0 * c(1, i);
}

void QuadSeg(VectorXd& xq, VectorXd& wq) {
    /* 5-point Gauss-Legendre quadrature on [0,1] */
    xq.resize(5);
    wq.resize(5);
    xq << 0.046910077030668, 0.230765344947158, 0.500000000000000,
        0.769234655052842, 0.953089922969332;
    wq << 0.118463442528095, 0.239314335249683, 0.284444444444444,
        0.239314335249683, 0.118463442528095;
}

Spline2D Spline2DCreate(const MatrixXd& points) {
    int n = static_cast<int>(points.cols());
    if (n < 2) {
        throw std::invalid_argument("Need at least 2 points for spline");
    }

    Spline2D spline;

    /* Initial arc-length estimate */
    VectorXd s = VectorXd::Zero(n);
    for (int i = 1; i < n; ++i) {
        s(i) = s(i - 1) + (points.col(i) - points.col(i - 1)).norm();
    }

    /* Build initial splines */
    spline.X.s = s;
    spline.X.c = BuildCubicSpline(s, points.row(0));
    spline.Y.s = s;
    spline.Y.c = BuildCubicSpline(s, points.row(1));

    /* Re-integrate to true arc-length via several passes */
    VectorXd xq, wq;
    QuadSeg(xq, wq);

    VectorXd s_new = VectorXd::Zero(n);
    /* Always ten passes, as in mfoil's spline2d. SplineArcLengthDeriv
     * differentiates exactly ten passes, so stopping early would leave the
     * primal and its derivative describing different knots (by ~1e-12). */
    constexpr int kArcLengthPasses = 10;
    for (int pass = 0; pass < kArcLengthPasses; ++pass) {
        s_new(0) = s(0);

        for (int i = 0; i < n - 1; ++i) {
            double ds = s(i + 1) - s(i);

            /* Integrate arc length over segment using quadrature */
            double sint = 0.0;
            for (int q = 0; q < 5; ++q) {
                double sq = s(i) + xq(q) * ds;
                double xs = spline.X.EvalDeriv(sq);
                double ys = spline.Y.EvalDeriv(sq);
                sint += wq(q) * std::sqrt(xs * xs + ys * ys);
            }
            sint *= ds;

            s_new(i + 1) = s_new(i) + sint;
        }

        s = s_new;
        spline.X.s = s;
        spline.X.c = BuildCubicSpline(s, points.row(0));
        spline.Y.s = s;
        spline.Y.c = BuildCubicSpline(s, points.row(1));
    }

    return spline;
}

MatrixXd SplineEval(const Spline2D& spline, const VectorXd& s_vals) {
    int n = static_cast<int>(s_vals.size());
    MatrixXd result(2, n);
    for (int i = 0; i < n; ++i) {
        result(0, i) = spline.X.Eval(s_vals(i));
        result(1, i) = spline.Y.Eval(s_vals(i));
    }
    return result;
}

MatrixXd SplineTangent(const Spline2D& spline, const VectorXd& s_vals) {
    int n = static_cast<int>(s_vals.size());
    MatrixXd result(2, n);
    for (int i = 0; i < n; ++i) {
        result(0, i) = spline.X.EvalDeriv(s_vals(i));
        result(1, i) = spline.Y.EvalDeriv(s_vals(i));
    }
    return result;
}

CubicSpline1D Spline1DCreate(const VectorXd& x, const VectorXd& y) {
    CubicSpline1D spline;
    spline.s = x;
    spline.c = BuildCubicSpline(x, y);
    return spline;
}

}  // namespace scirocco
