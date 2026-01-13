// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/spline.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cppfoil {

namespace {

// Find segment index for parameter value
int FindSegment(const VectorXd& s, double sval) {
    int n = static_cast<int>(s.size());
    if (sval <= s(0)) return 0;
    if (sval >= s(n - 1)) return n - 2;

    // Binary search
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

// Build natural cubic spline from data points
// Returns coefficient matrix [4 x (N-1)]
MatrixXd BuildCubicSpline(const VectorXd& s, const VectorXd& y) {
    int n = static_cast<int>(s.size());
    if (n < 2) {
        throw std::invalid_argument("Need at least 2 points for spline");
    }

    // Build tridiagonal system for second derivatives
    VectorXd h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h(i) = s(i + 1) - s(i);
    }

    // Natural spline: second derivative = 0 at endpoints
    VectorXd d2y = VectorXd::Zero(n);

    if (n > 2) {
        // Set up tridiagonal system
        int m = n - 2;
        VectorXd diag(m), upper(m - 1), lower(m - 1), rhs(m);

        for (int i = 0; i < m; ++i) {
            diag(i) = 2.0 * (h(i) + h(i + 1));
            rhs(i) = 6.0 * ((y(i + 2) - y(i + 1)) / h(i + 1) -
                            (y(i + 1) - y(i)) / h(i));
        }
        for (int i = 0; i < m - 1; ++i) {
            upper(i) = h(i + 1);
            lower(i) = h(i + 1);
        }

        // Thomas algorithm for tridiagonal system
        VectorXd cp(m), dp(m);
        if (m == 1) {
            // Special case: single interior point
            d2y(1) = rhs(0) / diag(0);
        } else {
            cp(0) = upper(0) / diag(0);
            dp(0) = rhs(0) / diag(0);
            for (int i = 1; i < m; ++i) {
                double denom = diag(i) - lower(i - 1) * cp(i - 1);
                if (i < m - 1) {
                    cp(i) = upper(i) / denom;
                }
                dp(i) = (rhs(i) - lower(i - 1) * dp(i - 1)) / denom;
            }

            // Back substitution
            d2y(m) = dp(m - 1);
            for (int i = m - 2; i >= 0; --i) {
                d2y(i + 1) = dp(i) - cp(i) * d2y(i + 2);
            }
        }
    }

    // Compute coefficients for each segment
    // f(t) = c0*t^3 + c1*t^2 + c2*t + c3
    MatrixXd c(4, n - 1);
    for (int i = 0; i < n - 1; ++i) {
        double hi = h(i);
        c(3, i) = y(i);                             // constant term
        c(2, i) = (y(i + 1) - y(i)) / hi -
                  hi * (2.0 * d2y(i) + d2y(i + 1)) / 6.0;  // linear term
        c(1, i) = d2y(i) / 2.0;                     // quadratic term
        c(0, i) = (d2y(i + 1) - d2y(i)) / (6.0 * hi);     // cubic term
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
    // 5-point Gauss-Legendre quadrature on [0,1]
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

    // Initial arc-length estimate
    VectorXd s = VectorXd::Zero(n);
    for (int i = 1; i < n; ++i) {
        s(i) = s(i - 1) + (points.col(i) - points.col(i - 1)).norm();
    }

    // Build initial splines
    spline.X.s = s;
    spline.X.c = BuildCubicSpline(s, points.row(0));
    spline.Y.s = s;
    spline.Y.c = BuildCubicSpline(s, points.row(1));

    // Re-integrate to true arc-length via several passes
    VectorXd xq, wq;
    QuadSeg(xq, wq);

    VectorXd s_new = VectorXd::Zero(n);
    for (int pass = 0; pass < 10; ++pass) {
        double serr = 0.0;
        s_new(0) = s(0);

        for (int i = 0; i < n - 1; ++i) {
            double ds = s(i + 1) - s(i);

            // Integrate arc length over segment using quadrature
            double sint = 0.0;
            for (int q = 0; q < 5; ++q) {
                double sq = s(i) + xq(q) * ds;
                double xs = spline.X.EvalDeriv(sq);
                double ys = spline.Y.EvalDeriv(sq);
                sint += wq(q) * std::sqrt(xs * xs + ys * ys);
            }
            sint *= ds;

            serr = std::max(serr, std::abs(sint - ds));
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

}  // namespace cppfoil
