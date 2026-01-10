// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unit tests for spline module.

#include "xfoil/spline.h"

#include <cmath>

#include <gtest/gtest.h>

namespace xfoil {
namespace {

TEST(SplineTest, LinearInterpolation) {
    int n = 5;
    VectorXd x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
        y(i) = 2.0 * x(i) + 1.0;
    }
    s = x;

    spline(y, s, ys);

    for (int i = 0; i < n - 1; ++i) {
        double ss = x(i) + 0.5;
        double yi = seval(ss, y, ys, s);
        double y_expected = 2.0 * ss + 1.0;
        EXPECT_NEAR(yi, y_expected, 1e-10);
    }

    for (int i = 0; i < n - 1; ++i) {
        double ss = x(i) + 0.5;
        double dyi = deval(ss, y, ys, s);
        EXPECT_NEAR(dyi, 2.0, 1e-10);
    }
}

TEST(SplineTest, QuadraticInterpolation) {
    int n = 11;
    VectorXd x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i) / 10.0;
        y(i) = x(i) * x(i);
    }

    s = x;
    spline(y, s, ys);

    double ss = 0.35;
    double yi = seval(ss, y, ys, s);
    double y_expected = ss * ss;
    EXPECT_NEAR(yi, y_expected, 1e-4);

    double dyi = deval(ss, y, ys, s);
    double dy_expected = 2.0 * ss;
    EXPECT_NEAR(dyi, dy_expected, 1e-3);
}

TEST(SplineTest, ArcLengthCalculation) {
    int n = 5;
    VectorXd x(n), y(n), s(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
        y(i) = 0.0;
    }

    scalc(x, y, s);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(s(i), x(i), 1e-10);
    }

    for (int i = 0; i < n; ++i) {
        y(i) = static_cast<double>(i);
    }

    scalc(x, y, s);

    for (int i = 0; i < n; ++i) {
        double expected = std::sqrt(2.0) * static_cast<double>(i);
        EXPECT_NEAR(s(i), expected, 1e-10);
    }
}

TEST(SplineTest, SimpleAveragingSpline) {
    int n = 5;
    VectorXd x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
    }
    y << 0.0, 0.0, 1.0, 1.0, 1.0;
    s = x;

    splina(y, s, ys);

    EXPECT_GE(ys(1), 0.0);
    EXPECT_GE(ys(2), 0.0);
}

TEST(SplineTest, CircleCurvature) {
    int n = 37;
    VectorXd x(n), y(n), s(n), xs(n), ys(n);

    double radius = 1.0;
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / (n - 1);
        x(i) = radius * std::cos(theta);
        y(i) = radius * std::sin(theta);
    }

    scalc(x, y, s);
    spline(x, s, xs);
    spline(y, s, ys);

    double expected_curv = 1.0 / radius;
    double ss = s(n / 2);
    double kappa = curv(ss, x, xs, y, ys, s);

    EXPECT_NEAR(std::abs(kappa), expected_curv, 0.1);
}

TEST(SplineTest, InverseSpline) {
    int n = 11;
    VectorXd x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i) / 10.0;
        y(i) = x(i) * x(i);
    }

    s = x;
    spline(y, s, ys);

    double si = 0.4;
    double xi = 0.25;
    bool converged = sinvrt(si, xi, y, ys, s);

    EXPECT_TRUE(converged);
    EXPECT_NEAR(si, 0.5, 1e-4);
}

TEST(SplineTest, TridiagonalSolver) {
    VectorXd a(3), b(3), c(3), d(3);

    a << 2.0, 2.0, 2.0;
    b << 0.0, 1.0, 1.0;
    c << 1.0, 1.0, 0.0;
    d << 1.0, 2.0, 1.0;

    trisol(a, b, c, d);

    EXPECT_NEAR(d(0) + d(1), 1.0, 1e-10);
    EXPECT_NEAR(d(0) + 2 * d(1) + d(2), 2.0, 1e-10);
    EXPECT_NEAR(d(1) + 2 * d(2), 1.0, 1e-10);
}

}  // namespace
}  // namespace xfoil
