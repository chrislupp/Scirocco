// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unit tests for spline module.

#include "xfoil/spline.h"

#include <cmath>

#include <gtest/gtest.h>

namespace xfoil {
namespace {

// Test linear function - spline should interpolate exactly
TEST(SplineTest, LinearInterpolation) {
    // Create linear data y = 2*x + 1
    int n = 5;
    VectorX<double> x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
        y(i) = 2.0 * x(i) + 1.0;
    }

    // Use x as independent variable (s = x)
    s = x;

    // Compute spline coefficients
    spline(y, s, ys);

    // Test interpolation at midpoints
    for (int i = 0; i < n - 1; ++i) {
        double ss = x(i) + 0.5;
        double yi = seval(ss, y, ys, s);
        double y_expected = 2.0 * ss + 1.0;
        EXPECT_NEAR(yi, y_expected, 1e-10);
    }

    // Test derivative - should be constant 2.0
    for (int i = 0; i < n - 1; ++i) {
        double ss = x(i) + 0.5;
        double dyi = deval(ss, y, ys, s);
        EXPECT_NEAR(dyi, 2.0, 1e-10);
    }
}

// Test quadratic function
TEST(SplineTest, QuadraticInterpolation) {
    // Create quadratic data y = x^2
    int n = 11;
    VectorX<double> x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i) / 10.0;
        y(i) = x(i) * x(i);
    }

    s = x;
    spline(y, s, ys);

    // Test interpolation
    double ss = 0.35;
    double yi = seval(ss, y, ys, s);
    double y_expected = ss * ss;
    EXPECT_NEAR(yi, y_expected, 1e-4);

    // Test derivative dy/dx = 2x
    double dyi = deval(ss, y, ys, s);
    double dy_expected = 2.0 * ss;
    EXPECT_NEAR(dyi, dy_expected, 1e-3);
}

// Test arc length calculation
TEST(SplineTest, ArcLengthCalculation) {
    // Create a horizontal line
    int n = 5;
    VectorX<double> x(n), y(n), s(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
        y(i) = 0.0;
    }

    scalc(x, y, s);

    // Arc length should equal x for horizontal line
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(s(i), x(i), 1e-10);
    }

    // Test diagonal line
    for (int i = 0; i < n; ++i) {
        y(i) = static_cast<double>(i);
    }

    scalc(x, y, s);

    // Arc length for 45-degree line should be sqrt(2) * x
    for (int i = 0; i < n; ++i) {
        double expected = std::sqrt(2.0) * static_cast<double>(i);
        EXPECT_NEAR(s(i), expected, 1e-10);
    }
}

// Test simple averaging spline (SPLINA)
TEST(SplineTest, SimpleAveragingSpline) {
    // Create step function to test non-oscillatory behavior
    int n = 5;
    VectorX<double> x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i);
    }
    y << 0.0, 0.0, 1.0, 1.0, 1.0;
    s = x;

    splina(y, s, ys);

    // Derivatives should be non-negative for monotonic data
    // and the averaging should produce smooth transitions
    EXPECT_GE(ys(1), 0.0);
    EXPECT_GE(ys(2), 0.0);
}

// Test curvature calculation for circle
TEST(SplineTest, CircleCurvature) {
    // Create circle points
    int n = 37;  // Every 10 degrees
    VectorX<double> x(n), y(n), s(n), xs(n), ys(n);

    double radius = 1.0;
    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / (n - 1);
        x(i) = radius * std::cos(theta);
        y(i) = radius * std::sin(theta);
    }

    // Calculate arc length and splines
    scalc(x, y, s);
    spline(x, s, xs);
    spline(y, s, ys);

    // Curvature of circle should be 1/radius
    double expected_curv = 1.0 / radius;
    double ss = s(n / 2);  // Test at midpoint
    double kappa = curv(ss, x, xs, y, ys, s);

    // Allow some error due to discretization
    EXPECT_NEAR(std::abs(kappa), expected_curv, 0.1);
}

// Test inverse spline (SINVRT)
TEST(SplineTest, InverseSpline) {
    // Create monotonic data
    int n = 11;
    VectorX<double> x(n), y(n), s(n), ys(n);

    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i) / 10.0;
        y(i) = x(i) * x(i);  // y = x^2
    }

    s = x;
    spline(y, s, ys);

    // Find s where y = 0.25 (should be x = 0.5)
    double si = 0.4;  // Initial guess
    double xi = 0.25;
    bool converged = sinvrt(si, xi, y, ys, s);

    EXPECT_TRUE(converged);
    EXPECT_NEAR(si, 0.5, 1e-4);
}

// Test tridiagonal solver
TEST(SplineTest, TridiagonalSolver) {
    // Simple 3x3 system
    VectorX<double> a(3), b(3), c(3), d(3);

    // Set up: [2 1 0] [x0]   [1]
    //         [1 2 1] [x1] = [2]
    //         [0 1 2] [x2]   [1]
    a << 2.0, 2.0, 2.0;
    b << 0.0, 1.0, 1.0;  // b(0) not used
    c << 1.0, 1.0, 0.0;  // c(n-1) not used
    d << 1.0, 2.0, 1.0;

    trisol(a, b, c, d);

    // Verify solution (d now contains x)
    // For this symmetric system, solution should be approximately [0, 1, 0]
    EXPECT_NEAR(d(0) + d(1), 1.0, 1e-10);
    EXPECT_NEAR(d(0) + 2 * d(1) + d(2), 2.0, 1e-10);
    EXPECT_NEAR(d(1) + 2 * d(2), 1.0, 1e-10);
}

}  // namespace
}  // namespace xfoil
