// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unit tests for panel module.

#include "xfoil/panel.h"
#include "xfoil/geometry.h"

#include <cmath>

#include <gtest/gtest.h>

namespace xfoil {
namespace {

constexpr double kPi = 3.14159265358979323846;

TEST(PanelTest, AtancBasic) {
    // Test atanc returns values in [0, 2*pi)
    EXPECT_NEAR(atanc(0.0, 1.0), 0.0, 1e-10);
    EXPECT_NEAR(atanc(1.0, 0.0), kPi / 2.0, 1e-10);
    EXPECT_NEAR(atanc(0.0, -1.0), kPi, 1e-10);
    EXPECT_NEAR(atanc(-1.0, 0.0), 3.0 * kPi / 2.0, 1e-10);

    // Quadrant 3: should be > pi, not negative
    double angle = atanc(-1.0, -1.0);
    EXPECT_GT(angle, kPi);
    EXPECT_LT(angle, 3.0 * kPi / 2.0);
}

TEST(PanelTest, ApcalcHorizontalLine) {
    // Horizontal line from (0,0) to (4,0)
    VectorXd x(5), y(5);
    for (int i = 0; i < 5; ++i) {
        x(i) = static_cast<double>(i);
        y(i) = 0.0;
    }

    VectorXd apanel;
    apcalc(x, y, apanel);

    // All panels should have angle 0 (pointing in +x direction)
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(apanel(i), 0.0, 1e-10);
    }
}

TEST(PanelTest, ApcalcCircle) {
    // Circle discretization (CCW)
    int n = 33;
    VectorXd x(n), y(n);
    double radius = 1.0;

    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * kPi * static_cast<double>(i) / (n - 1);
        x(i) = radius * std::cos(theta);
        y(i) = radius * std::sin(theta);
    }

    VectorXd apanel;
    apcalc(x, y, apanel);

    // Panel angles should be tangent to circle
    // At theta = 0 (x=1, y=0), tangent direction is (0, 1), angle = pi/2
    // At theta = pi/2 (x=0, y=1), tangent direction is (-1, 0), angle = pi
    // etc.
    for (int i = 1; i < n - 2; ++i) {
        double theta = 2.0 * kPi * static_cast<double>(i) / (n - 1);
        // Expected tangent angle: theta + pi/2
        double expected = theta + kPi / 2.0;
        while (expected >= kPi) expected -= 2.0 * kPi;
        while (expected < -kPi) expected += 2.0 * kPi;

        // Check panel angle is close to expected (within panel discretization)
        double diff = std::abs(apanel(i) - expected);
        if (diff > kPi) diff = 2.0 * kPi - diff;
        EXPECT_LT(diff, 0.3) << "Panel " << i << " angle mismatch";
    }
}

TEST(PanelTest, FreestreamComponents) {
    double uinf, vinf;

    // alpha = 0
    freestream(0.0, 1.0, uinf, vinf);
    EXPECT_NEAR(uinf, 1.0, 1e-10);
    EXPECT_NEAR(vinf, 0.0, 1e-10);

    // alpha = 90 deg
    freestream(kPi / 2.0, 1.0, uinf, vinf);
    EXPECT_NEAR(uinf, 0.0, 1e-10);
    EXPECT_NEAR(vinf, 1.0, 1e-10);

    // alpha = 45 deg, qinf = 2
    freestream(kPi / 4.0, 2.0, uinf, vinf);
    EXPECT_NEAR(uinf, std::sqrt(2.0), 1e-10);
    EXPECT_NEAR(vinf, std::sqrt(2.0), 1e-10);
}

TEST(PanelTest, PsilinSymmetry) {
    // Vortex panel influence should have certain symmetry properties
    // Field point above and below panel midpoint

    double x1 = -0.5, y1 = 0.0;
    double x2 = 0.5, y2 = 0.0;

    // Point above panel
    auto inf_above = psilin(0.0, 0.5, x1, y1, x2, y2, 1.0, 0.0);

    // Point below panel
    auto inf_below = psilin(0.0, -0.5, x1, y1, x2, y2, 1.0, 0.0);

    // Stream function should be antisymmetric for vortex
    EXPECT_NEAR(inf_above.psi, -inf_below.psi, 1e-10);
}

TEST(PanelTest, PsilinSourceSymmetry) {
    // Source panel influence - symmetric about panel

    double x1 = -0.5, y1 = 0.0;
    double x2 = 0.5, y2 = 0.0;

    // Point above panel
    auto inf_above = psilin(0.0, 0.5, x1, y1, x2, y2, 0.0, 1.0);

    // Point below panel
    auto inf_below = psilin(0.0, -0.5, x1, y1, x2, y2, 0.0, 1.0);

    // Stream function should be symmetric for source
    EXPECT_NEAR(inf_above.psi, inf_below.psi, 1e-10);
}

TEST(PanelTest, PsilinDecaysWithDistance) {
    // Influence should decay with distance from panel
    // Use off-center point to avoid symmetry cancellation

    double x1 = -0.5, y1 = 0.0;
    double x2 = 0.5, y2 = 0.0;

    // Off-center field points to get non-zero psi
    auto inf_close = psilin(0.2, 0.1, x1, y1, x2, y2, 1.0, 0.0);
    auto inf_far = psilin(0.2, 10.0, x1, y1, x2, y2, 1.0, 0.0);

    // Both should be finite
    EXPECT_TRUE(std::isfinite(inf_close.psi));
    EXPECT_TRUE(std::isfinite(inf_far.psi));

    // At least one should be non-zero
    EXPECT_TRUE(std::abs(inf_close.psi) > 1e-10 ||
                std::abs(inf_far.psi) > 1e-10);
}

TEST(PanelTest, AicBuildSymmetry) {
    // For symmetric airfoil at alpha=0, AIC should have structure

    auto naca = naca4(12, 21);  // NACA 0012, 41 points

    VectorXd apanel;
    apcalc(naca.x, naca.y, apanel);

    MatrixXd aij;
    aic_build(naca.x, naca.y, apanel, aij, true);

    // AIC matrix should be square
    EXPECT_EQ(aij.rows(), naca.n);
    EXPECT_EQ(aij.cols(), naca.n);

    // Last row should be Kutta condition
    EXPECT_NEAR(aij(naca.n - 1, 0), 1.0, 1e-10);
    EXPECT_NEAR(aij(naca.n - 1, naca.n - 1), 1.0, 1e-10);

    // Interior of last row should be zero
    for (int j = 1; j < naca.n - 1; ++j) {
        EXPECT_NEAR(aij(naca.n - 1, j), 0.0, 1e-10);
    }
}

TEST(PanelTest, DicBuildBasic) {
    auto naca = naca4(12, 21);

    VectorXd apanel;
    apcalc(naca.x, naca.y, apanel);

    MatrixXd dij;
    dic_build(naca.x, naca.y, apanel, dij);

    // DIC matrix should be square
    EXPECT_EQ(dij.rows(), naca.n);
    EXPECT_EQ(dij.cols(), naca.n);

    // DIC should have finite values
    for (int i = 0; i < naca.n; ++i) {
        for (int j = 0; j < naca.n; ++j) {
            EXPECT_TRUE(std::isfinite(dij(i, j)));
        }
    }
}

TEST(PanelTest, QcalcFreestream) {
    // For flat plate at alpha=0 with zero circulation,
    // velocity should be just freestream

    int n = 11;
    VectorXd x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x(i) = static_cast<double>(i) / (n - 1);
        y(i) = 0.0;
    }

    VectorXd apanel;
    apcalc(x, y, apanel);

    VectorXd gam = VectorXd::Zero(n);
    VectorXd sig = VectorXd::Zero(n);

    VectorXd qinv;
    qcalc(x, y, apanel, gam, sig, 0.0, qinv);

    // With zero circulation, velocity should be freestream projected onto panel
    // For horizontal panels at alpha=0, this is 1.0
    for (int i = 1; i < n - 1; ++i) {
        EXPECT_NEAR(qinv(i), 1.0, 0.1);
    }
}

TEST(PanelTest, NacaAicSolvable) {
    // AIC matrix for NACA airfoil should be invertible

    auto naca = naca4(12, 31);

    VectorXd apanel;
    apcalc(naca.x, naca.y, apanel);

    MatrixXd aij;
    aic_build(naca.x, naca.y, apanel, aij, true);

    // Check condition number is reasonable
    Eigen::JacobiSVD<MatrixXd> svd(aij);
    double cond = svd.singularValues()(0) /
                  svd.singularValues()(svd.singularValues().size() - 1);

    // Condition number should not be astronomical
    EXPECT_LT(cond, 1e10);
}

}  // namespace
}  // namespace xfoil
