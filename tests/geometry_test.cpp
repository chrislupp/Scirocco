// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unit tests for geometry module.

#include "xfoil/geometry.h"
#include "xfoil/spline.h"

#include <cmath>

#include <gtest/gtest.h>

namespace xfoil {
namespace {

TEST(GeometryTest, Naca4Generation) {
    auto result = naca4(12, 51);

    EXPECT_EQ(result.n, 101);
    EXPECT_EQ(result.name, "NACA 0012");

    EXPECT_NEAR(result.x(0), 1.0, 1e-6);
    EXPECT_NEAR(result.x(result.n - 1), 1.0, 1e-6);

    int ile = result.n / 2;
    EXPECT_NEAR(result.x(ile), 0.0, 1e-6);

    EXPECT_NEAR(result.y(0), 0.0, 0.01);
    EXPECT_NEAR(result.y(result.n - 1), 0.0, 0.01);

    double max_thick = 0.0;
    for (int i = 0; i < result.n / 2; ++i) {
        if (result.y(i) > max_thick) {
            max_thick = result.y(i);
        }
    }
    EXPECT_NEAR(max_thick, 0.06, 0.005);
}

TEST(GeometryTest, Naca4Cambered) {
    auto result = naca4(2412, 51);

    EXPECT_EQ(result.name, "NACA 2412");
    EXPECT_GT(result.y(0), result.y(result.n - 1));

    double y_upper = 0.0, y_lower = 0.0;
    for (int i = 0; i < result.n; ++i) {
        if (std::abs(result.x(i) - 0.5) < 0.02) {
            if (i < result.n / 2) {
                y_upper = std::max(y_upper, result.y(i));
            } else {
                y_lower = std::min(y_lower, result.y(i));
            }
        }
    }
    EXPECT_GT(y_upper, -y_lower);
}

TEST(GeometryTest, Naca5Generation) {
    auto result = naca5(23012, 51);

    EXPECT_EQ(result.n, 101);
    EXPECT_EQ(result.name, "NACA 23012");

    EXPECT_NEAR(result.x(0), 1.0, 1e-6);
    EXPECT_NEAR(result.x(result.n - 1), 1.0, 1e-6);
}

TEST(GeometryTest, Naca5Invalid) {
    auto result = naca5(12345, 51);

    EXPECT_EQ(result.n, 0);
}

TEST(GeometryTest, LeadingEdgeFind) {
    auto naca = naca4(12, 51);

    VectorXd s, xp, yp;
    scalc(naca.x, naca.y, s);
    spline(naca.x, s, xp);
    spline(naca.y, s, yp);

    double sle = lefind(naca.x, xp, naca.y, yp, s);

    double xle = seval(sle, naca.x, xp, s);
    double yle = seval(sle, naca.y, yp, s);

    EXPECT_NEAR(xle, 0.0, 1e-4);
    EXPECT_NEAR(yle, 0.0, 1e-4);
}

TEST(GeometryTest, GeometricParameters) {
    auto naca = naca4(12, 51);

    VectorXd s, xp, yp;
    scalc(naca.x, naca.y, s);
    spline(naca.x, s, xp);
    spline(naca.y, s, yp);

    auto gp = geopar(naca.x, xp, naca.y, yp, s);

    EXPECT_NEAR(gp.chord, 1.0, 0.01);
    EXPECT_NEAR(gp.xle, 0.0, 1e-4);
    EXPECT_NEAR(gp.xte, 1.0, 1e-4);
    EXPECT_GT(gp.radle, 0.0);
}

TEST(GeometryTest, Normalization) {
    auto naca = naca4(12, 51);

    VectorXd x = naca.x * 2.0;
    VectorXd y = naca.y * 2.0;
    x.array() += 1.0;
    y.array() += 0.5;
    VectorXd s;

    normalize(x, y, s);

    VectorXd xp, yp;
    scalc(x, y, s);
    spline(x, s, xp);
    spline(y, s, yp);
    double sle = lefind(x, xp, y, yp, s);
    double xle = seval(sle, x, xp, s);
    double yle = seval(sle, y, yp, s);

    EXPECT_NEAR(xle, 0.0, 1e-4);
    EXPECT_NEAR(yle, 0.0, 1e-4);

    double xte = 0.5 * (x(0) + x(x.size() - 1));
    EXPECT_NEAR(xte, 1.0, 0.01);
}

TEST(GeometryTest, NormalCalculation) {
    auto naca = naca4(12, 51);

    VectorXd s;
    scalc(naca.x, naca.y, s);

    VectorXd nx, ny;
    ncalc(naca.x, naca.y, s, nx, ny);

    for (int i = 0; i < naca.n; ++i) {
        double mag = std::sqrt(nx(i) * nx(i) + ny(i) * ny(i));
        EXPECT_NEAR(mag, 1.0, 1e-10);
    }

    int imid = naca.n / 4;
    EXPECT_GT(ny(imid), 0.0);

    int imid_lower = 3 * naca.n / 4;
    EXPECT_LT(ny(imid_lower), 0.0);
}

TEST(GeometryTest, TrailingEdgeSharp) {
    auto naca = naca4(12, 51);

    VectorXd s, xp, yp;
    scalc(naca.x, naca.y, s);
    segspl(naca.x, s, xp);
    segspl(naca.y, s, yp);

    auto te = tecalc(naca.x, xp, naca.y, yp, 1.0);

    EXPECT_LT(te.dste, 0.01);

    if (te.sharp) {
        EXPECT_NEAR(te.scs, 1.0, 1e-6);
        EXPECT_NEAR(te.sds, 0.0, 1e-6);
    }
}

TEST(GeometryTest, TrailingEdgeBlunt) {
    auto naca = naca4(12, 51);

    naca.y(naca.n - 1) -= 0.02;

    VectorXd s, xp, yp;
    scalc(naca.x, naca.y, s);
    segspl(naca.x, s, xp);
    segspl(naca.y, s, yp);

    auto te = tecalc(naca.x, xp, naca.y, yp, 1.0);

    EXPECT_FALSE(te.sharp);
    EXPECT_GT(te.dste, 0.01);

    EXPECT_TRUE(std::isfinite(te.scs));
    EXPECT_TRUE(std::isfinite(te.sds));
}

TEST(GeometryTest, CircleNormals) {
    int n = 65;
    VectorXd x(n), y(n);
    double radius = 0.5;

    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / (n - 1);
        x(i) = 0.5 + radius * std::cos(theta);
        y(i) = radius * std::sin(theta);
    }

    VectorXd s, nx, ny;
    scalc(x, y, s);
    ncalc(x, y, s, nx, ny);

    for (int i = 1; i < n - 1; ++i) {
        double rx = x(i) - 0.5;
        double ry = y(i);
        double rmag = std::sqrt(rx * rx + ry * ry);
        rx /= rmag;
        ry /= rmag;

        double dot = nx(i) * rx + ny(i) * ry;
        EXPECT_GT(std::abs(dot), 0.95);
    }
}

TEST(GeometryTest, PanelGeneration) {
    auto naca = naca4(12, 101);

    PanelParams params;
    params.npan = 80;

    auto result = pangen(naca.x, naca.y, params);

    EXPECT_TRUE(result.converged);
    EXPECT_EQ(result.n, 80);

    // Check TE points
    EXPECT_NEAR(result.x(0), 1.0, 0.01);
    EXPECT_NEAR(result.x(result.n - 1), 1.0, 0.01);

    // Check LE is found
    EXPECT_GT(result.sle, 0.0);
    EXPECT_LT(result.sle, result.s(result.n - 1));
}

}  // namespace
}  // namespace xfoil
