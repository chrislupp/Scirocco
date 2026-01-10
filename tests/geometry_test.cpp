// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Unit tests for geometry module.

#include "xfoil/geometry.h"

#include <cmath>

#include <gtest/gtest.h>

namespace xfoil {
namespace {

// Test NACA 4-digit airfoil generation
TEST(GeometryTest, Naca4Generation) {
    // Generate NACA 0012 (symmetric airfoil)
    auto result = naca4<double>(12, 51);

    EXPECT_EQ(result.n, 101);  // 2*51 - 1
    EXPECT_EQ(result.name, "NACA 0012");

    // First and last points should be at trailing edge (x=1)
    EXPECT_NEAR(result.x(0), 1.0, 1e-6);
    EXPECT_NEAR(result.x(result.n - 1), 1.0, 1e-6);

    // Middle point should be at leading edge (x=0)
    int ile = result.n / 2;
    EXPECT_NEAR(result.x(ile), 0.0, 1e-6);

    // Symmetric airfoil should have y=0 at TE
    EXPECT_NEAR(result.y(0), 0.0, 0.01);  // Small gap due to thickness
    EXPECT_NEAR(result.y(result.n - 1), 0.0, 0.01);

    // Maximum thickness should be ~12% at x=0.3
    double max_thick = 0.0;
    for (int i = 0; i < result.n / 2; ++i) {
        // Upper surface (first half)
        if (result.y(i) > max_thick) {
            max_thick = result.y(i);
        }
    }
    // Half-thickness should be ~6%
    EXPECT_NEAR(max_thick, 0.06, 0.005);
}

// Test NACA 4-digit cambered airfoil
TEST(GeometryTest, Naca4Cambered) {
    // Generate NACA 2412
    auto result = naca4<double>(2412, 51);

    EXPECT_EQ(result.name, "NACA 2412");

    // For cambered airfoil, TE should have slight gap
    // Upper surface y > lower surface y at TE
    EXPECT_GT(result.y(0), result.y(result.n - 1));

    // Camber should make upper surface thicker than lower
    // Check at mid-chord
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
    // Upper surface should be higher due to camber
    EXPECT_GT(y_upper, -y_lower);
}

// Test NACA 5-digit airfoil generation
TEST(GeometryTest, Naca5Generation) {
    // Generate NACA 23012
    auto result = naca5<double>(23012, 51);

    EXPECT_EQ(result.n, 101);
    EXPECT_EQ(result.name, "NACA 23012");

    // Basic sanity checks
    EXPECT_NEAR(result.x(0), 1.0, 1e-6);
    EXPECT_NEAR(result.x(result.n - 1), 1.0, 1e-6);
}

// Test invalid NACA 5-digit
TEST(GeometryTest, Naca5Invalid) {
    // Invalid 5-digit code (not in table)
    auto result = naca5<double>(12345, 51);

    EXPECT_EQ(result.n, 0);  // Should return empty
}

// Test leading edge finding
TEST(GeometryTest, LeadingEdgeFind) {
    // Generate NACA 0012
    auto naca = naca4<double>(12, 51);

    // Compute arc length and splines
    VectorX<double> s, xp, yp;
    scalc(naca.x, naca.y, s);
    spline(naca.x, s, xp);
    spline(naca.y, s, yp);

    // Find leading edge
    double sle = lefind(naca.x, xp, naca.y, yp, s);

    // Evaluate LE position
    double xle = seval(sle, naca.x, xp, s);
    double yle = seval(sle, naca.y, yp, s);

    // LE should be at x=0, y=0 for symmetric airfoil
    EXPECT_NEAR(xle, 0.0, 1e-4);
    EXPECT_NEAR(yle, 0.0, 1e-4);
}

// Test geometric parameters
TEST(GeometryTest, GeometricParameters) {
    // Generate NACA 0012
    auto naca = naca4<double>(12, 51);

    // Compute arc length and splines
    VectorX<double> s, xp, yp;
    scalc(naca.x, naca.y, s);
    spline(naca.x, s, xp);
    spline(naca.y, s, yp);

    // Get geometric parameters
    auto gp = geopar(naca.x, xp, naca.y, yp, s);

    // Chord should be ~1.0
    EXPECT_NEAR(gp.chord, 1.0, 0.01);

    // LE at x=0
    EXPECT_NEAR(gp.xle, 0.0, 1e-4);

    // TE at x=1
    EXPECT_NEAR(gp.xte, 1.0, 1e-4);

    // LE radius should be positive
    EXPECT_GT(gp.radle, 0.0);
}

// Test normalization
TEST(GeometryTest, Normalization) {
    // Generate airfoil and shift/scale it
    auto naca = naca4<double>(12, 51);

    // Modify: shift and scale
    VectorX<double> x = naca.x * 2.0;  // Double the size
    VectorX<double> y = naca.y * 2.0;
    x.array() += 1.0;  // Shift x by 1
    y.array() += 0.5;  // Shift y by 0.5
    VectorX<double> s;

    // Normalize
    normalize(x, y, s);

    // After normalization, LE should be at origin
    VectorX<double> xp, yp;
    scalc(x, y, s);
    spline(x, s, xp);
    spline(y, s, yp);
    double sle = lefind(x, xp, y, yp, s);
    double xle = seval(sle, x, xp, s);
    double yle = seval(sle, y, yp, s);

    EXPECT_NEAR(xle, 0.0, 1e-4);
    EXPECT_NEAR(yle, 0.0, 1e-4);

    // Chord should be 1.0
    double xte = 0.5 * (x(0) + x(x.size() - 1));
    EXPECT_NEAR(xte, 1.0, 0.01);
}

// Test normal vector calculation
TEST(GeometryTest, NormalCalculation) {
    // Generate NACA 0012
    auto naca = naca4<double>(12, 51);

    // Compute arc length
    VectorX<double> s;
    scalc(naca.x, naca.y, s);

    // Compute normals
    VectorX<double> nx, ny;
    ncalc(naca.x, naca.y, s, nx, ny);

    // All normals should be unit vectors
    for (int i = 0; i < naca.n; ++i) {
        double mag = std::sqrt(nx(i) * nx(i) + ny(i) * ny(i));
        EXPECT_NEAR(mag, 1.0, 1e-10);
    }

    // On upper surface (first half), normals should point up (ny > 0)
    // except near TE where they might point backward
    int imid = naca.n / 4;  // Quarter way around (mid upper surface)
    EXPECT_GT(ny(imid), 0.0);

    // On lower surface (second half), normals should point down (ny < 0)
    int imid_lower = 3 * naca.n / 4;
    EXPECT_LT(ny(imid_lower), 0.0);
}

// Test trailing edge calculation for sharp TE
TEST(GeometryTest, TrailingEdgeSharp) {
    // Generate NACA 0012 (has sharp TE in our generation)
    auto naca = naca4<double>(12, 51);

    // Compute splines
    VectorX<double> s, xp, yp;
    scalc(naca.x, naca.y, s);
    segspl(naca.x, s, xp);
    segspl(naca.y, s, yp);

    // Get TE parameters
    auto te = tecalc(naca.x, xp, naca.y, yp, 1.0);

    // NACA 4-digit has small TE gap but not exactly zero
    // The thickness polynomial gives finite thickness at x=1
    EXPECT_LT(te.dste, 0.01);  // Small TE gap

    // For nearly sharp TE, direction factors should be ~(1, 0)
    if (te.sharp) {
        EXPECT_NEAR(te.scs, 1.0, 1e-6);
        EXPECT_NEAR(te.sds, 0.0, 1e-6);
    }
}

// Test trailing edge calculation for blunt TE
TEST(GeometryTest, TrailingEdgeBlunt) {
    // Create an airfoil with blunt TE by modifying NACA
    auto naca = naca4<double>(12, 51);

    // Make TE blunt by shifting the last point
    naca.y(naca.n - 1) -= 0.02;  // Move lower TE down

    // Compute splines
    VectorX<double> s, xp, yp;
    scalc(naca.x, naca.y, s);
    segspl(naca.x, s, xp);
    segspl(naca.y, s, yp);

    // Get TE parameters
    auto te = tecalc(naca.x, xp, naca.y, yp, 1.0);

    // Now TE should not be sharp
    EXPECT_FALSE(te.sharp);
    EXPECT_GT(te.dste, 0.01);

    // Direction factors should be normalized properly
    double mag = std::sqrt(te.scs * te.scs + te.sds * te.sds);
    // scs and sds are not unit normalized, they're projections
    // Just check they're finite
    EXPECT_TRUE(std::isfinite(te.scs));
    EXPECT_TRUE(std::isfinite(te.sds));
}

// Test circle geometry (for curvature validation)
TEST(GeometryTest, CircleNormals) {
    // Create a circle
    int n = 65;
    VectorX<double> x(n), y(n);
    double radius = 0.5;

    for (int i = 0; i < n; ++i) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / (n - 1);
        x(i) = 0.5 + radius * std::cos(theta);  // Center at (0.5, 0)
        y(i) = radius * std::sin(theta);
    }

    // Compute arc length and normals
    VectorX<double> s, nx, ny;
    scalc(x, y, s);
    ncalc(x, y, s, nx, ny);

    // For a circle, normals should point radially outward
    // At each point, normal should be parallel to (x-center, y-center)
    for (int i = 1; i < n - 1; ++i) {
        double rx = x(i) - 0.5;
        double ry = y(i);
        double rmag = std::sqrt(rx * rx + ry * ry);
        rx /= rmag;
        ry /= rmag;

        // Normal should be close to radial direction (within sign)
        double dot = nx(i) * rx + ny(i) * ry;
        EXPECT_GT(std::abs(dot), 0.95);  // Allow some tolerance
    }
}

}  // namespace
}  // namespace xfoil
