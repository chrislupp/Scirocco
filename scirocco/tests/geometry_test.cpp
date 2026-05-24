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
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "scirocco/core/geometry.h"

#include <gtest/gtest.h>

namespace scirocco {
namespace {

TEST(GeometryTest, NacaPoints4Digit) {
    Geom geom;
    NacaPoints(geom, "0012");

    EXPECT_EQ(geom.name, "NACA 0012");
    EXPECT_GT(geom.npoint, 0);
    EXPECT_NEAR(geom.chord, 1.0, 1e-6);

    // Points should form a closed loop (first and last x should be similar)
    EXPECT_NEAR(geom.xpoint(0, 0), geom.xpoint(0, geom.npoint - 1), 0.1);
}

TEST(GeometryTest, NacaPointsCambered) {
    Geom geom;
    NacaPoints(geom, "2412");

    EXPECT_EQ(geom.name, "NACA 2412");
    EXPECT_GT(geom.npoint, 0);

    // Cambered airfoil should have y != 0 at max camber location
    double max_camber = 0.0;
    for (int i = 0; i < geom.npoint; ++i) {
        max_camber = std::max(max_camber, std::abs(geom.xpoint(1, i)));
    }
    EXPECT_GT(max_camber, 0.01);  // Should have some thickness/camber
}

TEST(GeometryTest, SetCoords) {
    Geom geom;
    MatrixXd coords(2, 3);
    coords << 1.0, 0.0, 1.0, 0.0, 0.0, 0.0;

    SetCoords(geom, coords);

    EXPECT_EQ(geom.npoint, 3);
    EXPECT_EQ(geom.name, "custom");
    EXPECT_NEAR(geom.chord, 1.0, 1e-6);
}

TEST(GeometryTest, TEInfoClosed) {
    // Create closed airfoil (TE gap = 0)
    MatrixXd x(2, 5);
    x << 1.0, 0.5, 0.0, 0.5, 1.0, 0.0, -0.1, 0.0, 0.1, 0.0;

    TEInfo info = GetTEInfo(x);

    EXPECT_NEAR(info.hTE, 0.0, 1e-6);
}

TEST(GeometryTest, TEInfoBlunt) {
    // Create blunt TE airfoil
    MatrixXd x(2, 5);
    x << 1.0, 0.5, 0.0, 0.5, 1.0, -0.01, -0.1, 0.0, 0.1, 0.01;

    TEInfo info = GetTEInfo(x);

    EXPECT_GT(info.hTE, 0.0);
    EXPECT_LT(info.hTE, 0.1);
}

}  // namespace
}  // namespace scirocco
