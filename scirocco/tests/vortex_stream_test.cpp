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
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "scirocco/core/panel_influence_deriv.h"

using namespace scirocco;
TEST(VortexStreamTest, HighPrecisionIntegralAndCoordinateDerivatives) {
    std::ifstream input(VORTEX_STREAM_REFERENCE_FILE);
    ASSERT_TRUE(input.good());
    std::string line;
    int cases = 0;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::array<double, 20> ref{};
        for (double& value : ref)
            ASSERT_TRUE(static_cast<bool>(fields >> value));
        SCOPED_TRACE(cases++);
        Matrix2d panel;
        panel << ref[0], ref[2], ref[1], ref[3];
        Vector2d xi(ref[4], ref[5]);
        double a, b;
        PanelLinVortexStream(panel, xi, a, b);
        const auto deriv = PanelLinVortexStreamDeriv(panel, xi);
        const double length = (panel.col(1) - panel.col(0)).norm();
        EXPECT_NEAR(a, ref[6], 5e-16 * length + 2e-13 * std::abs(ref[6]));
        EXPECT_NEAR(b, ref[7], 5e-16 * length + 2e-13 * std::abs(ref[7]));
        EXPECT_DOUBLE_EQ(a, deriv.a);
        EXPECT_DOUBLE_EQ(b, deriv.b);
        for (int k = 0; k < 6; ++k) {
            SCOPED_TRACE(k);
            const double da = k < 4 ? deriv.a_Xj(k) : deriv.a_xi(k - 4);
            const double db = k < 4 ? deriv.b_Xj(k) : deriv.b_xi(k - 4);
            EXPECT_NEAR(da, ref[8 + k], 5e-14 + 2e-13 * std::abs(ref[8 + k]));
            EXPECT_NEAR(db, ref[14 + k], 5e-14 + 2e-13 * std::abs(ref[14 + k]));
        }
    }
    EXPECT_EQ(cases, 7);
}

TEST(SourceStreamTest, HighPrecisionIntegralAndCoordinateDerivatives) {
    std::ifstream input(SOURCE_STREAM_REFERENCE_FILE);
    ASSERT_TRUE(input.good());
    std::string line;
    int cases = 0;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::array<double, 27> ref{};
        for (double& value : ref)
            ASSERT_TRUE(static_cast<bool>(fields >> value));
        SCOPED_TRACE(cases++);
        Matrix2d panel;
        panel << ref[0], ref[2], ref[1], ref[3];
        Vector2d xi(ref[4], ref[5]);
        double a, b;
        PanelLinSourceStream(panel, xi, a, b);
        const auto lin = PanelLinSourceStreamDeriv(panel, xi);
        const auto con = PanelConstSourceStreamDeriv(panel, xi);
        const std::array<double, 3> values{PanelConstSourceStream(panel, xi), a,
                                           b};
        const double length = (panel.col(1) - panel.col(0)).norm();
        const std::array<Vector4d, 3> endpoint{con.val_Xj, lin.a_Xj, lin.b_Xj};
        const std::array<Vector2d, 3> point{con.val_xi, lin.a_xi, lin.b_xi};
        EXPECT_DOUBLE_EQ(values[0], con.val);
        EXPECT_DOUBLE_EQ(a, lin.a);
        EXPECT_DOUBLE_EQ(b, lin.b);
        for (int j = 0; j < 3; ++j) {
            SCOPED_TRACE(j);
            EXPECT_NEAR(values[j], ref[6 + j],
                        5e-16 * length + 2e-13 * std::abs(ref[6 + j]));
            for (int k = 0; k < 6; ++k) {
                SCOPED_TRACE(k);
                const double got = k < 4 ? endpoint[j](k) : point[j](k - 4);
                const double expected = ref[9 + 6 * j + k];
                EXPECT_NEAR(got, expected, 5e-14 + 2e-13 * std::abs(expected));
            }
        }
    }
    EXPECT_EQ(cases, 11);
}
