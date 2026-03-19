/*
    cppfoil: C++ Airfoil Analysis Library

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
#include "cppfoil/core/spline.h"

#include <gtest/gtest.h>

namespace cppfoil {
namespace {

TEST(SplineTest, QuadSegReturnsCorrectWeights) {
    VectorXd xq, wq;
    QuadSeg(xq, wq);

    EXPECT_EQ(xq.size(), 5);
    EXPECT_EQ(wq.size(), 5);

    // Weights should sum to 1
    EXPECT_NEAR(wq.sum(), 1.0, 1e-10);
}

TEST(SplineTest, Spline2DCreateBasic) {
    // Simple line from (0,0) to (1,1)
    MatrixXd points(2, 3);
    points << 0.0, 0.5, 1.0, 0.0, 0.5, 1.0;

    Spline2D spline = Spline2DCreate(points);

    // Midpoint should be approximately (0.5, 0.5)
    double smax = spline.X.s(spline.X.s.size() - 1);
    VectorXd s_mid(1);
    s_mid << smax / 2.0;
    MatrixXd result = SplineEval(spline, s_mid);

    EXPECT_NEAR(result(0, 0), 0.5, 0.1);
    EXPECT_NEAR(result(1, 0), 0.5, 0.1);
}

TEST(SplineTest, SplineTangentNonZero) {
    // Circle arc
    MatrixXd points(2, 5);
    for (int i = 0; i < 5; ++i) {
        double angle = i * M_PI / 4.0;
        points(0, i) = std::cos(angle);
        points(1, i) = std::sin(angle);
    }

    Spline2D spline = Spline2DCreate(points);

    VectorXd s_vals =
        VectorXd::LinSpaced(3, 0, spline.X.s(spline.X.s.size() - 1));
    MatrixXd tangents = SplineTangent(spline, s_vals);

    // Tangents should be non-zero
    for (int i = 0; i < tangents.cols(); ++i) {
        double norm = tangents.col(i).norm();
        EXPECT_GT(norm, 0.1);
    }
}

}  // namespace
}  // namespace cppfoil
