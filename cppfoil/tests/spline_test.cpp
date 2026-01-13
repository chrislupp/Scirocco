// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include <gtest/gtest.h>
#include "cppfoil/core/spline.h"

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
    points << 0.0, 0.5, 1.0,
              0.0, 0.5, 1.0;

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

    VectorXd s_vals = VectorXd::LinSpaced(3, 0, spline.X.s(spline.X.s.size() - 1));
    MatrixXd tangents = SplineTangent(spline, s_vals);

    // Tangents should be non-zero
    for (int i = 0; i < tangents.cols(); ++i) {
        double norm = tangents.col(i).norm();
        EXPECT_GT(norm, 0.1);
    }
}

}  // namespace
}  // namespace cppfoil
