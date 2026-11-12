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
/* Guards the cost of the full viscous gradient.
 *
 * The adjoint for cl, cd, cm and cdpi with respect to every node, alpha, Re
 * and Ma should cost a fraction of the viscous solve it differentiates: about
 * 0.2 of it on NACA 0012 at 149 panels. Losing the batched and reversed
 * chains of the gradient (as happened when PRs #26-#28 never reached
 * develop) pushed that to about 1.2. Comparing two phases of one run keeps
 * the check independent of the machine's speed. Built for optimized,
 * uninstrumented builds only. */
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>

#include "scirocco/scirocco.h"

namespace {
using Clock = std::chrono::steady_clock;

double Seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}
}  // namespace

TEST(GradientCost, FullGradientCostsLessThanTheViscousSolve) {
    constexpr int kRuns = 3;
    constexpr double kMaxRatio = 0.6;
    double solve = 1e30;
    double gradient = 1e30;
    for (int run = 0; run < kRuns; ++run) {
        scirocco::Scirocco foil("0012", 149);
        foil.SetOper(5.0, -1.0, 1.0e6, 0.0);
        auto start = Clock::now();
        ASSERT_TRUE(foil.Solve());
        solve = std::min(solve, Seconds(start));
        start = Clock::now();
        static_cast<void>(foil.ComputeSensitivities());
        gradient = std::min(gradient, Seconds(start));
    }
    EXPECT_LT(gradient / solve, kMaxRatio)
        << "gradient " << gradient * 1e3 << " ms, viscous solve " << solve * 1e3
        << " ms";
}
