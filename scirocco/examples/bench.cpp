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
*/

/**
 * @file bench.cpp
 * @brief Wall-clock timings for the numbers quoted in the documentation.
 *
 * Each case is a cold start: a fresh Scirocco object, paneling, operating
 * point, and solve. Phases are reported separately as medians over the
 * repetitions.
 *
 * Usage: bench [repetitions]   (default 11)
 */

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "scirocco/scirocco.h"

namespace {

using Clock = std::chrono::steady_clock;

/// Milliseconds elapsed since @p start.
double Since(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start)
        .count();
}

/// Phase timings for one cold start, in milliseconds.
struct Phases {
    double paneling = 0.0;
    double solve = 0.0;
    double gradient = 0.0;
};

/// Median of each phase over @p reps runs of @p fn.
Phases Median(int reps, const std::function<Phases()>& fn) {
    static_cast<void>(fn());  // warm caches and the allocator
    std::vector<double> p, s, g;
    for (int i = 0; i < reps; ++i) {
        const Phases t = fn();
        p.push_back(t.paneling);
        s.push_back(t.solve);
        g.push_back(t.gradient);
    }
    auto med = [](std::vector<double>& v) {
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    };
    return {med(p), med(s), med(g)};
}

void Require(bool ok, const char* what) {
    if (!ok) {
        throw std::runtime_error(std::string("did not converge: ") + what);
    }
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 11;
    if (argc > 1) {
        char* end = nullptr;
        const long n = std::strtol(argv[1], &end, 10);
        if (end == argv[1] || *end != '\0' || n < 1 || n > 1000) {
            std::cerr << "usage: bench [repetitions 1-1000]\n";
            return 2;
        }
        reps = static_cast<int>(n);
    }

    std::cout << "Scirocco benchmark: NACA 0012, alpha 5 deg, Re 1e6.\n"
              << "Median of " << reps << " cold starts; milliseconds.\n\n"
              << "  panels   paneling   inviscid    viscous   gradient\n";

    try {
        for (const int npanel : {99, 149, 199}) {
            // One cold start: panel, solve at (alpha, Re), optionally
            // differentiate. Each phase is timed separately.
            auto run = [npanel](double re, bool gradient) {
                Phases t;
                auto start = Clock::now();
                scirocco::Scirocco foil("0012", npanel);
                t.paneling = Since(start);
                foil.SetOper(5.0, -1.0, re, 0.0);
                start = Clock::now();
                Require(foil.Solve(), re > 0.0 ? "viscous" : "inviscid");
                t.solve = Since(start);
                if (gradient) {
                    start = Clock::now();
                    static_cast<void>(foil.ComputeSensitivities());
                    t.gradient = Since(start);
                }
                return t;
            };
            const Phases inviscid =
                Median(reps, [&run] { return run(0.0, false); });
            const Phases viscous =
                Median(reps, [&run] { return run(1.0e6, true); });
            std::cout << std::fixed << std::setprecision(1) << std::setw(8)
                      << npanel << std::setw(11) << viscous.paneling
                      << std::setw(11) << inviscid.solve << std::setw(11)
                      << viscous.solve << std::setw(11) << viscous.gradient
                      << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "bench failed: " << error.what() << '\n';
        return 1;
    }
    std::cout << "\n  Columns are separate phases. gradient is the "
                 "adjoint for cl, cd, cm\n  w.r.t. every node coordinate "
                 "plus alpha, Re, and Ma.\n";
    return 0;
}
