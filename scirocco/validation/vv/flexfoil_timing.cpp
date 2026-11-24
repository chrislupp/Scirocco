/* Scirocco viscous solve times on prescribed panel nodes, for the flexfoil
 * comparison on the documentation's Performance page (see README.md).
 *
 * Usage: flexfoil_timing NODES RE ALPHA REPS
 *   NODES  one node per line: x y s tx ty, as vv_polar reads them
 *
 * Prints the median cold-start time [ms], the median warm-start time [ms]
 * (a solve at ALPHA after one at ALPHA - 0.5 on the same object), the cold
 * start's Newton steps, whether every solve converged, Cl and Cd.
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "scirocco/scirocco.h"

using namespace scirocco;
using Clock = std::chrono::steady_clock;

namespace {

std::vector<std::array<double, 5>> nodes;

Scirocco Make() {
    Scirocco foil("0012", static_cast<int>(nodes.size()) - 1);
    auto& st = foil.GetState();
    for (size_t i = 0; i < nodes.size(); ++i) {
        st.foil.x.col(i) << nodes[i][0], nodes[i][1];
        st.foil.s(i) = nodes[i][2];
        st.foil.t.col(i) << nodes[i][3], nodes[i][4];
    }
    st.param.verb = 0;
    st.param.ncrit = 9.0;
    return foil;
}

double Median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

double Ms(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: flexfoil_timing NODES RE ALPHA REPS\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::array<double, 5> row{};
    while (in >> row[0] >> row[1] >> row[2] >> row[3] >> row[4]) nodes.push_back(row);
    const double re = std::stod(argv[2]), alpha = std::stod(argv[3]);
    const int reps = std::stoi(argv[4]);

    std::vector<double> cold, warm;
    bool ok = true;
    int steps = 0;
    double cl = 0, cd = 0;
    for (int k = 0; k < reps; ++k) {
        auto foil = Make();
        foil.SetOper(alpha, -1, re, 0.0);
        foil.GetState().glob.iteration_observer = [&](const FoilState&, int i) { steps = i + 1; };
        auto t0 = Clock::now();
        ok = foil.Solve() && ok;
        cold.push_back(Ms(t0));
        cl = foil.GetCl();
        cd = foil.GetCd();

        auto next = Make();
        next.SetOper(alpha - 0.5, -1, re, 0.0);
        next.Solve();
        next.SetOper(alpha, -1, re, 0.0);
        t0 = Clock::now();
        ok = next.Solve() && ok;
        warm.push_back(Ms(t0));
    }
    std::cout << Median(cold) << " " << Median(warm) << " " << steps << " " << ok << " "
              << cl << " " << cd << "\n";
    return 0;
}
