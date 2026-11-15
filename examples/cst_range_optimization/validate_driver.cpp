// Independent cold-start function-value checks, with no optimization server.
// c++ -O3 -std=c++17 -I scirocco/include -isystem /opt/homebrew/include/eigen3 \
//   examples/cst_range_optimization/validate_driver.cpp BUILD/libscirocco.a -o DRIVER
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "scirocco/scirocco.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/paneling.h"

int main(int argc, char** argv) {
    using namespace scirocco;
    if (argc != 7 && argc != 8) return 2; // xy-file, alpha, Re, Mach, Ncrit, result-file, [verbosity]
    std::ifstream input(argv[1]);
    std::vector<Vector2d> nodes;
    double x, y;
    while (input >> x >> y) nodes.emplace_back(x, y);
    if (nodes.size() < 5) return 3;
    try {
        Scirocco foil("0012", static_cast<int>(nodes.size()) - 1);
        auto& state = foil.GetState();
        for (size_t i = 0; i < nodes.size(); ++i) state.foil.x.col(i) = nodes[i];
        RefreshPanelGeometry(state.foil);
        state.param.verb = argc == 8 ? std::stoi(argv[7]) : 0;
        if (state.param.verb > 0) std::cout << std::unitbuf;
        state.param.ncrit = std::stod(argv[5]);
        state.param.niglob = 200;
        foil.SetOper(std::stod(argv[2]), -1, std::stod(argv[3]), std::stod(argv[4]));
        const bool ok = foil.Solve();
        const VectorXd ue = state.glob.U.row(3).transpose();
        const VectorXd coupling = ue - GetUeInv(state) - state.vsol.ue_m *
            state.glob.U.row(1).transpose().cwiseProduct(ue);
        const double residual = std::sqrt(state.glob.R.squaredNorm() + coupling.squaredNorm());
        if (!std::isfinite(residual) || !std::isfinite(foil.GetCd())) return 4;
        std::ofstream output(argv[6]);
        output << std::setprecision(17) << "{\"converged\":" << (ok ? "true" : "false")
               << ",\"alpha\":" << state.oper.alpha << ",\"cl\":" << foil.GetCl()
               << ",\"cd\":" << foil.GetCd() << ",\"cm\":" << foil.GetCm()
               << ",\"full_residual\":" << residual << "}\n";
        return ok ? 0 : 5;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
