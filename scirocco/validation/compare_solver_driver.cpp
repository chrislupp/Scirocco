// Standalone comparison driver; does not alter the production solver.
// Build against the existing Release libscirocco.a (see compare_solvers.py).
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/forces.h"
#include "scirocco/scirocco.h"

int main(int argc, char** argv) {
    using namespace scirocco;
    if (argc < 7 || argc > 9) return 2;
    std::ifstream input(argv[1]);
    std::vector<std::array<double, 5>> nodes;
    std::array<double, 5> row{};
    while (input >> row[0] >> row[1] >> row[2] >> row[3] >> row[4])
        nodes.push_back(row);
    if (nodes.size() < 3) return 3;
    Scirocco foil("0012", static_cast<int>(nodes.size()) - 1);
    auto& state = foil.GetState();
    for (size_t i = 0; i < nodes.size(); ++i) {
        state.foil.x.col(i) << nodes[i][0], nodes[i][1];
        state.foil.s(i) = nodes[i][2];
        state.foil.t.col(i) << nodes[i][3], nodes[i][4];
    }
    state.param.verb = 1;
    state.param.ncrit = 9;
    state.param.niglob = 200;
    if (argc >= 8) state.param.use_ptc = std::stoi(argv[7]) != 0;
    if (argc >= 9) state.param.scale_coupled_system = std::stoi(argv[8]) >= 1;
    if (argc >= 9) state.param.use_trust_region = std::stoi(argv[8]) >= 2;
    if (argc >= 9) state.param.trust_update_topology = std::stoi(argv[8]) >= 3;
    foil.SetOper(std::stod(argv[2]), -1, std::stod(argv[3]),
                 std::stod(argv[4]));
    const std::string prefix(argv[5]);
    std::ofstream summary(prefix + "/scirocco.json");
    summary << std::setprecision(17);
    auto number = [&](const char* key, double value) {
        summary << ",\"" << key << "\":";
        if (std::isfinite(value))
            summary << value;
        else
            summary << "null";
    };
    try {
        const bool ok = foil.Solve();
        summary << "{\"converged\":" << (ok ? "true" : "false");
        number("cl", foil.GetCl());
        number("cd", foil.GetCd());
        number("cm", foil.GetCm());
        number("cdf", foil.GetCdf());
        number("cdpi", state.post.cdpi);
        number("rtol", state.param.rtol);
        number("ptc_enabled", state.param.use_ptc ? 1 : 0);
        number("trust_accepted", state.glob.trust_accepted);
        number("trust_rejected", state.glob.trust_rejected);
        number("residual", state.oper.viscous ? state.glob.R.norm() : 0.0);
        if (state.oper.viscous) {
            number("stored_xtr_upper", state.vsol.Xt(1, 1));
            number("stored_xtr_lower", state.vsol.Xt(0, 1));
            // The production solve currently leaves the stored x coordinates
            // unset. Recover the actual transition from its final residual on
            // a copy, without modifying the solution used for comparisons.
            for (int si = 0; si < 2; ++si) {
                double xtransition = 1.0;
                const auto& indices = state.vsol.Is[si];
                for (size_t k = 1; k < indices.size(); ++k) {
                    const int i0 = indices[k - 1], i1 = indices[k];
                    if (!state.vsol.turb[i0] && state.vsol.turb[i1]) {
                        auto copy = state;
                        auto param = BuildParam(copy, si);
                        param.simi = false;
                        param.turb = false;
                        Matrix42d u;
                        u.col(0) = state.glob.U.col(i0);
                        u.col(1) = state.glob.U.col(i1);
                        Vector2d xi(state.isol.xi(i0), state.isol.xi(i1));
                        ResidualTransition(copy, param, xi, u,
                                           Vector2d::Zero());
                        const double weight =
                            (copy.vsol.xt - xi(0)) / (xi(1) - xi(0));
                        xtransition = (1 - weight) * state.foil.x(0, i0) +
                                      weight * state.foil.x(0, i1);
                        break;
                    }
                }
                number(si == 0 ? "xtr_lower" : "xtr_upper", xtransition);
            }
            GetDistributions(state);
        }
        std::ofstream surface(prefix + "/scirocco_surface.csv");
        surface << "x,y,cp,cf,theta,dstar\n" << std::setprecision(17);
        for (size_t i = 0; i < nodes.size(); ++i) {
            surface << state.foil.x(0, i) << ',' << state.foil.x(1, i) << ','
                    << state.post.cp(i);
            if (state.oper.viscous)
                surface << ',' << state.post.cf(i) << ',' << state.post.th(i)
                        << ',' << state.post.ds(i);
            else
                surface << ",0,0,0";
            surface << '\n';
        }
        if (ok && std::stoi(argv[6])) {
            const auto gradient = foil.ComputeSensitivities();
            number("dcl_dalpha", gradient.dcl_dalpha);
            number("dcd_dalpha", gradient.dcd_dalpha);
            number("dcm_dalpha", gradient.dcm_dalpha);
        }
        summary << "}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
