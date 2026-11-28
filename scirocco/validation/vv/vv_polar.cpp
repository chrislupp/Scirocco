/* Scirocco polar on prescribed panel nodes, for the verification and
 * validation pages of the documentation (see README.md).
 *
 * Usage: vv_polar NODES RE MA NCRIT ALPHAS DIST_ALPHAS OUT.json [TRIPS
 *                 [mfoil] [captured]]
 *   NODES        one node per line: x y s tx ty (lower TE -> LE -> upper TE),
 *                as mfoil panels it; Scirocco uses the nodes, arc length
 *                and tangents unchanged. With only x y per line, Scirocco
 *                panels the coordinates itself, with as many nodes.
 *   RE           Reynolds number (0 for an inviscid solve)
 *   ALPHAS       comma-separated angles of attack [deg], solved in order on
 *                one object, so each point warm-starts from the previous one
 *   DIST_ALPHAS  comma-separated subset of ALPHAS whose surface
 *                distributions are written ("-" for none)
 *   TRIPS        optional forced transition: XTRIP_UPPER,XTRIP_LOWER (x/c),
 *                or "-" for none
 *   mfoil        optional: mfoil's closure (Param::cf_laminar_floor = false)
 *   captured     optional: TransitionModel::Captured
 */
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/forces.h"
#include "scirocco/scirocco.h"

using namespace scirocco;

namespace {

std::vector<double> ParseList(const std::string& text) {
    std::vector<double> out;
    if (text == "-") return out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) out.push_back(std::stod(item));
    return out;
}

std::string Num(double v) {
    if (!std::isfinite(v)) return "null";
    std::ostringstream os;
    os << std::setprecision(10) << v;
    return os.str();
}

/* Transition x/c on surface si (0 lower, 1 upper); 1 where the surface
 * stays laminar. */
double TransitionX(const FoilState& state, int si) {
    return state.vsol.turb[state.vsol.Is[si].back()] ? state.vsol.Xt(si, 1)
                                                     : 1.0;
}

}  // namespace

int main(int argc, char** argv) {
    bool mfoil_closure = false, captured = false, usage = argc < 8;
    for (int k = 9; k < argc; ++k) {
        const std::string flag = argv[k];
        if (flag == "mfoil")
            mfoil_closure = true;
        else if (flag == "captured")
            captured = true;
        else
            usage = true;
    }
    if (usage) {
        std::cerr << "usage: vv_polar NODES RE MA NCRIT ALPHAS DIST_ALPHAS "
                     "OUT.json [XTRIP_UPPER,XTRIP_LOWER|- [mfoil] "
                     "[captured]]\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    std::vector<std::vector<double>> nodes;
    for (std::string line; std::getline(input, line);) {
        std::istringstream fields(line);
        std::vector<double> row;
        for (double v; fields >> v;) row.push_back(v);
        if (!row.empty()) nodes.push_back(row);
    }
    if (nodes.size() < 3) return 3;
    const bool panelled = nodes[0].size() == 5;
    for (const auto& row : nodes)
        if (row.size() != (panelled ? 5U : 2U)) return 3;
    const double re = std::stod(argv[2]), ma = std::stod(argv[3]);
    const auto alphas = ParseList(argv[5]);
    const auto dist = ParseList(argv[6]);
    const std::set<double> dist_set(dist.begin(), dist.end());

    const int npanel = static_cast<int>(nodes.size()) - 1;
    Eigen::MatrixXd coords(2, nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i)
        coords.col(i) << nodes[i][0], nodes[i][1];
    Scirocco foil =
        panelled ? Scirocco("0012", npanel) : Scirocco(coords, npanel);
    auto& state = foil.GetState();
    for (size_t i = 0; panelled && i < nodes.size(); ++i) {
        state.foil.x.col(i) = coords.col(i);
        state.foil.s(i) = nodes[i][2];
        state.foil.t.col(i) << nodes[i][3], nodes[i][4];
    }
    state.param.verb = 0;
    state.param.warm_start = true; /* a polar sweep, as a user would run it */
    state.param.ncrit = std::stod(argv[4]);
    if (mfoil_closure) state.param.cf_laminar_floor = false;
    if (captured) state.param.transition = TransitionModel::Captured;
    if (argc >= 9 && std::string(argv[8]) != "-") {
        const auto trip = ParseList(argv[8]);
        if (trip.size() != 2) return 2;
        state.param.xtrip_upper = trip[0];
        state.param.xtrip_lower = trip[1];
    }

    std::ofstream out(argv[7]);
    out << "{\"points\":[";
    std::ostringstream dists;
    bool first_dist = true;
    for (size_t k = 0; k < alphas.size(); ++k) {
        const double alpha = alphas[k];
        foil.SetOper(alpha, -1, re, ma);
        bool ok = false;
        try {
            ok = foil.Solve();
        } catch (const std::exception& e) {
            std::cerr << "alpha " << alpha << ": " << e.what() << '\n';
        }
        auto& st = foil.GetState();
        out << (k ? "," : "") << "{\"alpha\":" << Num(alpha)
            << ",\"converged\":" << (ok ? "true" : "false");
        if (ok) {
            out << ",\"cl\":" << Num(foil.GetCl())
                << ",\"cd\":" << Num(foil.GetCd())
                << ",\"cdf\":" << Num(foil.GetCdf())
                << ",\"cm\":" << Num(foil.GetCm());
            if (re > 0)
                out << ",\"xtr_upper\":" << Num(TransitionX(st, 1))
                    << ",\"xtr_lower\":" << Num(TransitionX(st, 0));
        }
        out << "}";
        if (ok && dist_set.count(alpha)) {
            if (re > 0) GetDistributions(st);
            dists << (first_dist ? "" : ",") << "{\"alpha\":" << Num(alpha);
            first_dist = false;
            auto column = [&](const char* key, auto value) {
                dists << ",\"" << key << "\":[";
                for (int i = 0; i < st.foil.N; ++i)
                    dists << (i ? "," : "") << Num(value(i));
                dists << "]";
            };
            column("x", [&](int i) { return st.foil.x(0, i); });
            column("y", [&](int i) { return st.foil.x(1, i); });
            column("cp", [&](int i) { return st.post.cp(i); });
            if (re > 0) {
                column("cf", [&](int i) { return st.post.cf(i); });
                column("theta", [&](int i) { return st.post.th(i); });
                column("dstar", [&](int i) { return st.post.ds(i); });
            }
            dists << "}";
        }
    }
    out << "],\"distributions\":[" << dists.str() << "]}\n";
    return 0;
}
