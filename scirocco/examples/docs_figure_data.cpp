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
 * @file docs_figure_data.cpp
 * @brief Writes the solver data plotted on the documentation site.
 *
 * Usage: docs_figure_data <output-dir>
 *
 * Produces cp.json, polar.json, convergence.json, and panels.json. Regenerate
 * them whenever solver results change:
 *
 *   docs_figure_data docs/src/data
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "scirocco/scirocco.h"

namespace {

using scirocco::Scirocco;

/// Minimal JSON object writer: numbers, bools, and number arrays.
class Json {
   public:
    Json& Num(const std::string& key, double v) { return Raw(key, Fmt(v)); }
    Json& Bool(const std::string& key, bool v) {
        return Raw(key, v ? "true" : "false");
    }
    Json& Str(const std::string& key, const std::string& v) {
        return Raw(key, "\"" + v + "\"");
    }
    template <typename Vec>
    Json& Arr(const std::string& key, const Vec& v) {
        std::string s = "[";
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(v.size()); ++i) {
            s += (i ? "," : "") + Fmt(v[i]);
        }
        return Raw(key, s + "]");
    }
    Json& Obj(const std::string& key, const Json& o) {
        return Raw(key, o.Str());
    }
    Json& Objs(const std::string& key, const std::vector<Json>& v) {
        std::string s = "[";
        for (std::size_t i = 0; i < v.size(); ++i) {
            s += (i ? "," : "") + v[i].Str();
        }
        return Raw(key, s + "]");
    }
    [[nodiscard]] std::string Str() const { return "{" + body_ + "}"; }

   private:
    static std::string Fmt(double v) {
        if (!std::isfinite(v)) {
            return "null";
        }
        std::ostringstream os;
        os.precision(7);
        os << v;
        return os.str();
    }
    Json& Raw(const std::string& key, const std::string& value) {
        body_ += (body_.empty() ? "\"" : ",\"") + key + "\":" + value;
        return *this;
    }
    std::string body_;
};

void Write(const std::string& dir, const std::string& name, const Json& j) {
    std::ofstream out(dir + "/" + name);
    if (!out) {
        throw std::runtime_error("cannot write " + dir + "/" + name);
    }
    out << j.Str() << '\n';
    std::cout << "wrote " << dir << "/" << name << '\n';
}

/// Solve with verbose output captured; returns the |R| history.
std::pair<bool, std::vector<double>> SolveWithHistory(Scirocco& foil) {
    foil.GetState().param.verb = 1;
    std::ostringstream log;
    std::streambuf* old = std::cout.rdbuf(log.rdbuf());
    bool ok = false;
    try {
        ok = foil.Solve();
    } catch (...) {
        std::cout.rdbuf(old);
        throw;
    }
    std::cout.rdbuf(old);
    foil.GetState().param.verb = 0;

    std::vector<double> history;
    std::istringstream lines(log.str());
    std::string line;
    const std::string tag = "|R| = ";
    while (std::getline(lines, line)) {
        if (line.find("Iteration ") == std::string::npos) {
            continue;
        }
        const auto pos = line.find(tag);
        if (pos != std::string::npos) {
            history.push_back(std::stod(line.substr(pos + tag.size())));
        }
    }
    return {ok, history};
}

Json Coords(const Eigen::MatrixXd& x) {
    Json j;
    j.Arr("x", Eigen::VectorXd(x.row(0).transpose()));
    j.Arr("z", Eigen::VectorXd(x.row(1).transpose()));
    return j;
}

}  // namespace

int main(int argc, char** argv) {
    constexpr double kRe = 1.0e6;
    try {
        if (argc != 2) {
            std::cerr << "usage: docs_figure_data <output-dir>\n";
            return 2;
        }
        const std::string dir = argv[1];

        // Pressure and boundary layer: NACA 2412, alpha 4, Re 1e6.
        {
            Scirocco foil("2412", 160);
            foil.SetOper(4.0, -1.0, kRe, 0.0);
            if (!foil.Solve()) {
                throw std::runtime_error("cp case did not converge");
            }
            const auto& post = foil.GetState().post;
            const Eigen::Index n = foil.GetFoilCoords().cols();
            Json j;
            j.Str("airfoil", "NACA 2412")
                .Num("alpha", 4.0)
                .Num("Re", kRe)
                .Num("cl", foil.GetCl())
                .Num("cd", foil.GetCd())
                .Num("cm", foil.GetCm())
                .Obj("coords", Coords(foil.GetFoilCoords()))
                .Arr("cp", post.cp)
                .Arr("cpi", post.cpi)
                .Arr("cf", Eigen::VectorXd(post.cf.head(n)))
                .Arr("hk", Eigen::VectorXd(post.Hk.head(n)))
                .Arr("theta", Eigen::VectorXd(post.th.head(n)))
                .Arr("dstar", Eigen::VectorXd(post.ds.head(n)))
                .Obj("wake", Coords(foil.GetWakeCoords()));
            Write(dir, "cp.json", j);
        }

        // Polar: NACA 2412, Re 1e6, viscous and inviscid lift.
        {
            std::vector<Json> points;
            for (int a = -4; a <= 14; ++a) {
                Scirocco foil("2412", 160);
                foil.SetOper(a, -1.0, kRe, 0.0);
                const bool ok = foil.Solve();
                Scirocco inv("2412", 160);
                inv.SetOper(a, -1.0, 0.0, 0.0);
                static_cast<void>(inv.Solve());
                Json p;
                p.Num("alpha", a)
                    .Bool("converged", ok)
                    .Num("cl_inviscid", inv.GetCl());
                if (ok) {
                    p.Num("cl", foil.GetCl())
                        .Num("cd", foil.GetCd())
                        .Num("cm", foil.GetCm());
                }
                points.push_back(p);
            }
            Json j;
            j.Str("airfoil", "NACA 2412").Num("Re", kRe).Objs("points", points);
            Write(dir, "polar.json", j);
        }

        // Newton convergence with and without pseudo-transient continuation.
        {
            Json j;
            j.Str("airfoil", "NACA 2412").Num("alpha", 4.0).Num("Re", kRe);
            for (const bool ptc : {true, false}) {
                Scirocco foil("2412", 160);
                foil.GetState().param.use_ptc = ptc;
                foil.SetOper(4.0, -1.0, kRe, 0.0);
                auto [ok, hist] = SolveWithHistory(foil);
                Json h;
                h.Bool("converged", ok).Arr("residual", hist);
                j.Obj(ptc ? "ptc" : "newton", h);
            }
            // A harder start: NACA 0012 at 199 panels needs many damped
            // steps before Newton takes over.
            for (const bool ptc : {true, false}) {
                Scirocco foil("0012", 199);
                foil.GetState().param.use_ptc = ptc;
                foil.SetOper(5.0, -1.0, kRe, 0.0);
                bool ok = false;
                std::vector<double> hist;
                try {
                    std::tie(ok, hist) = SolveWithHistory(foil);
                } catch (const std::exception&) {
                    ok = false;
                }
                Json h;
                h.Bool("converged", ok).Arr("residual", hist);
                j.Obj(ptc ? "hard_ptc" : "hard_newton", h);
            }
            Write(dir, "convergence.json", j);
        }

        // Panel layout at a count small enough to see individual nodes.
        {
            Scirocco foil("2412", 60);
            Json j;
            j.Str("airfoil", "NACA 2412")
                .Obj("coords", Coords(foil.GetFoilCoords()));
            Write(dir, "panels.json", j);
        }
    } catch (const std::exception& error) {
        std::fputs("docs_figure_data failed: ", stderr);
        std::fputs(error.what(), stderr);
        std::fputs("\n", stderr);
        return 1;
    } catch (...) {
        return 1;
    }
    return 0;
}
