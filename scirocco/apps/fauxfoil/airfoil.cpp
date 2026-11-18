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
#include "airfoil.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

#include "format.h"
#include "scirocco/core/spline.h"

namespace fauxfoil {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Parse up to @p max numbers from a line; stops at the first non-number.
std::vector<double> ParseNumbers(const std::string& line, std::size_t max) {
    std::vector<double> values;
    std::string text = line;
    std::replace(text.begin(), text.end(), ',', ' ');
    std::istringstream stream(text);
    std::string token;
    while (values.size() < max && stream >> token) {
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end == token.c_str() || *end != '\0' || !std::isfinite(value))
            break;
        values.push_back(value);
    }
    return values;
}

/// Spline parameter of the leading edge (XFOIL's LEFIND): where
/// (X - Xte) . dX/ds = 0 next to the point of minimum x.
double LeadingEdgeParam(const scirocco::Spline2D& spline, const MatrixXd& x,
                        int ile) {
    const int n = static_cast<int>(x.cols());
    const scirocco::Vector2d te = 0.5 * (x.col(0) + x.col(n - 1));
    auto f = [&](double s) {
        return (spline.X.Eval(s) - te(0)) * spline.X.EvalDeriv(s) +
               (spline.Y.Eval(s) - te(1)) * spline.Y.EvalDeriv(s);
    };
    double a = spline.X.s(ile - 1);
    double b = spline.X.s(ile + 1);
    double fa = f(a);
    if (fa * f(b) > 0.0) return spline.X.s(ile);
    for (int it = 0; it < 80; ++it) {
        const double m = 0.5 * (a + b);
        const double fm = f(m);
        if (fa * fm <= 0.0) {
            b = m;
        } else {
            a = m;
            fa = fm;
        }
    }
    return 0.5 * (a + b);
}

/// Parameter in [a, b] where g(s) = target, by bisection; g is assumed
/// monotonic on the interval.
template <class G>
double Solve(const G& g, double a, double b, double target) {
    double ga = g(a) - target;
    for (int it = 0; it < 60; ++it) {
        const double m = 0.5 * (a + b);
        const double gm = g(m) - target;
        if ((ga <= 0.0) == (gm <= 0.0)) {
            a = m;
            ga = gm;
        } else {
            b = m;
        }
    }
    return 0.5 * (a + b);
}

}  // namespace

int LeadingEdgeIndex(const MatrixXd& x) {
    Eigen::Index index = 0;
    x.row(0).minCoeff(&index);
    return static_cast<int>(index);
}

double TrailingEdgeGap(const MatrixXd& x) {
    return (x.col(0) - x.col(x.cols() - 1)).norm();
}

scirocco::Vector2d LeadingEdge(const MatrixXd& x) {
    const int ile = LeadingEdgeIndex(x);
    if (ile <= 0 || ile >= x.cols() - 1) return x.col(ile);
    const scirocco::Spline2D spline = scirocco::Spline2DCreate(x);
    const double s = LeadingEdgeParam(spline, x, ile);
    return {spline.X.Eval(s), spline.Y.Eval(s)};
}

double Chord(const MatrixXd& x) {
    const scirocco::Vector2d te = 0.5 * (x.col(0) + x.col(x.cols() - 1));
    return (te - LeadingEdge(x)).norm();
}

bool IsCounterclockwise(const MatrixXd& x) {
    double area = 0.0;
    const Eigen::Index n = x.cols();
    for (Eigen::Index i = 0; i < n; ++i) {
        const Eigen::Index j = (i + 1) % n;
        area += x(0, i) * x(1, j) - x(0, j) * x(1, i);
    }
    return area > 0.0;
}

MatrixXd Reversed(const MatrixXd& x) { return x.rowwise().reverse(); }

MatrixXd Normalized(const MatrixXd& x) {
    const scirocco::Vector2d le = LeadingEdge(x);
    const double chord = Chord(x);
    MatrixXd result = x.colwise() - le;
    if (chord > 0.0) result /= chord;
    return result;
}

ThicknessCamber GetThicknessCamber(const MatrixXd& x) {
    ThicknessCamber result;
    const int ile = LeadingEdgeIndex(x);
    const int n = static_cast<int>(x.cols());
    if (ile <= 0 || ile >= n - 1) return result;

    const scirocco::Spline2D spline = scirocco::Spline2DCreate(x);
    const double s_le = LeadingEdgeParam(spline, x, ile);
    const double s_first = spline.X.s(0);
    const double s_last = spline.X.s(n - 1);

    /* Chord-line frame: leading edge at the origin, trailing edge on +x */
    const scirocco::Vector2d le(spline.X.Eval(s_le), spline.Y.Eval(s_le));
    const scirocco::Vector2d te = 0.5 * (x.col(0) + x.col(n - 1));
    scirocco::Vector2d along = te - le;
    if (along.norm() <= 0.0) return result;
    along.normalize();
    const scirocco::Vector2d normal(-along(1), along(0));
    auto point = [&](double s) -> scirocco::Vector2d {
        return scirocco::Vector2d(spline.X.Eval(s), spline.Y.Eval(s)) - le;
    };
    auto xc = [&](double s) { return point(s).dot(along); };
    auto yc = [&](double s) { return point(s).dot(normal); };

    /* Each surface runs from the leading edge to one trailing-edge point */
    const double xte = std::min(xc(s_first), xc(s_last));
    constexpr int kStations = 1000;
    for (int k = 1; k <= kStations; ++k) {
        /* Cosine spacing resolves the leading edge */
        const double xs = xte * 0.5 * (1.0 - std::cos(kPi * k / kStations));
        const double y_first = yc(Solve(xc, s_le, s_first, xs));
        const double y_second = yc(Solve(xc, s_le, s_last, xs));
        const double thickness = std::fabs(y_second - y_first);
        const double camber = 0.5 * (y_first + y_second);
        if (thickness > result.thickness) {
            result.thickness = thickness;
            result.x_thickness = xs;
        }
        if (std::fabs(camber) > std::fabs(result.camber)) {
            result.camber = camber;
            result.x_camber = xs;
        }
    }
    return result;
}

double MaxCornerAngle(const MatrixXd& x, int* node) {
    double max_angle = 0.0;
    int max_node = 0;
    for (Eigen::Index i = 1; i + 1 < x.cols(); ++i) {
        const scirocco::Vector2d a = x.col(i) - x.col(i - 1);
        const scirocco::Vector2d b = x.col(i + 1) - x.col(i);
        const double cross = a(0) * b(1) - a(1) * b(0);
        const double angle =
            std::fabs(std::atan2(cross, a.dot(b))) * 180.0 / kPi;
        if (angle > max_angle) {
            max_angle = angle;
            max_node = static_cast<int>(i);
        }
    }
    if (node) *node = max_node;
    return max_angle;
}

bool ReadAirfoilFile(const std::string& path, AirfoilFile& file,
                     std::string& error) {
    std::ifstream stream(path);
    if (!stream) {
        error = "open";
        return false;
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }

    std::size_t first = 0;
    while (first < lines.size() &&
           lines[first].find_first_not_of(" \t") == std::string::npos)
        ++first;
    if (first == lines.size()) {
        error = "empty";
        return false;
    }

    file = AirfoilFile();
    if (ParseNumbers(lines[first], 2).size() < 2) {
        file.labeled = true;
        const std::size_t begin = lines[first].find_first_not_of(" \t");
        const std::size_t end = lines[first].find_last_not_of(" \t");
        file.name = lines[first].substr(begin, end - begin + 1);
        ++first;
        /* ISES/MSES files carry a domain-box line of four numbers */
        if (first < lines.size() && ParseNumbers(lines[first], 4).size() == 4)
            ++first;
    }

    std::vector<double> xs;
    std::vector<double> ys;
    for (std::size_t i = first; i < lines.size(); ++i) {
        const std::vector<double> values = ParseNumbers(lines[i], 2);
        if (values.size() < 2) {
            if (lines[i].find_first_not_of(" \t") == std::string::npos)
                continue;
            break;
        }
        /* Multi-element separators (999 999) end the first element */
        if (values[0] == 999.0 && values[1] == 999.0) break;
        xs.push_back(values[0]);
        ys.push_back(values[1]);
    }
    if (xs.size() < 3) {
        error = "points";
        return false;
    }
    file.x.resize(2, static_cast<Eigen::Index>(xs.size()));
    for (std::size_t i = 0; i < xs.size(); ++i) {
        file.x(0, static_cast<Eigen::Index>(i)) = xs[i];
        file.x(1, static_cast<Eigen::Index>(i)) = ys[i];
    }
    return true;
}

bool WriteAirfoilFile(const std::string& path, const std::string* name,
                      const MatrixXd& x) {
    std::ofstream stream(path);
    if (!stream) return false;
    if (name) stream << *name << '\n';
    for (Eigen::Index i = 0; i < x.cols(); ++i)
        stream << ' ' << G(x(0, i), 15, 7) << G(x(1, i), 15, 7) << '\n';
    return static_cast<bool>(stream);
}

}  // namespace fauxfoil
