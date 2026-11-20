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
#include "plots.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

/* Layout constants are page positions measured on XFOIL 6.99's default
   window (1245 x 962), counted from the top edge as on screen; Top()
   converts them to page coordinates. */

namespace fauxfoil {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Page y of a position measured down from the top of the window.
double Top(double y) { return Scene::kHeight - y; }

/// Decimals needed to print multiples of @p del exactly.
int DecimalsFor(double del) {
    del = std::fabs(del);
    if (del <= 0.0) return 1;
    for (int d = 0; d < 7; ++d) {
        const double scaled = del * std::pow(10.0, d);
        if (std::fabs(scaled - std::round(scaled)) < 1e-6 * scaled) return d;
    }
    return 7;
}

/// Label spacing that gives at most @p intervals intervals over [lo, hi].
Axis NiceAxis(double lo, double hi, int intervals = 8) {
    if (!std::isfinite(lo) || !std::isfinite(hi)) return {0.0, 1.0, 0.2};
    if (hi - lo < 1e-12) {
        const double pad = std::max(std::fabs(hi) * 0.1, 1e-3);
        lo -= pad;
        hi += pad;
    }
    const double raw = (hi - lo) / intervals;
    const double decade = std::pow(10.0, std::floor(std::log10(raw)) - 1.0);
    for (int k = 0; k < 4; ++k) {
        for (double m : {1.0, 2.0, 4.0, 5.0}) {
            const double step = m * decade * std::pow(10.0, k);
            const double a = std::floor(lo / step + 1e-9) * step;
            const double b = std::ceil(hi / step - 1e-9) * step;
            if ((b - a) / step <= intervals + 1e-9) return {a, b, step};
        }
    }
    return {lo, hi, (hi - lo) / intervals};
}

/// Values a, a+del, ... up to b (inclusive, with tolerance).
std::vector<double> Steps(double a, double b, double del) {
    std::vector<double> values;
    if (del <= 0.0) return values;
    const int n = static_cast<int>(std::floor((b - a) / del + 1e-6));
    for (int k = 0; k <= n; ++k) values.push_back(a + k * del);
    return values;
}

/// "1.000×10^6" as XFOIL writes Reynolds numbers on plots.
std::string TimesTen(double value) {
    if (value <= 0.0) return "0";
    const int exponent = static_cast<int>(std::floor(std::log10(value)));
    const double mantissa = value / std::pow(10.0, exponent);
    return Fixed(mantissa, 3) + "×10^{" + std::to_string(exponent) + "}";
}

/// Text height of a label centred on a tick: baseline half a height down.
void TickLabel(Scene& scene, double x, double y, double height,
               const std::string& text, double align) {
    scene.Text(x, y - 0.5 * height, height, text, Ink::Fore, align);
}

/* ---------------------------------------------------------------- */
/* Cp(x)                                                             */
/* ---------------------------------------------------------------- */

/// Leading edge and chord of the paneled airfoil, for scaling.
struct Chord {
    double xle = 0.0;
    double yle = 0.0;
    double chord = 1.0;
};

Chord ChordOf(const std::vector<Point>& airfoil) {
    Chord c;
    if (airfoil.size() < 3) return c;
    std::size_t le = 0;
    for (std::size_t i = 1; i < airfoil.size(); ++i)
        if (airfoil[i].x < airfoil[le].x) le = i;
    c.xle = airfoil[le].x;
    c.yle = airfoil[le].y;
    const Point& a = airfoil.front();
    const Point& b = airfoil.back();
    const double xte = 0.5 * (a.x + b.x);
    const double yte = 0.5 * (a.y + b.y);
    c.chord = std::hypot(xte - c.xle, yte - c.yle);
    if (c.chord <= 0.0) c.chord = 1.0;
    return c;
}

/// Page mapping of the Cp(x) plot. The constants are for an airfoil whose
/// thickness puts the Cp = 1 tick at 803; the whole Cp part moves up or
/// down by @c shift so the airfoil always fits below it.
struct CpFrame {
    Chord c;
    double cp_min = -2.0;
    double cp_max = 1.0;
    double shift = 0.0;        ///< Downward offset of the Cp part
    double airfoil_y = 880.0;  ///< Chord line, from the top
    static constexpr double kAxisX = 105.0;
    static constexpr double kAxisTop = 183.0;
    static constexpr double kAxisBottom = 803.0;
    static constexpr double kLeadingEdgeX = 195.0;
    static constexpr double kChord = 750.0;
    static constexpr double kAirfoilBottom = 912.0;
    static constexpr double kAirfoilGap = 20.0;

    [[nodiscard]] double X(double x) const {
        return kLeadingEdgeX + kChord * (x - c.xle) / c.chord;
    }
    [[nodiscard]] double Ycp(double cp) const {
        return Top(shift + kAxisTop +
                   (cp - cp_min) * (kAxisBottom - kAxisTop) /
                       (cp_max - cp_min));
    }
    [[nodiscard]] Point Foil(const Point& p) const {
        return {X(p.x), Top(airfoil_y) + kChord * (p.y - c.yle) / c.chord};
    }
};

void CpAxes(Scene& scene, const CpFrame& f, const CpOptions& options) {
    constexpr double kNumber = 17.0;
    scene.Segment(CpFrame::kAxisX, f.Ycp(f.cp_min), CpFrame::kAxisX,
                  f.Ycp(f.cp_max));
    for (double cp :
         Steps(std::ceil(f.cp_min / 0.5 - 1e-9) * 0.5, f.cp_max, 0.5)) {
        const double y = f.Ycp(cp);
        scene.Segment(CpFrame::kAxisX - 6.0, y, CpFrame::kAxisX + 2.0, y);
        TickLabel(scene, 91.0, y, kNumber, Fixed(cp, 1), 1.0);
    }
    scene.Text(33.0, Top(350.0 + f.shift), 25.0, "C_{p}");

    /* x/c axis along Cp = 0 */
    const double y0 = f.Ycp(0.0);
    scene.Segment(CpFrame::kAxisX, y0, 1230.0, y0);
    for (int k = 1; k <= 10; ++k) {
        const double x = f.X(f.c.xle + 0.1 * k * f.c.chord);
        const double half = k == 10 ? 7.0 : 4.0;
        scene.Segment(x, y0 - half, x, y0 + half);
    }

    if (options.grid) {
        for (int k = 0; k <= 10; ++k) {
            const double x = f.X(f.c.xle + 0.1 * k * f.c.chord);
            scene.Segment(x, f.Ycp(f.cp_min), x, f.Ycp(f.cp_max), Ink::Fore,
                          Dash::Dotted);
        }
        for (double cp :
             Steps(std::ceil(f.cp_min / 0.1 - 1e-9) * 0.1, f.cp_max, 0.1))
            scene.Segment(CpFrame::kAxisX, f.Ycp(cp), f.X(f.c.xle + f.c.chord),
                          f.Ycp(cp), Ink::Fore, Dash::Dotted);
    }

    scene.Text(112.0, Top(193.0 + f.shift), 11.8, "fauXfoil");
    scene.Text(112.0, Top(214.0 + f.shift), 9.4, "V " + options.version);
}

std::vector<Point> CpCurve(const CpFrame& f, const std::vector<Point>& at,
                           const std::vector<double>& cp,
                           const std::vector<int>* nodes = nullptr) {
    std::vector<Point> curve;
    if (nodes != nullptr) {
        for (int i : *nodes)
            if (i >= 0 && static_cast<std::size_t>(i) < cp.size() &&
                static_cast<std::size_t>(i) < at.size())
                curve.push_back({f.X(at[i].x), f.Ycp(cp[i])});
        return curve;
    }
    for (std::size_t i = 0; i < at.size() && i < cp.size(); ++i)
        curve.push_back({f.X(at[i].x), f.Ycp(cp[i])});
    return curve;
}

std::vector<Point> FoilCurve(const CpFrame& f, const std::vector<Point>& at) {
    std::vector<Point> curve;
    curve.reserve(at.size());
    for (const Point& p : at) curve.push_back(f.Foil(p));
    return curve;
}

CpFrame FrameFor(const Solution& solution, const CpOptions& options) {
    CpFrame f;
    f.c = ChordOf(solution.airfoil);
    f.cp_min = std::min(options.cp_min, -0.5);
    if (!solution.airfoil.empty()) {
        double lo = 0.0;
        double hi = 0.0;
        for (const Point& p : solution.airfoil) {
            lo = std::min(lo, (p.y - f.c.yle) / f.c.chord);
            hi = std::max(hi, (p.y - f.c.yle) / f.c.chord);
        }
        f.airfoil_y = CpFrame::kAirfoilBottom + CpFrame::kChord * lo;
        const double cp_bottom = CpFrame::kAirfoilBottom -
                                 CpFrame::kChord * (hi - lo) -
                                 CpFrame::kAirfoilGap;
        f.shift = cp_bottom - CpFrame::kAxisBottom;
    }
    return f;
}

/// One legend line: label, "=" and value in fixed columns.
struct LegendColumns {
    double label;
    double equals;
    double value;
};

void LegendLine(Scene& scene, const LegendColumns& c, double baseline,
                double height, const std::string& label,
                const std::string& value, Ink ink = Ink::Fore) {
    scene.Text(c.label, baseline, height, label, ink);
    scene.Text(c.equals, baseline, height, "=", ink);
    scene.Text(c.value, baseline, height, value, ink);
}

}  // namespace

std::string Fixed(double value, int decimals) {
    char buffer[64];
    std::snprintf(buffer, sizeof buffer, "%.*f", decimals, value);
    std::string text = buffer;
    if (text[0] == '-' && text.find_first_not_of("-0.") == std::string::npos)
        text.erase(0, 1);
    return text;
}

Scene CpPlot(const Solution& s, const CpOptions& options) {
    Scene scene;
    const CpFrame f = FrameFor(s, options);
    CpAxes(scene, f, options);

    /* Airfoil, with displacement surfaces when viscous */
    scene.Line(FoilCurve(f, s.airfoil));
    if (s.viscous) {
        scene.Line(FoilCurve(f, s.dstar_top), Ink::Blue);
        scene.Line(FoilCurve(f, s.dstar_bottom), Ink::Red);
        scene.Line(FoilCurve(f, s.wake_upper));
        scene.Line(FoilCurve(f, s.wake_lower));
    }

    /* Cp */
    if (s.viscous) {
        if (options.inviscid) {
            scene.Line(CpCurve(f, s.airfoil, s.cpi), Ink::Fore, Dash::Dashed);
            scene.Line(CpCurve(f, s.wake, s.cpi_wake), Ink::Fore, Dash::Dashed);
        }
        scene.Line(CpCurve(f, s.airfoil, s.cp, &s.upper), Ink::Blue);
        scene.Line(CpCurve(f, s.airfoil, s.cp, &s.lower), Ink::Red);
        scene.Line(CpCurve(f, s.wake, s.cp_wake));
    } else {
        scene.Line(CpCurve(f, s.airfoil, s.cp));
    }

    /* Operating point */
    scene.Text(898.0, Top(202.0 + f.shift), 21.6, s.name);
    const LegendColumns columns{892.0, 938.0, 968.0};
    constexpr double kHeight = 17.5;
    double baseline = Top(240.0 + f.shift);
    auto line = [&](const std::string& label, const std::string& value) {
        LegendLine(scene, columns, baseline, kHeight, label, value);
        baseline -= 33.5;
    };
    if (s.mach > 0.0) line("Ma", Fixed(s.mach, 3));
    if (s.viscous) line("Re", TimesTen(s.re));
    line("α", Fixed(s.alpha, 4) + "°");
    line("C_{L}", Fixed(s.cl, 4));
    line("C_{M}", Fixed(s.cm, 4));
    if (s.viscous) {
        line("C_{D}", Fixed(s.cd, 5));
        /* XFOIL letters "L/D" smaller so that it fits the label column */
        scene.Text(columns.label, baseline, 0.72 * kHeight, "L/D");
        scene.Text(columns.equals, baseline, kHeight, "=");
        scene.Text(columns.value, baseline, kHeight,
                   s.cd != 0.0 ? Fixed(s.cl / s.cd, 2) : "");
        baseline -= 33.5;
        line("N_{cr}", Fixed(s.ncrit, 2));
    } else {
        line("C_{Dp}", Fixed(s.cdp, 5));
    }
    return scene;
}

Scene CpSequencePlot(const std::vector<Solution>& solutions,
                     const CpOptions& options) {
    Scene scene;
    if (solutions.empty()) return scene;
    const Solution& last = solutions.back();
    const CpFrame f = FrameFor(last, options);
    CpAxes(scene, f, options);
    scene.Line(FoilCurve(f, last.airfoil));
    for (const Solution& s : solutions) {
        scene.Line(CpCurve(f, s.airfoil, s.cp));
        if (s.viscous) scene.Line(CpCurve(f, s.wake, s.cp_wake));
    }

    /* Fixed parameters */
    scene.Text(870.0, Top(202.0 + f.shift), 21.6, last.name, Ink::Fore, 0.5);
    const LegendColumns columns{780.0, 827.0, 857.0};
    double baseline = Top(235.0 + f.shift);
    auto line = [&](const std::string& label, const std::string& value) {
        LegendLine(scene, columns, baseline, 17.5, label, value);
        baseline -= 31.5;
    };
    line("Ma", Fixed(last.mach, 3));
    if (last.viscous) {
        line("Re", TimesTen(last.re));
        line("N_{cr}", Fixed(last.ncrit, 3));
    }

    /* Table of results */
    struct Column {
        double x;
        int chars;
        const char* head;
    };
    std::vector<Column> table = {
        {645.0, 6, "α"}, {727.0, 6, "C_{L}"}, {822.0, 6, "C_{M}"}};
    if (last.viscous) {
        table.push_back({917.0, 7, "C_{D}"});
        table.push_back({1022.0, 5, "Top Xtr"});
        table.push_back({1105.0, 5, "Bot Xtr"});
    }
    constexpr double kHead = 13.0;
    constexpr double kRow = 14.5;
    const double head_y = Top(330.0 + f.shift);
    for (const Column& c : table) {
        const double width = c.chars * kAdvance * kRow;
        const bool words = std::string(c.head).find(' ') != std::string::npos;
        const double h = words ? 0.8 * kHead : kHead;
        scene.Text(c.x + 0.5 * width, head_y, h, c.head, Ink::Fore, 0.5);
        scene.Segment(c.x - 2.0, Top(337.0 + f.shift), c.x + width + 2.0,
                      Top(337.0 + f.shift));
    }
    const std::size_t first = solutions.size() > 16 ? solutions.size() - 16 : 0;
    double row = Top(362.0 + f.shift);
    for (std::size_t k = first; k < solutions.size(); ++k) {
        const Solution& s = solutions[k];
        const std::string values[6] = {
            Fixed(s.alpha, 3), Fixed(s.cl, 4),      Fixed(s.cm, 3),
            Fixed(s.cd, 5),    Fixed(s.xtr_top, 3), Fixed(s.xtr_bot, 3)};
        for (std::size_t c = 0; c < table.size(); ++c)
            scene.Text(table[c].x, row, kRow, values[c]);
        row -= 25.0;
    }
    return scene;
}

Scene CpVectorPlot(const Solution& s) {
    Scene scene;
    if (s.airfoil.size() < 3) return scene;
    const Chord c = ChordOf(s.airfoil);
    constexpr double kLeX = 262.0;
    constexpr double kLeY = 713.0;
    constexpr double kScale = 886.0;
    const double a = s.alpha * kPi / 180.0;
    const double ca = std::cos(a);
    const double sa = std::sin(a);
    std::vector<Point> page;
    for (const Point& p : s.airfoil) {
        const double x = (p.x - c.xle) / c.chord;
        const double y = (p.y - c.yle) / c.chord;
        page.push_back({kLeX + kScale * (ca * x + sa * y),
                        Top(kLeY) + kScale * (-sa * x + ca * y)});
    }
    scene.Line(page);

    /* Pressure vectors: suction points away from the surface, pressure
       into it; length 0.3 chord per unit Cp */
    const std::size_t n = page.size();
    for (std::size_t i = 0; i < n && i < s.cp.size(); ++i) {
        const Point& prev = page[i == 0 ? 0 : i - 1];
        const Point& next = page[i + 1 < n ? i + 1 : n - 1];
        double tx = next.x - prev.x;
        double ty = next.y - prev.y;
        const double t = std::hypot(tx, ty);
        if (t <= 0.0) continue;
        tx /= t;
        ty /= t;
        const double nx = -ty;
        const double ny = tx;
        const double length = 0.3 * kScale * std::fabs(s.cp[i]);
        if (length < 1.0) continue;
        const Point surface = page[i];
        const Point outer{surface.x + nx * length, surface.y + ny * length};
        const bool suction = s.cp[i] < 0.0;
        const Point tip = suction ? outer : surface;
        const Point tail = suction ? surface : outer;
        scene.Segment(tail.x, tail.y, tip.x, tip.y);
        const double dx = (tip.x - tail.x) / length;
        const double dy = (tip.y - tail.y) / length;
        const double head = std::min(0.1 * length + 3.0, 30.0);
        const double half = 0.22 * head;
        const Point base{tip.x - dx * head, tip.y - dy * head};
        scene.Line({tip,
                    {base.x - dy * half, base.y + dx * half},
                    {base.x + dy * half, base.y - dx * half},
                    tip});
    }
    return scene;
}

/* ---------------------------------------------------------------- */
/* Boundary-layer variables                                          */
/* ---------------------------------------------------------------- */

namespace {

struct BlCurve {
    std::vector<Point> points;  ///< (x or s, value)
    Ink ink = Ink::Fore;
    Dash dash = Dash::Solid;
};

/// Value of @p variable at a station; false where it is not plotted.
bool BlValue(const BlStation& st, BlVariable variable, bool second,
             double& value) {
    switch (variable) {
        case BlVariable::Hk:
            value = st.hk;
            return true;
        case BlVariable::DstarThetaTop:
        case BlVariable::DstarThetaBottom:
            value = second ? st.th : st.ds;
            return true;
        case BlVariable::Ue:
            value = st.ue;
            return true;
        case BlVariable::Cf:
            value = st.cf;
            return true;
        case BlVariable::Dissipation:
            value = st.cd;
            return true;
        case BlVariable::Amplification:
            value = st.sa;
            return !st.turbulent;
        case BlVariable::Ctau:
            value = second ? st.cteq : st.sa;
            return st.turbulent;
        case BlVariable::RetTheta:
            value = st.ret;
            return true;
        case BlVariable::LogRetTheta:
            if (st.ret <= 0.0) return false;
            value = std::log10(st.ret);
            return true;
    }
    return false;
}

std::vector<BlCurve> BlCurves(const Solution& s, BlVariable variable,
                              bool vs_s) {
    struct Part {
        std::vector<const std::vector<BlStation>*> stations;
        Ink ink;
    };
    std::vector<Part> parts;
    if (variable == BlVariable::DstarThetaTop)
        parts.push_back({{&s.bl_top, &s.bl_wake}, Ink::Blue});
    else if (variable == BlVariable::DstarThetaBottom)
        parts.push_back({{&s.bl_bottom, &s.bl_wake}, Ink::Red});
    else if (variable == BlVariable::Amplification)
        parts = {{{&s.bl_top}, Ink::Blue}, {{&s.bl_bottom}, Ink::Red}};
    else
        parts = {{{&s.bl_top}, Ink::Blue},
                 {{&s.bl_bottom, &s.bl_wake}, Ink::Red}};

    const bool two = variable == BlVariable::DstarThetaTop ||
                     variable == BlVariable::DstarThetaBottom ||
                     variable == BlVariable::Ctau;
    std::vector<BlCurve> curves;
    for (const Part& part : parts) {
        for (int k = 0; k < (two ? 2 : 1); ++k) {
            BlCurve curve;
            curve.ink = part.ink;
            if (variable == BlVariable::Ctau && k == 1)
                curve.dash = Dash::Dashed;
            for (const auto* list : part.stations) {
                for (const BlStation& st : *list) {
                    double value = 0.0;
                    if (!BlValue(st, variable, k == 1, value)) {
                        if (curve.points.size() > 1) curves.push_back(curve);
                        curve.points.clear();
                        continue;
                    }
                    curve.points.push_back({vs_s ? st.s : st.x, value});
                }
            }
            if (curve.points.size() > 1) curves.push_back(curve);
        }
    }
    return curves;
}

}  // namespace

BlAxes DefaultBlAxes(const Solution& s, BlVariable variable,
                     const BlAxes& current) {
    BlAxes axes = current;
    if (!axes.x_set) {
        if (axes.vs_s) {
            double smax = 0.0;
            for (const auto* list : {&s.bl_top, &s.bl_bottom})
                if (!list->empty()) smax = std::max(smax, list->back().s);
            axes.x = NiceAxis(0.0, 1.4 * (smax > 0.0 ? smax : 1.0), 7);
        } else {
            axes.x = {0.0, 1.4, 0.2};
        }
    }
    if (axes.y_set) return axes;

    double lo = std::numeric_limits<double>::infinity();
    double hi = -lo;
    for (const BlCurve& curve : BlCurves(s, variable, axes.vs_s))
        for (const Point& p : curve.points)
            if (p.x >= axes.x.min && p.x <= axes.x.max && std::isfinite(p.y)) {
                lo = std::min(lo, p.y);
                hi = std::max(hi, p.y);
            }
    if (!std::isfinite(lo)) {
        lo = 0.0;
        hi = 1.0;
    }
    switch (variable) {
        case BlVariable::Hk:
            /* XFOIL's fixed 0..6; separated flow widens it */
            axes.y = hi <= 7.0 ? Axis{0.0, 6.0, 1.0}
                               : NiceAxis(0.0, std::min(hi, 20.0), 8);
            break;
        case BlVariable::Amplification: {
            const double top = std::max(hi, s.ncrit);
            const double del = top > 8.0 ? 4.0 : 2.0;
            axes.y = {0.0, del * std::ceil((top + 0.5) / del), del};
            break;
        }
        case BlVariable::LogRetTheta:
            axes.y = NiceAxis(std::floor(lo), std::ceil(hi), 6);
            break;
        case BlVariable::Cf:
            /* Small negative values (wake, separation onset) do not widen
               the axis below zero */
            axes.y = NiceAxis(lo < -0.1 * hi ? lo : 0.0, hi, 5);
            break;
        case BlVariable::DstarThetaTop:
        case BlVariable::DstarThetaBottom:
            axes.y = NiceAxis(0.0, hi, 5);
            break;
        default:
            axes.y = NiceAxis(std::min(0.0, lo), hi, 8);
            break;
    }
    return axes;
}

Scene BlPlot(const Solution& s, BlVariable variable, const BlAxes& in) {
    Scene scene;
    const BlAxes axes = DefaultBlAxes(s, variable, in);
    constexpr double kLeft = 137.0;
    constexpr double kRight = 1153.0;
    const double bottom = Top(860.0);
    const double top = Top(240.0);
    auto X = [&](double x) {
        return kLeft +
               (x - axes.x.min) * (kRight - kLeft) / (axes.x.max - axes.x.min);
    };
    auto Y = [&](double y) {
        return bottom +
               (y - axes.y.min) * (top - bottom) / (axes.y.max - axes.y.min);
    };

    /* Header */
    scene.Text(135.0, Top(125.0), 23.4, s.name);
    constexpr double kHeader = 19.5;
    const double row1 = Top(172.0);
    const double row2 = Top(214.0);
    LegendLine(scene, {135.0, 190.0, 220.0}, row1, kHeader, "Ma",
               Fixed(s.mach, 4));
    LegendLine(scene, {135.0, 190.0, 220.0}, row2, kHeader, "Re",
               TimesTen(s.re));
    LegendLine(scene, {405.0, 460.0, 490.0}, row1, kHeader, "α",
               Fixed(s.alpha, 4) + "°");
    LegendLine(scene, {405.0, 460.0, 490.0}, row2, kHeader, "N_{cr}",
               Fixed(s.ncrit, 2));
    LegendLine(scene, {643.0, 700.0, 728.0}, row1, kHeader, "C_{L}",
               Fixed(s.cl, 4));
    LegendLine(scene, {643.0, 700.0, 728.0}, row2, kHeader, "C_{D}",
               Fixed(s.cd, 5));
    LegendLine(scene, {897.0, 1020.0, 1050.0}, row1, kHeader, "T:x_{tr}/c",
               Fixed(s.xtr_top, 4), Ink::Blue);
    LegendLine(scene, {897.0, 1020.0, 1050.0}, row2, kHeader, "B:x_{tr}/c",
               Fixed(s.xtr_bot, 4), Ink::Red);

    /* Axes, numbers and grid */
    constexpr double kNumber = 18.5;
    scene.Segment(kLeft, bottom, kRight, bottom);
    scene.Segment(kLeft, bottom, kLeft, top);
    const int xd = std::max(1, DecimalsFor(axes.x.del));
    for (double x : Steps(axes.x.min, axes.x.max, axes.x.del)) {
        scene.Segment(X(x), bottom, X(x), bottom - 6.0);
        scene.Text(X(x), Top(895.0), kNumber, Fixed(x, xd), Ink::Fore, 0.5);
    }
    const int yd = std::max(1, DecimalsFor(axes.y.del));
    for (double y : Steps(axes.y.min, axes.y.max, axes.y.del)) {
        scene.Segment(kLeft - 6.0, Y(y), kLeft, Y(y));
        TickLabel(scene, 122.0, Y(y), kNumber, Fixed(y, yd), 1.0);
    }
    if (axes.grid) {
        for (double x : Steps(axes.x.min, axes.x.max, 0.5 * axes.x.del))
            if (x > axes.x.min + 1e-9)
                scene.Segment(X(x), bottom, X(x), top, Ink::Fore, Dash::Dotted);
        for (double y : Steps(axes.y.min, axes.y.max, 0.5 * axes.y.del))
            if (y > axes.y.min + 1e-12)
                scene.Segment(kLeft, Y(y), kRight, Y(y), Ink::Fore,
                              Dash::Dotted);
    }
    scene.Text(930.0, Top(918.0), 20.0, axes.vs_s ? "s" : "X");

    /* Axis title */
    switch (variable) {
        case BlVariable::Hk:
            scene.Text(68.0, Top(400.0), 22.0, "H_{k}");
            break;
        case BlVariable::DstarThetaTop:
        case BlVariable::DstarThetaBottom: {
            const bool topside = variable == BlVariable::DstarThetaTop;
            const Ink ink = topside ? Ink::Blue : Ink::Red;
            scene.Text(topside ? 58.0 : 30.0, Top(305.0), 20.0,
                       topside ? "Top" : "Bottom", ink);
            scene.Text(68.0, Top(427.0), 20.0, "δ*", ink);
            scene.Text(68.0, Top(552.0), 20.0, "θ", ink);
            break;
        }
        case BlVariable::Ue:
            scene.Text(44.0, Top(356.0), 17.0, "U_{e}/V_{∞}");
            break;
        case BlVariable::Cf:
            scene.Text(78.0, Top(426.0), 22.0, "C_{f}");
            break;
        case BlVariable::Dissipation:
            scene.Text(70.0, Top(426.0), 22.0, "C_{D}");
            break;
        case BlVariable::Amplification:
            scene.Text(72.0, Top(346.0), 20.0, "ln");
            scene.Text(60.0, Top(550.0), 20.0, "A/A_{0}");
            scene.Segment(kLeft, Y(s.ncrit), kRight, Y(s.ncrit), Ink::Fore,
                          Dash::Dotted);
            break;
        case BlVariable::Ctau:
            scene.Text(62.0, Top(400.0), 20.0, "√C_{τ}");
            break;
        case BlVariable::RetTheta:
            scene.Text(62.0, Top(400.0), 20.0, "Re_{θ}");
            break;
        case BlVariable::LogRetTheta:
            scene.Text(22.0, Top(400.0), 18.0, "log Re_{θ}");
            break;
    }

    /* Curves */
    for (const BlCurve& curve : BlCurves(s, variable, axes.vs_s)) {
        std::vector<Point> page;
        page.reserve(curve.points.size());
        for (const Point& p : curve.points) page.push_back({X(p.x), Y(p.y)});
        scene.Line(page, curve.ink, curve.dash);
    }
    return scene;
}

/* ---------------------------------------------------------------- */
/* Polars                                                            */
/* ---------------------------------------------------------------- */

namespace {

/// Colours of successive polars.
Ink PolarInk(std::size_t k) {
    constexpr Ink kInks[] = {Ink::Red,    Ink::Blue,    Ink::Green,
                             Ink::Orange, Ink::Magenta, Ink::Cyan};
    return kInks[k % (sizeof kInks / sizeof kInks[0])];
}

/// Expand [lo, hi] in steps of @p del to cover the data.
void Cover(double& lo, double& hi, double del, double dmin, double dmax) {
    if (!std::isfinite(dmin)) return;
    lo = std::min(lo, std::floor(dmin / del + 1e-9) * del);
    hi = std::max(hi, std::ceil(dmax / del - 1e-9) * del);
}

struct Box {
    double left, right, bottom, top;  ///< Page coordinates
};

void Frame(Scene& scene, const Box& b, Dash dash = Dash::Solid) {
    scene.Line({{b.left, b.bottom},
                {b.right, b.bottom},
                {b.right, b.top},
                {b.left, b.top},
                {b.left, b.bottom}},
               Ink::Fore, dash);
}

}  // namespace

Scene PolarPlot(const std::vector<Polar>& polars, const std::string& version) {
    Scene scene;

    /* Axis limits: XFOIL's defaults, widened to cover the points */
    double cl_lo = 0.0, cl_hi = 1.5;
    double cd_lo = 0.0, cd_hi = 0.04;
    double a_lo = -4.0, a_hi = 10.0;
    double cm_lo = -0.25;
    {
        const double inf = std::numeric_limits<double>::infinity();
        double clmin = inf, clmax = -inf, cdmax = -inf, amin = inf, amax = -inf,
               cmmin = inf;
        for (const Polar& p : polars)
            for (const PolarPoint& q : p.points) {
                clmin = std::min(clmin, q.cl);
                clmax = std::max(clmax, q.cl);
                cdmax = std::max(cdmax, q.cd);
                amin = std::min(amin, q.alpha);
                amax = std::max(amax, q.alpha);
                cmmin = std::min(cmmin, q.cm);
            }
        Cover(cl_lo, cl_hi, 0.5, clmin, clmax);
        Cover(cd_lo, cd_hi, 0.01, 0.0, cdmax);
        Cover(a_lo, a_hi, 2.0, amin, amax);
        if (std::isfinite(cmmin))
            cm_lo = std::min(cm_lo, std::floor(cmmin / 0.05 - 1e-9) * 0.05);
    }

    const double bottom = Top(843.0);
    const double top = Top(285.0);
    auto Ycl = [&](double cl) {
        return bottom + (cl - cl_lo) * (top - bottom) / (cl_hi - cl_lo);
    };
    /* CM shares the CL scale at 4:1, negative upward, zero at CL = 0 */
    auto Ycm = [&](double cm) { return Ycl(-4.0 * cm); };

    constexpr double kNumber = 12.0;
    constexpr double kTitle = 22.0;
    const double numbers_y = Top(868.0);

    /* Panel 1: CL vs CD */
    const Box b1{86.0, 543.0, bottom, top};
    auto Xcd = [&](double cd) {
        return b1.left + (cd - cd_lo) * (b1.right - b1.left) / (cd_hi - cd_lo);
    };
    for (double cd : Steps(cd_lo, cd_hi, 0.002))
        scene.Segment(Xcd(cd), bottom, Xcd(cd), top, Ink::Fore, Dash::Dotted);
    for (double cl : Steps(cl_lo, cl_hi, 0.1))
        scene.Segment(b1.left, Ycl(cl), b1.right, Ycl(cl), Ink::Fore,
                      Dash::Dotted);
    for (double cd : Steps(cd_lo, cd_hi, 0.01)) {
        scene.Segment(Xcd(cd), bottom, Xcd(cd), top);
        if (std::fabs(cd) > 1e-12)
            scene.Text(Xcd(cd), numbers_y, kNumber, Fixed(1e4 * cd, 0),
                       Ink::Fore, 0.5);
    }
    for (double cl : Steps(cl_lo, cl_hi, 0.5)) {
        scene.Segment(b1.left, Ycl(cl), b1.right, Ycl(cl));
        TickLabel(scene, 77.0, Ycl(cl), kNumber, Fixed(cl, 1), 1.0);
    }
    Frame(scene, b1);
    scene.Text(0.5 * (b1.left + b1.right), Top(898.0), kTitle, "10^{4} × C_{D}",
               Ink::Fore, 0.5);
    scene.Text(40.0, Top(388.0), kTitle, "C_{L}");

    /* Panel 2: CL and CM vs alpha */
    const Box b2{594.0, 848.0, bottom, top};
    auto Xa = [&](double a) {
        return b2.left + (a - a_lo) * (b2.right - b2.left) / (a_hi - a_lo);
    };
    for (double a : Steps(a_lo, a_hi, 1.0))
        scene.Segment(Xa(a), bottom, Xa(a), top, Ink::Fore, Dash::Dotted);
    for (double cl : Steps(cl_lo, cl_hi, 0.1))
        scene.Segment(b2.left, Ycl(cl), b2.right, Ycl(cl), Ink::Fore,
                      Dash::Dotted);
    const double axis_a = std::clamp(0.0, a_lo, a_hi);
    scene.Segment(Xa(axis_a), bottom, Xa(axis_a), top);
    for (double a : Steps(a_lo, a_hi, 2.0))
        if (std::fabs(a - axis_a) > 1e-9)
            scene.Text(Xa(a), numbers_y, kNumber, Fixed(a, 0), Ink::Fore, 0.5);
    for (double cl : Steps(cl_lo, cl_hi, 0.5))
        if (cl > 0.0)
            TickLabel(scene, Xa(axis_a) + 11.0, Ycl(cl), kNumber, Fixed(cl, 1),
                      0.0);
    for (int k = 1; -0.05 * k >= cm_lo - 1e-9; ++k) {
        const double cm = -0.05 * k;
        if (Ycm(cm) <= top + 1e-6)
            TickLabel(scene, Xa(axis_a) - 10.0, Ycm(cm), kNumber, Fixed(cm, 2),
                      1.0);
    }
    scene.Text(Xa(axis_a) + 26.0, Top(388.0), kTitle, "C_{L}");
    scene.Text(Xa(axis_a) - 47.0, Top(518.0), kTitle, "C_{M}");
    scene.Text(0.5 * (b2.left + b2.right) + 72.0, Top(900.0), kTitle, "α",
               Ink::Fore, 0.5);

    /* Panel 3: CL vs transition */
    const Box b3{898.0, 1102.0, bottom, top};
    auto Xt = [&](double x) { return b3.left + x * (b3.right - b3.left); };
    for (double x : Steps(0.0, 1.0, 0.1))
        scene.Segment(Xt(x), bottom, Xt(x), top, Ink::Fore, Dash::Dotted);
    for (double cl : Steps(cl_lo, cl_hi, 0.1))
        scene.Segment(b3.left, Ycl(cl), b3.right, Ycl(cl), Ink::Fore,
                      Dash::Dotted);
    scene.Segment(Xt(0.5), bottom, Xt(0.5), top);
    for (double cl : Steps(cl_lo, cl_hi, 0.5))
        scene.Segment(b3.left, Ycl(cl), b3.right, Ycl(cl));
    for (double x : {0.0, 0.5, 1.0})
        scene.Text(Xt(x), numbers_y, kNumber, Fixed(x, 1), Ink::Fore, 0.5);
    Frame(scene, b3);
    scene.Text(1025.0, Top(898.0), 17.0, "x_{tr} /c");
    scene.Text(b3.right, Top(281.0), 9.9, "fauXfoil V " + version, Ink::Fore,
               1.0);

    /* Polars */
    for (std::size_t k = 0; k < polars.size(); ++k) {
        const Polar& polar = polars[k];
        const Ink ink = PolarInk(k);
        std::string title = polar.name;
        if (polar.re > 0.0) {
            char re[32];
            std::snprintf(re, sizeof re, "%.0f", polar.re);
            title += std::string("   Re = ") + re;
        }
        title += "   Ma = " + Fixed(polar.mach, 3);
        if (polar.re > 0.0) title += "   Ncrit = " + Fixed(polar.ncrit_top, 3);
        scene.Text(
            85.0,
            Top(262.0) + 20.0 * static_cast<double>(polars.size() - 1 - k),
            14.5, title, ink);

        std::vector<PolarPoint> points = polar.points;
        std::stable_sort(points.begin(), points.end(),
                         [](const PolarPoint& a, const PolarPoint& b) {
                             return a.alpha < b.alpha;
                         });
        std::vector<Point> drag, lift, moment, top_tr, bot_tr;
        for (const PolarPoint& q : points) {
            const double y = Ycl(q.cl);
            if (polar.re > 0.0) drag.push_back({Xcd(q.cd), y});
            lift.push_back({Xa(q.alpha), y});
            moment.push_back({Xa(q.alpha), Ycm(q.cm)});
            if (polar.re > 0.0) {
                top_tr.push_back({Xt(q.xtr_top), y});
                bot_tr.push_back({Xt(q.xtr_bot), y});
            }
        }
        scene.Line(drag, ink);
        scene.Line(lift, ink);
        scene.Line(moment, ink);
        scene.Line(top_tr, ink);
        scene.Line(bot_tr, ink);
    }
    return scene;
}

/* ---------------------------------------------------------------- */
/* Airfoil geometry                                                  */
/* ---------------------------------------------------------------- */

Scene AirfoilPlot(const std::vector<Point>& contour, const GeometryInfo& info,
                  const AirfoilOptions& options) {
    Scene scene;
    if (contour.size() < 2) return scene;
    double xmin = contour[0].x, xmax = xmin, ymin = contour[0].y, ymax = ymin;
    for (const Point& p : contour) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
    }
    const double span = std::max(xmax - xmin, 1e-6);
    const double del = 0.1 * std::pow(10.0, std::round(std::log10(span)));
    const double ax0 = std::floor(xmin / del + 1e-6) * del - del;
    const double ax1 = std::ceil(xmax / del - 1e-6) * del + del;
    double ay0 = std::floor(ymin / del - 1e-6) * del;
    double ay1 = std::ceil(ymax / del + 1e-6) * del;
    if (ay1 <= ay0) ay1 = ay0 + del;

    constexpr double kLeft = 96.0;
    constexpr double kRight = 1224.0;
    const double scale = (kRight - kLeft) / (ax1 - ax0);
    const double bottom = Top(880.0);
    const double top = bottom + (ay1 - ay0) * scale;
    auto X = [&](double x) { return kLeft + (x - ax0) * scale; };
    auto Y = [&](double y) { return bottom + (y - ay0) * scale; };

    constexpr double kNumber = 14.0;
    const int d = DecimalsFor(del);
    if (options.grid) {
        for (double x : Steps(ax0, ax1, del))
            if (x > ax0 + 1e-9 && x < ax1 - 1e-9)
                scene.Segment(X(x), bottom, X(x), top, Ink::Fore, Dash::Dotted);
        for (double y : Steps(ay0, ay1, del))
            if (y > ay0 + 1e-9 && y < ay1 - 1e-9)
                scene.Segment(kLeft, Y(y), kRight, Y(y), Ink::Fore,
                              Dash::Dotted);
    }
    Frame(scene, {kLeft, kRight, bottom, top});
    for (double x : Steps(ax0, ax1, del))
        scene.Text(X(x), Top(906.0), kNumber, Fixed(x, d), Ink::Fore, 0.5);
    for (double y : Steps(ay0, ay1, del))
        TickLabel(scene, 88.0, Y(y), kNumber, Fixed(y, d), 1.0);

    std::vector<Point> page;
    page.reserve(contour.size());
    for (const Point& p : contour) page.push_back({X(p.x), Y(p.y)});
    scene.Line(page);

    if (options.ticks) {
        const std::size_t n = page.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Point& a = page[i == 0 ? 0 : i - 1];
            const Point& b = page[i + 1 < n ? i + 1 : n - 1];
            const double tx = b.x - a.x;
            const double ty = b.y - a.y;
            const double t = std::hypot(tx, ty);
            if (t <= 0.0) continue;
            const double nx = -ty / t * 5.0;
            const double ny = tx / t * 5.0;
            scene.Segment(page[i].x - nx, page[i].y - ny, page[i].x + nx,
                          page[i].y + ny);
        }
    }

    if (options.parameters) {
        constexpr double kText = 13.5;
        scene.Text(202.0, top + 153.0, kText, info.name);
        const LegendColumns columns{155.0, 178.0, 202.0};
        const struct {
            const char* label;
            std::string value;
        } rows[] = {{"area", Fixed(info.area, 5)},
                    {"thick.", Fixed(info.thickness, 5)},
                    {"camber", Fixed(info.camber, 5)},
                    {"r_{LE}", Fixed(info.le_radius, 5)},
                    {"Δθ_{TE}", Fixed(info.te_angle, 2) + "°"}};
        double baseline = top + 129.0;
        for (const auto& row : rows) {
            scene.Text(columns.label, baseline, kText, row.label, Ink::Fore,
                       1.0);
            scene.Text(columns.equals, baseline, kText, "=");
            scene.Text(columns.value, baseline, kText, row.value);
            baseline -= 23.5;
        }
    }
    return scene;
}

}  // namespace fauxfoil
