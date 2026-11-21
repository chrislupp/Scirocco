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
#include "bl_profiles.h"

#include <algorithm>
#include <cmath>

namespace fauxfoil {

namespace {

/* Pohlhausen quartic: u/Ue = 2e - 2e^3 + e^4 + (L/6) e (1-e)^3 */
double PohlhausenDsOverDelta(double lambda) { return 0.3 - lambda / 120.0; }
double PohlhausenThOverDelta(double lambda) {
    return 37.0 / 315.0 - lambda / 945.0 - lambda * lambda / 9072.0;
}
double PohlhausenH(double lambda) {
    return PohlhausenDsOverDelta(lambda) / PohlhausenThOverDelta(lambda);
}

/// Pressure-gradient parameter whose profile has shape factor @p hk. H
/// falls as lambda rises over [-30, 12].
double LambdaFor(double hk) {
    double lo = -30.0;
    double hi = 12.0;
    if (hk >= PohlhausenH(lo)) return lo;
    if (hk <= PohlhausenH(hi)) return hi;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (PohlhausenH(mid) > hk)
            lo = mid;
        else
            hi = mid;
    }
    return 0.5 * (lo + hi);
}

/// Per-node boundary-layer values, gathered from the surface stations.
struct NodeValues {
    std::vector<const BlStation*> station;  ///< Null where none
    std::vector<int> side;                  ///< 1 upper, 0 lower
};

NodeValues Nodes(const Solution& s) {
    NodeValues v;
    v.station.assign(s.airfoil.size(), nullptr);
    v.side.assign(s.airfoil.size(), 1);
    for (const BlStation& st : s.bl_top)
        if (st.node >= 0 &&
            static_cast<std::size_t>(st.node) < v.station.size())
            v.station[st.node] = &st;
    for (const BlStation& st : s.bl_bottom)
        if (st.node >= 0 &&
            static_cast<std::size_t>(st.node) < v.station.size()) {
            v.station[st.node] = &st;
            v.side[st.node] = 0;
        }
    return v;
}

/// Arc length at each node, from the upper trailing edge (the last node in
/// Scirocco's clockwise order).
std::vector<double> ArcFromUpperTe(const Solution& s) {
    const std::size_t n = s.airfoil.size();
    std::vector<double> arc(n, 0.0);
    for (std::size_t k = n - 1; k-- > 0;) {
        const Point& a = s.airfoil[k + 1];
        const Point& b = s.airfoil[k];
        arc[k] = arc[k + 1] + std::hypot(b.x - a.x, b.y - a.y);
    }
    return arc;
}

bool Build(const Solution& s, std::size_t i, double w, double weight,
           BlProfile& profile) {
    const NodeValues nodes = Nodes(s);
    const std::size_t n = s.airfoil.size();
    const std::size_t j = std::min(i + 1, n - 1);
    const BlStation* a = nodes.station[i];
    const BlStation* b = nodes.station[j];
    if (a == nullptr || b == nullptr) return false;
    auto mix = [w](double p, double q) { return (1.0 - w) * p + w * q; };
    profile.wall = {mix(s.airfoil[i].x, s.airfoil[j].x),
                    mix(s.airfoil[i].y, s.airfoil[j].y)};
    profile.ds = mix(a->ds, b->ds);
    profile.ret = mix(a->ret, b->ret);
    profile.hk = mix(a->hk, b->hk);
    profile.turbulent = w < 0.5 ? a->turbulent : b->turbulent;
    const double ue = mix(a->ue, b->ue);

    /* Wall direction: tangent along the node order, normal outward */
    const std::size_t ia = i == j ? (i > 0 ? i - 1 : 0) : i;
    double tx = s.airfoil[j].x - s.airfoil[ia].x;
    double ty = s.airfoil[j].y - s.airfoil[ia].y;
    const double t = std::hypot(tx, ty);
    if (t <= 0.0) return false;
    tx /= t;
    ty /= t;
    const double nx = -ty;
    const double ny = tx;
    /* Downstream runs along the node order on the upper surface (towards
       its trailing edge) and against it on the lower */
    const double flow = nodes.side[i] == 1 ? 1.0 : -1.0;

    double delta_over_ds = 1.0;
    ProfileU(0.5, profile.hk, profile.turbulent, &delta_over_ds);
    const double delta = delta_over_ds * profile.ds;
    profile.normal = {
        profile.wall,
        {profile.wall.x + nx * delta, profile.wall.y + ny * delta}};
    profile.curve.clear();
    constexpr int kSteps = 24;
    for (int k = 0; k <= kSteps; ++k) {
        const double eta = static_cast<double>(k) / kSteps;
        const double u = ue * ProfileU(eta, profile.hk, profile.turbulent);
        const double along = flow * weight * u;
        profile.curve.push_back(
            {profile.wall.x + nx * eta * delta + tx * along,
             profile.wall.y + ny * eta * delta + ty * along});
    }
    return true;
}

}  // namespace

double ProfileU(double eta, double hk, bool turbulent, double* delta_over_ds) {
    eta = std::clamp(eta, 0.0, 1.0);
    if (turbulent) {
        const double h = std::clamp(hk, 1.05, 4.0);
        const double n = 2.0 / (h - 1.0);
        if (delta_over_ds != nullptr) *delta_over_ds = n + 1.0;
        return std::pow(eta, 1.0 / n);
    }
    const double lambda = LambdaFor(hk);
    if (delta_over_ds != nullptr)
        *delta_over_ds = 1.0 / PohlhausenDsOverDelta(lambda);
    const double m = 1.0 - eta;
    return 2.0 * eta - 2.0 * eta * eta * eta + eta * eta * eta * eta +
           lambda / 6.0 * eta * m * m * m;
}

double ArcLength(const Solution& solution) {
    if (solution.airfoil.size() < 2) return 0.0;
    return ArcFromUpperTe(solution).front();
}

bool ProfileAtArc(const Solution& solution, double s, double weight,
                  BlProfile& profile) {
    const std::size_t n = solution.airfoil.size();
    if (!solution.viscous || n < 3) return false;
    const std::vector<double> arc = ArcFromUpperTe(solution);
    /* arc decreases with the node index */
    s = std::clamp(s, 0.0, arc.front());
    std::size_t i = 0;
    while (i + 1 < n && arc[i + 1] > s) ++i;
    const double span = arc[i] - arc[std::min(i + 1, n - 1)];
    const double w = span > 0.0 ? (arc[i] - s) / span : 0.0;
    return Build(solution, i, std::clamp(w, 0.0, 1.0), weight, profile);
}

bool ProfileNear(const Solution& solution, const Point& point, double weight,
                 BlProfile& profile) {
    const std::size_t n = solution.airfoil.size();
    if (!solution.viscous || n < 3) return false;
    double best = 1e300;
    std::size_t best_i = 0;
    double best_w = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const Point& a = solution.airfoil[i];
        const Point& b = solution.airfoil[i + 1];
        const double dx = b.x - a.x;
        const double dy = b.y - a.y;
        const double l2 = dx * dx + dy * dy;
        double w =
            l2 > 0.0 ? ((point.x - a.x) * dx + (point.y - a.y) * dy) / l2 : 0.0;
        w = std::clamp(w, 0.0, 1.0);
        const double d =
            std::hypot(a.x + w * dx - point.x, a.y + w * dy - point.y);
        if (d < best) {
            best = d;
            best_i = i;
            best_w = w;
        }
    }
    return Build(solution, best_i, best_w, weight, profile);
}

void DrawProfile(Scene& scene, const BlProfile& profile) {
    if (!scene.foil.set) return;
    auto page = [&](const std::vector<Point>& points) {
        std::vector<Point> out;
        out.reserve(points.size());
        for (const Point& p : points) out.push_back(scene.foil.ToPage(p));
        return out;
    };
    scene.Line(page(profile.normal), Ink::Green);
    scene.Line(page(profile.curve), Ink::Green);
}

}  // namespace fauxfoil
