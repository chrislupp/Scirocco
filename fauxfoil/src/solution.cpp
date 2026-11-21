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
#include "solution.h"

#include <algorithm>
#include <cmath>

#include "airfoil.h"
#include "scirocco/core/boundary_layer.h"
#include "scirocco/data/state.h"
#include "scirocco/properties/bl_properties.h"

namespace fauxfoil {

using scirocco::FoilState;
using scirocco::Vector4d;

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Outward unit normal at node @p i of a clockwise contour.
Point Normal(const std::vector<Point>& x, std::size_t i) {
    const Point& a = x[i == 0 ? 0 : i - 1];
    const Point& b = x[i + 1 < x.size() ? i + 1 : x.size() - 1];
    const double tx = b.x - a.x;
    const double ty = b.y - a.y;
    const double t = std::hypot(tx, ty);
    if (t <= 0.0) return {0.0, 0.0};
    return {-ty / t, tx / t};
}

/// Cp from an incompressible edge speed, with the solver's compressibility
/// correction.
double CpFromUe(const FoilState& state, double ue) {
    return scirocco::GetCpFromUe(std::fabs(ue), state.param).value;
}

}  // namespace

Solution Snapshot(const FoilState& state, double xtr_top, double xtr_bot) {
    Solution s;
    s.name = state.geom.name;
    s.alpha = state.oper.alpha;
    s.cl = state.post.cl;
    s.cm = state.post.cm;
    s.cd = state.post.cd;
    s.cdp = state.post.cdpi;
    s.re = state.oper.Re;
    s.mach = state.oper.Ma;
    s.ncrit = state.param.ncrit;
    s.xtr_top = xtr_top;
    s.xtr_bot = xtr_bot;

    const int n = state.foil.N;
    for (int i = 0; i < n; ++i)
        s.airfoil.push_back({state.foil.x(0, i), state.foil.x(1, i)});
    for (int i = 0; i < n && i < state.post.cp.size(); ++i)
        s.cp.push_back(state.post.cp(i));

    const bool viscous = state.oper.viscous && state.glob.U.cols() > n &&
                         state.post.ds.size() == state.glob.U.cols();
    s.viscous = viscous;
    if (!viscous) {
        /* Surfaces split at the stagnation point, for completeness */
        const int split = state.isol.Istag[1];
        for (int i = split; i < n; ++i) s.upper.push_back(i);
        for (int i = split - 1; i >= 0; --i) s.lower.push_back(i);
        return s;
    }

    for (int i = 0; i < n && i < state.post.cpi.size(); ++i)
        s.cpi.push_back(state.post.cpi(i));
    s.upper.assign(state.vsol.Is[1].begin(), state.vsol.Is[1].end());
    s.lower.assign(state.vsol.Is[0].begin(), state.vsol.Is[0].end());

    /* Wake points and Cp */
    for (int k = 0; k < state.wake.N; ++k) {
        const int j = n + k;
        s.wake.push_back({state.wake.x(0, k), state.wake.x(1, k)});
        s.cp_wake.push_back(CpFromUe(state, state.glob.U(3, j)));
        if (j < state.post.uei.size())
            s.cpi_wake.push_back(CpFromUe(state, state.post.uei(j)));
    }

    /* Boundary-layer stations */
    std::vector<BlStation>* lists[3] = {&s.bl_bottom, &s.bl_top, &s.bl_wake};
    double wake_s0 = 0.0;
    for (int si : {1, 0, 2}) {
        scirocco::Param param = scirocco::BuildParam(state, si);
        Point previous{};
        double s_wake = wake_s0;
        for (std::size_t k = 0; k < state.vsol.Is[si].size(); ++k) {
            const int j = state.vsol.Is[si][k];
            scirocco::StationParam(state, param, j);
            const Vector4d U = state.glob.U.col(j);
            BlStation st;
            const Point at =
                j < n ? s.airfoil[j]
                      : Point{state.wake.x(0, j - n), state.wake.x(1, j - n)};
            st.x = at.x;
            if (si == 2) {
                if (k > 0)
                    s_wake += std::hypot(at.x - previous.x, at.y - previous.y);
                st.s = s_wake;
            } else {
                st.s = j < state.isol.xi.size() ? state.isol.xi(j) : 0.0;
            }
            previous = at;
            st.hk = state.post.Hk.size() > j ? state.post.Hk(j) : 0.0;
            st.ds = state.post.ds(j);
            st.th = state.post.th(j);
            st.ue =
                state.post.ue.size() > j ? std::fabs(state.post.ue(j)) : 0.0;
            st.cf = state.post.cf.size() > j ? state.post.cf(j) : 0.0;
            st.ret = state.post.Ret.size() > j ? state.post.Ret(j) : 0.0;
            st.sa = state.post.sa(j);
            st.turbulent =
                static_cast<std::size_t>(j) < state.vsol.turb.size() &&
                state.vsol.turb[j];
            if (si == 2) st.turbulent = true;
            st.node = j < n ? j : -1;
            st.cd = scirocco::GetCDi(U, param).value;
            st.cteq = st.turbulent ? scirocco::GetCteq(U, param).value : 0.0;
            lists[si]->push_back(st);
        }
        if (si == 0 && !s.bl_bottom.empty()) wake_s0 = s.bl_bottom.back().s;
    }

    /* Displacement surfaces */
    for (int si : {1, 0}) {
        auto& out = si == 1 ? s.dstar_top : s.dstar_bottom;
        for (int j : state.vsol.Is[si]) {
            const Point nrm = Normal(s.airfoil, static_cast<std::size_t>(j));
            const double ds = state.post.ds(j);
            out.push_back(
                {s.airfoil[j].x + nrm.x * ds, s.airfoil[j].y + nrm.y * ds});
        }
    }
    /* Wake displacement edges, drawn as XFOIL's CPDISP draws them: each
       edge starts at its side's trailing-edge displacement point, and the
       wake delta* (which carries both sides plus the trailing-edge gap) is
       shared out in proportion to each side's delta* plus half the gap.
       Near the trailing edge the surface and wake normals differ, so where
       the boundary layer is thick (separated flow at high alpha) the first
       wake edge points can fall behind the trailing-edge displacement
       point; they are left out so the edge does not fold back on itself. */
    const int te_top = state.vsol.Is[1].back();
    const int te_bot = state.vsol.Is[0].back();
    const double gap = std::hypot(s.airfoil[n - 1].x - s.airfoil[0].x,
                                  s.airfoil[n - 1].y - s.airfoil[0].y);
    const double ds_wake = state.post.ds(n);
    double share_top = 0.5;
    double share_bot = 0.5;
    if (ds_wake != 0.0) {
        share_top = (state.post.ds(te_top) + 0.5 * gap) / ds_wake;
        share_bot = (state.post.ds(te_bot) + 0.5 * gap) / ds_wake;
    }
    if (!s.dstar_top.empty()) s.wake_upper.push_back(s.dstar_top.back());
    if (!s.dstar_bottom.empty()) s.wake_lower.push_back(s.dstar_bottom.back());
    for (std::size_t k = 1; k < s.wake.size(); ++k) {
        const Point& a = s.wake[k - 1];
        const Point& b = s.wake[k + 1 < s.wake.size() ? k + 1 : k];
        double tx = b.x - a.x;
        double ty = b.y - a.y;
        const double t = std::hypot(tx, ty);
        if (t <= 0.0) continue;
        tx /= t;
        ty /= t;
        const double ds = state.post.ds(n + static_cast<int>(k));
        const Point upper{s.wake[k].x - ty * ds * share_top,
                          s.wake[k].y + tx * ds * share_top};
        const Point lower{s.wake[k].x + ty * ds * share_bot,
                          s.wake[k].y - tx * ds * share_bot};
        if (s.wake_upper.empty() || upper.x > s.wake_upper.back().x)
            s.wake_upper.push_back(upper);
        if (s.wake_lower.empty() || lower.x > s.wake_lower.back().x)
            s.wake_lower.push_back(lower);
    }
    return s;
}

GeometryInfo Geometry(const std::string& name,
                      const std::vector<Point>& contour) {
    GeometryInfo info;
    info.name = name;
    const std::size_t n = contour.size();
    if (n < 5) return info;

    double area = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Point& a = contour[i];
        const Point& b = contour[(i + 1) % n];
        area += a.x * b.y - b.x * a.y;
    }
    info.area = 0.5 * std::fabs(area);

    MatrixXd x(2, static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i) {
        x(0, static_cast<Eigen::Index>(i)) = contour[i].x;
        x(1, static_cast<Eigen::Index>(i)) = contour[i].y;
    }
    const ThicknessCamber tc = GetThicknessCamber(x);
    info.thickness = tc.thickness;
    info.camber = tc.camber;

    /* Leading-edge radius: circle through the leading-edge node and its
       neighbours */
    const int le = LeadingEdgeIndex(x);
    if (le > 0 && le + 1 < static_cast<int>(n)) {
        const Point& a = contour[le - 1];
        const Point& b = contour[le];
        const Point& c = contour[le + 1];
        const double ab = std::hypot(b.x - a.x, b.y - a.y);
        const double bc = std::hypot(c.x - b.x, c.y - b.y);
        const double ca = std::hypot(a.x - c.x, a.y - c.y);
        const double cross =
            std::fabs((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));
        if (cross > 0.0) info.le_radius = ab * bc * ca / (2.0 * cross);
    }

    /* Trailing-edge angle between the last panels of each surface */
    const Point& l0 = contour[0];
    const Point& l1 = contour[1];
    const Point& u0 = contour[n - 1];
    const Point& u1 = contour[n - 2];
    const double al = std::atan2(l1.y - l0.y, l1.x - l0.x);
    const double au = std::atan2(u1.y - u0.y, u1.x - u0.x);
    double angle = std::fabs(au - al) * 180.0 / kPi;
    if (angle > 180.0) angle = 360.0 - angle;
    info.te_angle = angle;
    return info;
}

}  // namespace fauxfoil
