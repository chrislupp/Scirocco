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
 * @file coupled_march.cpp
 * @brief Quasi-simultaneous boundary-layer march (see coupled_march.h).
 */
#include "scirocco/core/coupled_march.h"

#include <algorithm>
#include <cmath>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/properties/bl_properties.h"

namespace scirocco {

namespace {

constexpr int kStationNewton = 20;
constexpr double kStationTolerance = 1e-11;
/// Largest relative change of ue, theta or delta* accepted from a station
/// solve. Near the leading edge the interaction term is stiff and the local
/// problem can have spurious roots far from the Newton iterate.
constexpr double kMaxRelativeChange = 0.5;

/// True if @p after stays within kMaxRelativeChange of @p before.
bool CloseTo(const Vector4d& before, const Vector4d& after) {
    for (int k : {0, 1, 3})
        if (std::abs(after(k) - before(k)) > kMaxRelativeChange * before(k))
            return false;
    return true;
}

enum class StationKind { Laminar, Transition, Turbulent };

/// Everything one station solve needs besides the state.
struct Station {
    int surface;  /* 0 lower, 1 upper, 2 wake */
    int prev;     /* Global index of the upstream node */
    int node;     /* Global index of this node */
    Vector2d xi;  /* xi at prev and node */
    Vector2d aux; /* Wake gap at prev and node */
    double uenv;  /* ueinv + induced velocity from all other nodes */
    double self;  /* ue_m(node, node) */
};

/**
 * @brief Newton on the station's three boundary-layer equations plus its
 * edge-velocity equation, upstream node fixed. Updates U(:, node).
 * @return true if converged
 */
bool SolveStation(FoilState& state, const Station& st, StationKind kind) {
    auto& U = state.glob.U;
    Param param = BuildParam(state, st.surface);
    if (kind == StationKind::Transition) {
        param.turb = false;
    } else {
        StationParam(state, param, st.node);
    }
    for (int it = 0; it < kStationNewton; ++it) {
        Matrix42d pair;
        pair.col(0) = U.col(st.prev);
        pair.col(1) = U.col(st.node);
        const ResidualLin res =
            kind == StationKind::Transition
                ? ResidualTransition(state, param, st.xi, pair, st.aux)
                : ResidualStation(param, st.xi, pair, st.aux);
        const Vector4d u = U.col(st.node);
        Vector4d f;
        f.head<3>() = res.value;
        f(3) = u(3) - (st.uenv + st.self * u(1) * u(3));
        if (!f.allFinite()) return false;
        if (f.norm() < kStationTolerance) return true;

        Matrix4d a;
        a.topRows<3>() = res.lin_U.block<3, 4>(0, 4);
        a.row(3) << 0.0, -st.self * u(3), 0.0, 1.0 - st.self * u(1);
        Vector4d du = a.partialPivLu().solve(-f);
        if (!du.allFinite()) return false;

        /* Relative changes limited to 30%, as in the initial march. */
        double dm = std::max({std::abs(du(0)) / u(0), std::abs(du(1)) / u(1),
                              std::abs(du(3)) / u(3)});
        if (kind == StationKind::Laminar)
            dm = std::max(dm, std::abs(du(2)) / 10.0);
        else
            dm = std::max(dm, std::abs(du(2)) / std::max(u(2), 1e-7));
        if (dm > 0.3) du *= 0.3 / dm;

        Vector4d next = u + du;
        if (!(next(0) > 0.0) || !(next(1) > 0.0) || !(next(3) > 0.0))
            return false;
        if (kind != StationKind::Laminar)
            next(2) = std::clamp(next(2), 1e-7, 0.3);
        U.col(st.node) = next;
    }
    return false;
}

}  // namespace

int CoupledMarch(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    auto& U = state.glob.U;
    auto& turb = state.vsol.turb;
    const MatrixXd& ue_m = state.vsol.ue_m;
    const VectorXd ueinv = GetUeInv(state);
    const double ncrit = state.param.ncrit;

    VectorXd m = U.row(1).cwiseProduct(U.row(3)).transpose();
    VectorXd induced = ue_m * m; /* sum_j ue_m(i, j) m_j */
    int failures = 0;

    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        const int n = static_cast<int>(Is.size());
        if (n < 3) continue;
        int start = 1;
        if (si < 2) {
            const bool on_node =
                state.isol.xi(Is[0]) < 1e-8 * state.isol.xi(Is[n - 1]);
            start = (on_node ? 1 : 0) + 2;
        }
        for (int i = start; i < n; ++i) {
            const int j = Is[i];
            const int jp = Is[i - 1];
            if (j >= Nsys) continue;
            Station st{};
            st.surface = si;
            st.prev = jp;
            st.node = j;
            st.xi = Vector2d(state.isol.xi(jp), state.isol.xi(j));
            if (si == 2 && state.vsol.wgap.size() > i)
                st.aux = Vector2d(state.vsol.wgap(i - 1), state.vsol.wgap(i));
            else
                st.aux = Vector2d::Zero();
            st.self = ue_m(j, j);
            st.uenv = ueinv(j) + induced(j) - st.self * m(j);

            const Vector4d before = U.col(j);
            const bool was_turbulent = turb[j];
            const bool upstream_turbulent = turb[jp];
            bool ok = false;
            if (upstream_turbulent) {
                if (!was_turbulent) U(2, j) = U(2, jp); /* Shear stress */
                turb[j] = true;
                ok = SolveStation(state, st, StationKind::Turbulent) &&
                     CloseTo(before, U.col(j));
                if (!ok) {
                    U.col(j) = before;
                    if (!was_turbulent) U(2, j) = U(2, jp);
                }
            } else {
                /* Laminar upstream: solve laminar, and switch to the
                 * transition interval if the amplification passes ncrit. */
                turb[j] = false;
                if (was_turbulent) U(2, j) = U(2, jp); /* Amplification */
                ok = SolveStation(state, st, StationKind::Laminar);
                if (ok && U(2, j) > ncrit) {
                    turb[j] = true;
                    Param param = BuildParam(state, si);
                    param.turb = true;
                    U(2, j) = was_turbulent ? before(2)
                                            : GetCttr(U.col(j), param).value;
                    ok = SolveStation(state, st, StationKind::Transition);
                }
                ok = ok && CloseTo(before, U.col(j));
                if (!ok) {
                    U.col(j) = before;
                    turb[j] = was_turbulent;
                }
            }
            if (!ok) ++failures;

            const double dm = U(1, j) * U(3, j) - m(j);
            if (dm != 0.0) {
                induced += ue_m.col(j) * dm;
                m(j) += dm;
            }
        }
    }
    return failures;
}

}  // namespace scirocco
