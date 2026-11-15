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
 * @file bl_march.cpp
 * @brief Mixed-mode boundary-layer march (see bl_march.h), after MRCHDU in
 * XFOIL's xbl.f.
 */
#include "scirocco/core/bl_march.h"

#include <algorithm>
#include <cmath>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/transition.h"
#include "scirocco/properties/bl_properties.h"

namespace scirocco {

namespace {

constexpr int kStationNewton = 25;
/// Local convergence: largest relative change in theta, delta* (and shear).
constexpr double kStationTolerance = 5e-4;
/// A local solve that ends with changes below this is still usable.
constexpr double kUsableChange = 0.1;
/// Largest relative change per local Newton step.
constexpr double kMaxStep = 0.3;
/// Weight of ue against Hk in the mixed constraint (XFOIL's SENSWT).
constexpr double kSensWeight = 1000.0;
/// Laminar/transition switches allowed per station and march. A station at
/// the edge of the transition interval can test laminar in its turbulent
/// state and turbulent in its laminar one; switching back and forth would
/// reset it to the same values every march.
constexpr int kMaxSwitches = 1;

enum class Kind { Laminar, Transition, Turbulent };

struct Station {
    int surface;  /* 0 lower, 1 upper, 2 wake */
    int prev;     /* Global index of the upstream node */
    int node;     /* Global index of this node */
    int position; /* Index along the surface */
    Vector2d xi;  /* xi at prev and node */
    Vector2d aux; /* Wake gap at prev and node */
};

ResidualLin Residual(FoilState& state, const Station& st, Kind kind) {
    Matrix42d pair;
    pair.col(0) = state.glob.U.col(st.prev);
    pair.col(1) = state.glob.U.col(st.node);
    Param param = BuildParam(state, st.surface);
    if (kind == Kind::Transition) {
        param.turb = false;
        return ResidualTransition(state, param, st.xi, pair, st.aux, false);
    }
    StationParam(state, param, st.node);
    return ResidualStation(param, st.xi, pair, st.aux, false);
}

/// Shear stress for a station entering the transition interval.
double TransitionShear(const FoilState& state, const Station& st) {
    Param param = BuildParam(state, st.surface);
    param.turb = true;
    return GetCttr(state.glob.U.col(st.node), param).value;
}

/**
 * @brief Local Newton for one station (MRCHDU's loop over ITBL).
 * @param sens Ue-Hk constraint weight, carried from station to station
 * @return true if the local solve converged or ended usable
 */
bool SolveStation(FoilState& state, const Station& st, bool upstream_turbulent,
                  double& sens) {
    auto& U = state.glob.U;
    auto& turb = state.vsol.turb;
    const int j = st.node;
    const int jp = st.prev;
    const Vector4d before = U.col(j);
    const bool was_turbulent = turb[j];
    const double ncrit = state.param.ncrit;
    const double hklim = st.surface == 2 ? 1.00005 : 1.02;

    Param hk_param = BuildParam(state, st.surface);
    StationParam(state, hk_param, j);

    Kind kind = Kind::Turbulent;
    if (upstream_turbulent) {
        if (!was_turbulent) U(2, j) = U(2, jp); /* Shear stress */
        turb[j] = true;
    } else {
        kind = was_turbulent ? Kind::Transition : Kind::Laminar;
    }

    const double ueref = U(3, j);
    double hkref = GetHk(U.col(j), hk_param).value;
    int switches = 0;
    double dmax = 1.0;
    for (int it = 0; it < kStationNewton; ++it) {
        /* Transition test on the current station state, as in
         * UpdateTransition (XFOIL's TRCHEK). */
        if (!upstream_turbulent) {
            const double amp = AmplifyInterval(state, jp, j);
            const Kind want = amp > ncrit ? Kind::Transition : Kind::Laminar;
            if (want != kind && switches < kMaxSwitches) {
                ++switches;
                if (want == Kind::Transition) {
                    turb[j] = true;
                    U(2, j) =
                        was_turbulent ? before(2) : TransitionShear(state, st);
                } else {
                    turb[j] = false;
                    U(2, j) = amp;
                    /* Was turbulent, now laminar: take the reference Hk
                     * from upstream. */
                    if (it == 0) hkref = GetHk(U.col(jp), hk_param).value;
                }
                kind = want;
            }
        }

        const ResidualLin res = Residual(state, st, kind);
        const Vector4d u = U.col(j);
        const auto hk = GetHk(u, hk_param);
        if (!res.value.allFinite() || !std::isfinite(hk.value)) return false;

        Matrix4d a;
        a.topRows<3>() = res.lin_U.block<3, 4>(0, 4);
        Vector4d b;
        b.head<3>() = -res.value;

        /* dUe response to a unit dHk sets how ue and Hk are weighted. */
        if (it < 15) {
            Matrix4d c = a;
            c.row(3) = hk.lin.transpose();
            Vector4d z = b;
            z(3) = 1.0;
            const Vector4d dz = c.partialPivLu().solve(z);
            const double sennew = kSensWeight * dz(3) * hkref / ueref;
            if (std::isfinite(sennew))
                sens = it < 5 ? sennew : 0.5 * (sens + sennew);
        }
        a.row(3) = hkref * hk.lin.transpose();
        a(3, 3) += sens / ueref;
        b(3) = -(hkref * hkref * (hk.value / hkref - 1.0) +
                 sens * (u(3) / ueref - 1.0));

        const Vector4d du = a.partialPivLu().solve(b);
        if (!du.allFinite()) return false;
        dmax = std::max(std::abs(du(0) / u(0)), std::abs(du(1) / u(1)));
        if (kind != Kind::Laminar)
            dmax = std::max(dmax, std::abs(du(2) / (10.0 * u(2))));
        const double rlx = dmax > kMaxStep ? kMaxStep / dmax : 1.0;

        Vector4d next = u + rlx * du;
        if (!(next(0) > 0.0) || !(next(1) > 0.0) || !(next(3) > 0.0))
            return false;
        if (kind != Kind::Laminar) next(2) = std::clamp(next(2), 1e-7, 0.3);
        /* Keep Hk above its floor (XFOIL's DSLIM). */
        const auto hk_next = GetHk(next, hk_param);
        if (hk_next.value < hklim && hk_next.lin(1) > 0.0)
            next(1) += (hklim - hk_next.value) / hk_next.lin(1);
        U.col(j) = next;
        if (dmax <= kStationTolerance) return true;
    }
    return dmax <= kUsableChange;
}

/**
 * @brief Replace a failed station by an extrapolation from upstream
 * (MRCHDU's fallback), then set its laminar/turbulent state.
 */
void Extrapolate(FoilState& state, const Station& st, bool upstream_turbulent) {
    auto& U = state.glob.U;
    const int j = st.node;
    const int jp = st.prev;
    if (st.surface < 2) {
        const double r = std::sqrt(st.xi(1) / st.xi(0));
        U(0, j) = U(0, jp) * r;
        U(1, j) = U(1, jp) * r;
    } else {
        U(0, j) = U(0, jp);
        const double ratio = (st.xi(1) - st.xi(0)) / (10.0 * U(1, jp));
        U(1, j) = (U(1, jp) + U(0, j) * ratio) / (1.0 + ratio);
    }
    U(3, j) = U(3, jp);
    if (upstream_turbulent) {
        state.vsol.turb[j] = true;
        U(2, j) = U(2, jp);
        return;
    }
    const double amp = AmplifyInterval(state, jp, j);
    state.vsol.turb[j] = amp > state.param.ncrit;
    U(2, j) = state.vsol.turb[j] ? 0.05 : amp;
}

/**
 * @brief Solve the stagnation pair: theta, delta* and amplification at the
 * first two nodes of a surface, with ue held, from the similarity residual
 * at the first node and the station residual between them (the equations
 * BuildGlobSys assembles there). XFOIL marches its similarity station with
 * ue prescribed in the same way.
 * @return true if converged; the nodes are restored otherwise
 */
bool SolveStagnationPair(FoilState& state, int si, int n0, int n1) {
    auto& U = state.glob.U;
    const Vector4d before0 = U.col(n0);
    const Vector4d before1 = U.col(n1);
    const Vector2d xi(state.isol.xi(n0), state.isol.xi(n1));
    for (int it = 0; it < kStationNewton; ++it) {
        Matrix42d pair;
        pair.col(0) = U.col(n0);
        pair.col(1) = U.col(n1);
        const StagnationResult stag = StagnationState(pair, xi);

        Param param = BuildParam(state, si);
        param.turb = false;
        param.simi = true;
        Matrix42d stag_pair;
        stag_pair.col(0) = stag.Ust;
        stag_pair.col(1) = stag.Ust;
        const ResidualLin r1 =
            ResidualStation(param, Vector2d(stag.xst, stag.xst), stag_pair,
                            Vector2d::Zero(), false);
        const Matrix34d r1_ust =
            r1.lin_U.block<3, 4>(0, 0) + r1.lin_U.block<3, 4>(0, 4);
        const Matrix38d r1_u = r1_ust * stag.Ust_U;

        param.simi = false;
        StationParam(state, param, n1);
        const ResidualLin r2 =
            ResidualStation(param, xi, pair, Vector2d::Zero(), false);

        Eigen::Matrix<double, 6, 6> a;
        Eigen::Matrix<double, 6, 1> f;
        f << r1.value, r2.value;
        if (!f.allFinite()) break;
        const int cols[6] = {0, 1, 2, 4, 5, 6};
        for (int c = 0; c < 6; ++c) {
            a.block<3, 1>(0, c) = r1_u.col(cols[c]);
            a.block<3, 1>(3, c) = r2.lin_U.col(cols[c]);
        }
        const Eigen::Matrix<double, 6, 1> d = a.partialPivLu().solve(-f);
        if (!d.allFinite()) break;
        const double dmax =
            std::max({std::abs(d(0) / U(0, n0)), std::abs(d(1) / U(1, n0)),
                      std::abs(d(3) / U(0, n1)), std::abs(d(4) / U(1, n1))});
        const double rlx = dmax > kMaxStep ? kMaxStep / dmax : 1.0;
        U.col(n0).head<3>() += rlx * d.head<3>();
        U.col(n1).head<3>() += rlx * d.tail<3>();
        if (!(U(0, n0) > 0.0) || !(U(1, n0) > 0.0) || !(U(0, n1) > 0.0) ||
            !(U(1, n1) > 0.0))
            break;
        if (dmax <= kStationTolerance) return true;
    }
    U.col(n0) = before0;
    U.col(n1) = before1;
    return false;
}

}  // namespace

int BoundaryLayerMarch(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    auto& U = state.glob.U;
    auto& turb = state.vsol.turb;
    int failures = 0;
    double sens = 0.0;

    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        const int n = static_cast<int>(Is.size());
        if (n < 3) continue;
        int start = 1;
        if (si == 2) {
            /* First wake station: merge the two trailing edges as they now
             * stand (XFOIL's TESYS). The residual is explicit in it. */
            Param param = BuildParam(state, si);
            const auto wake = WakeSys(state, param);
            if (wake.R.allFinite()) U.col(Is[0]).head<3>() -= wake.R;
        }
        if (si < 2) {
            const bool on_node =
                state.isol.xi(Is[0]) < 1e-8 * state.isol.xi(Is[n - 1]);
            start = (on_node ? 1 : 0) + 2;
            if (!on_node && !turb[Is[0]] && !turb[Is[1]] &&
                !SolveStagnationPair(state, si, Is[0], Is[1])) {
                /* Start the pair over from stagnation flow (Thwaites, as
                 * in InitBoundaryLayer) and solve again. */
                const Vector4d keep0 = U.col(Is[0]);
                const Vector4d keep1 = U.col(Is[1]);
                const double K = U(3, Is[1]) / state.isol.xi(Is[1]);
                const auto [th, ds] =
                    ThwaitesInit(K, state.param.mu0 / state.param.rho0);
                for (int node : {Is[0], Is[1]}) {
                    U(0, node) = th;
                    U(1, node) = ds;
                    U(2, node) = 0.0;
                }
                if (!SolveStagnationPair(state, si, Is[0], Is[1])) {
                    U.col(Is[0]) = keep0;
                    U.col(Is[1]) = keep1;
                    ++failures;
                }
            }
        }
        for (int i = start; i < n; ++i) {
            const int j = Is[i];
            const int jp = Is[i - 1];
            if (j >= Nsys) continue;
            Station st{};
            st.surface = si;
            st.prev = jp;
            st.node = j;
            st.position = i;
            st.xi = Vector2d(state.isol.xi(jp), state.isol.xi(j));
            if (si == 2 && state.vsol.wgap.size() > i)
                st.aux = Vector2d(state.vsol.wgap(i - 1), state.vsol.wgap(i));
            else
                st.aux = Vector2d::Zero();

            const Vector4d before = U.col(j);
            const bool was_turbulent = turb[j];
            const bool upstream_turbulent = turb[jp];
            if (SolveStation(state, st, upstream_turbulent, sens)) continue;
            ++failures;
            if (i > start + 1) Extrapolate(state, st, upstream_turbulent);
            if (i <= start + 1 || !U.col(j).allFinite()) {
                U.col(j) = before;
                turb[j] = was_turbulent || upstream_turbulent;
                if (upstream_turbulent && !was_turbulent) U(2, j) = U(2, jp);
            }
        }
    }
    return failures;
}

}  // namespace scirocco
