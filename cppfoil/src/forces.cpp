/*
    cppfoil: C++ Airfoil Analysis Library

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
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "cppfoil/core/forces.h"

#include <cmath>

#include "cppfoil/core/boundary_layer.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

namespace {

// Get pressure coefficient from edge velocity
double GetCp(double ue, double Vinf, double KTb, double KTl) {
    double cp_inc = 1.0 - (ue / Vinf) * (ue / Vinf);

    // Karman-Tsien compressibility correction
    if (KTl > 0) {
        double denom = KTb + KTl * (1.0 + KTb) * 0.5 * cp_inc;
        if (std::abs(denom) > 1e-10) {
            return cp_inc / denom;
        }
    }
    return cp_inc;
}

}  // namespace

void CalcForce(FoilState& state) {
    double chord = state.geom.chord;
    Vector2d xref = state.geom.xref;
    double Vinf = state.param.Vinf > 0 ? state.param.Vinf : 1.0;
    double alpha = state.oper.alpha;
    int N = state.foil.N;

    // Get edge velocity
    VectorXd ue;
    if (state.oper.viscous && state.glob.U.cols() > 0) {
        ue = state.glob.U.row(3).transpose();
    } else {
        ue = GetUeInv(state);
    }

    // Calculate pressure coefficient at each node
    state.post.cp.resize(N);
    state.post.cpi.resize(N);

    VectorXd ueinv = GetUeInv(state);

    for (int i = 0; i < N; ++i) {
        state.post.cp(i) = GetCp(ue(i), Vinf, state.param.KTb, state.param.KTl);
        state.post.cpi(i) =
            GetCp(ueinv(i), Vinf, state.param.KTb, state.param.KTl);
    }

    // Integrate forces by panel integration
    double cl = 0.0, cm = 0.0, cdpi = 0.0;
    state.post.cl_ue.resize(N);
    state.post.cl_ue.setZero();
    double cl_alpha = 0.0;
    double cdpi_alpha = 0.0;

    for (int i0 = 1; i0 <= N; ++i0) {
        int i = (i0 == N) ? 0 : i0;
        int ip = (i0 == N) ? N - 1 : i0 - 1;

        Vector2d x1 = state.foil.x.col(ip);
        Vector2d x2 = state.foil.x.col(i);
        Vector2d dxv = x2 - x1;
        Vector2d dx1 = x1 - xref;
        Vector2d dx2 = x2 - xref;

        // Panel contributions (minus from CW node ordering)
        double dx = -dxv(0) * Cosd(alpha) - dxv(1) * Sind(alpha);
        double dz = dxv(1) * Cosd(alpha) - dxv(0) * Sind(alpha);

        double cp1 = state.post.cp(ip);
        double cp2 = state.post.cp(i);
        double cpbar = 0.5 * (cp1 + cp2);

        cl += dx * cpbar;

        // Linearization w.r.t. ue (through cp, with KT correction)
        auto [cp_ip_val, cp_ip_U] = GetCpFromUe(ue(ip), state.param);
        state.post.cl_ue(ip) += dx * 0.5 * cp_ip_U(3);
        auto [cp_i_val, cp_i_U] = GetCpFromUe(ue(i), state.param);
        state.post.cl_ue(i) += dx * 0.5 * cp_i_U(3);

        // Linearization w.r.t. alpha
        cl_alpha +=
            cpbar * (Sind(alpha) * dxv(0) - Cosd(alpha) * dxv(1)) * kDegToRad;

        // cdpi_alpha: d(dz)/dalpha = (-dxv(1)*sina - dxv(0)*cosa) * deg2rad ... wait
        // dz = dxv(1)*cos(a) - dxv(0)*sin(a)
        // d(dz)/da = (-dxv(1)*sin(a) - dxv(0)*cos(a)) * deg2rad
        cdpi_alpha +=
            cpbar * (-dxv(1) * Sind(alpha) - dxv(0) * Cosd(alpha)) * kDegToRad;

        // Moment (about xref)
        double dx1nds = dxv(0) * dx1(0) + dxv(1) * dx1(1);
        double dx2nds = dxv(0) * dx2(0) + dxv(1) * dx2(1);
        cm += cp1 * dx1nds / 3.0 + cp1 * dx2nds / 6.0 + cp2 * dx1nds / 6.0 +
              cp2 * dx2nds / 3.0;

        // Near-field pressure drag
        cdpi += dz * cpbar;
    }

    cl /= chord;
    cm /= (chord * chord);
    cdpi /= chord;
    cl_alpha /= chord;
    cdpi_alpha /= chord;
    state.post.cl_ue /= chord;

    state.post.cl = cl;
    state.post.cl_alpha = cl_alpha;
    state.post.cm = cm;
    state.post.cdpi = cdpi;
    state.post.cdpi_alpha = cdpi_alpha;

    // Viscous drag (if viscous solution exists)
    state.post.cd = 0.0;
    state.post.cdf = 0.0;
    state.post.cdp = 0.0;

    if (state.oper.viscous && state.glob.U.cols() > 0 &&
        !state.vsol.Is.empty()) {
        // Squire-Young relation for total drag (from wake trailing edge)
        if (state.vsol.Is.size() > 2 && !state.vsol.Is[2].empty()) {
            int iw = state.vsol.Is[2].back();  // End of wake
            if (iw < state.glob.U.cols()) {
                Vector4d Uw = state.glob.U.col(iw);
                auto [H, H_U] = GetH(Uw);
                auto [uk, uk_U] = GetUk(Uw(3), state.param);
                double uk_Vinf = uk / Vinf;

                state.post.cd =
                    2.0 * Uw(0) * std::pow(uk_Vinf, (5.0 + H) / 2.0);
            }
        }

        // Skin friction drag integration
        double Df = 0.0;
        double rho = state.oper.rho;
        double qinf = 0.5 * rho * Vinf * Vinf;

        for (int si = 0; si < 2; ++si) {  // Lower and upper surfaces
            if (si >= static_cast<int>(state.vsol.Is.size())) continue;
            const auto& Is = state.vsol.Is[si];
            if (Is.empty()) continue;

            Param param = state.param;
            param.wake = false;

            double cf1 = 0.0;
            double ue1 = 0.0;
            double rho1 = rho;
            Vector2d x1 = state.isol.xstag;

            for (size_t idx = 0; idx < Is.size(); ++idx) {
                int i = Is[idx];
                if (i >= state.glob.U.cols()) continue;

                // Update turbulence flag
                if (i < static_cast<int>(state.vsol.turb.size())) {
                    param.turb = state.vsol.turb[i];
                }

                Vector4d Ui = state.glob.U.col(i);
                auto [cf2, cf2_U] = GetCf(Ui, param);
                auto [ue2, ue2_U] = GetUk(Ui(3), param);
                auto [rho2, rho2_U] = GetRho(Ui, param);

                Vector2d x2 = state.foil.x.col(i);
                Vector2d dxv = x2 - x1;
                double dx_proj = dxv(0) * Cosd(alpha) + dxv(1) * Sind(alpha);

                Df += 0.25 * (rho1 * cf1 * ue1 * ue1 + rho2 * cf2 * ue2 * ue2) *
                      dx_proj;

                cf1 = cf2;
                ue1 = ue2;
                x1 = x2;
                rho1 = rho2;
            }
        }

        state.post.cdf = Df / (qinf * chord);
        state.post.cdp = state.post.cd - state.post.cdf;
    }
}

void GetDistributions(FoilState& state) {
    if (state.glob.U.cols() == 0) {
        return;
    }

    int N = static_cast<int>(state.glob.U.cols());

    state.post.th = state.glob.U.row(0).transpose();
    state.post.ds = state.glob.U.row(1).transpose();
    state.post.sa = state.glob.U.row(2).transpose();

    // Edge velocity (compressible, Karman-Tsien correction)
    state.post.ue.resize(N);
    for (int i = 0; i < N; ++i) {
        auto [uk, uk_U] = GetUk(state.glob.U(3, i), state.param);
        state.post.ue(i) = uk;
    }

    // Inviscid edge velocity
    state.post.uei = GetUeInv(state);

    // Compute cf, Ret, Hk distributions
    state.post.cf.resize(N);
    state.post.cf.setZero();
    state.post.Ret.resize(N);
    state.post.Ret.setZero();
    state.post.Hk.resize(N);
    state.post.Hk.setZero();

    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        Param param = BuildParam(state, si);
        for (int idx = 0; idx < static_cast<int>(Is.size()); ++idx) {
            int j = Is[idx];
            Vector4d Uj = state.glob.U.col(j);
            StationParam(state, param, j);
            auto [uk, uk_U] = GetUk(Uj(3), param);
            auto [cfloc, cf_U] = GetCf(Uj, param);
            state.post.cf(j) = cfloc * uk * uk / (param.Vinf * param.Vinf);
            auto [ret, Ret_U] = GetRet(Uj, param);
            state.post.Ret(j) = ret;
            auto [hk, Hk_U] = GetHk(Uj, param);
            state.post.Hk(j) = hk;
        }
    }
}

}  // namespace cppfoil
