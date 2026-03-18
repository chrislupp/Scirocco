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
#include "cppfoil/core/transition.h"

#include <cmath>
#include <cassert>
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/core/boundary_layer.h"

namespace cppfoil {

int MarchAmplification(FoilState& state, int si) {
    // Marches amplification equation on surface si
    // Returns: index of last laminar station before transition (in Is array)
    // Also updates state.glob.U[2,:] with amplification at laminar stations

    const auto& Is = state.vsol.Is[si];
    int N = static_cast<int>(Is.size());
    if (N == 0) return -1;

    Param param = BuildParam(state, si);
    param.turb = false;
    param.wake = false;

    // No amplification at first station
    state.glob.U(2, Is[0]) = 0.0;

    int i = 1;
    while (i < N) {
        int ip = Is[i - 1];
        int ii = Is[i];

        Vector4d U1 = state.glob.U.col(ip);
        Vector4d U2 = state.glob.U.col(ii);

        // Initialize amp if currently turbulent
        if (ii < static_cast<int>(state.vsol.turb.size()) && state.vsol.turb[ii]) {
            U2(2) = U1(2) * 1.01;
        }

        double dx = state.isol.xi(ii) - state.isol.xi(ip);

        // Newton iterations for amplification equation
        int nNewton = 20;
        for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
            // Amplification rate, averaged
            StationParam(state, param, ip);
            auto [damp1, damp1_U1] = GetDamp(U1, param);
            StationParam(state, param, ii);
            auto [damp2, damp2_U2] = GetDamp(U2, param);

            double damp = 0.5 * (damp1 + damp2);

            double Ramp = U2(2) - U1(2) - damp * dx;

            if (std::abs(Ramp) < 1e-12) break;  // Converged

            // Jacobian: dRamp/dU2[2] = 1 - 0.5 * damp2_U2[2] * dx
            double Ramp_U2 = 1.0 - 0.5 * damp2_U2(2) * dx;
            double dU = -Ramp / Ramp_U2;

            double omega = 1.0;
            double dmax = 0.5 * (1.01 - static_cast<double>(iNewton) / nNewton);
            if (std::abs(dU) > dmax) {
                omega = dmax / std::abs(dU);
            }
            U2(2) += omega * dU;
        }

        // Check for transition
        if (U2(2) > param.ncrit) {
            // Transition found - return last laminar station
            return i - 1;
        } else {
            // Store amplification
            state.glob.U(2, ii) = U2(2);
        }

        ++i;
    }

    // No transition - return last station as laminar
    return N - 1;
}

void UpdateTransition(FoilState& state) {
    // Updates transition location using current state
    // Updates turb flags and U[2,:] (amp/ctau) as needed

    for (int si = 0; si < 2; ++si) {  // Lower and upper only
        const auto& Is = state.vsol.Is[si];
        int N = static_cast<int>(Is.size());
        if (N == 0) continue;

        Param param = BuildParam(state, si);

        // Find current last laminar station (ilam0)
        int ilam0 = -1;
        for (int idx = 0; idx < N; ++idx) {
            int j = Is[idx];
            if (j < static_cast<int>(state.vsol.turb.size()) && state.vsol.turb[j]) {
                ilam0 = idx - 1;
                break;
            }
        }
        if (ilam0 < 0) ilam0 = N - 1;  // All laminar

        // Save current sa values (so we don't change them unnecessarily)
        VectorXd sa(N);
        for (int idx = 0; idx < N; ++idx) {
            sa(idx) = state.glob.U(2, Is[idx]);
        }

        // March amplification to get new last laminar station
        int ilam = MarchAmplification(state, si);

        if (ilam == ilam0) {
            // No change - restore original sa
            for (int idx = 0; idx < N; ++idx) {
                state.glob.U(2, Is[idx]) = sa(idx);
            }
            continue;
        }

        if (ilam < ilam0) {
            // Transition is now earlier: fill in turb between [ilam+1, ilam0]
            param.turb = true;
            StationParam(state, param, Is[ilam + 1]);
            auto [sa0, temp] = GetCttr(state.glob.U.col(Is[ilam + 1]), param);
            double sa1 = (ilam0 < N - 1) ? state.glob.U(2, Is[ilam0 + 1]) : sa0;

            double xi_start = state.isol.xi(Is[ilam + 1]);
            double xi_end = state.isol.xi(Is[std::min(ilam0 + 1, N - 1)]);
            double dx = xi_end - xi_start;

            for (int idx = ilam + 1; idx <= ilam0; ++idx) {
                double f = 0.0;
                if (dx != 0.0 && idx != ilam + 1) {
                    f = (state.isol.xi(Is[idx]) - xi_start) / dx;
                }
                if ((ilam + 1) == ilam0) f = 1.0;

                state.glob.U(2, Is[idx]) = sa0 + f * (sa1 - sa0);
                assert(state.glob.U(2, Is[idx]) > 0 && "negative ctau in UpdateTransition");
                state.vsol.turb[Is[idx]] = true;
            }
        } else {  // ilam > ilam0
            // Transition is now later: mark stations as laminar
            // Amplification already filled in by march_amplification
            for (int idx = ilam0; idx <= ilam && idx < N; ++idx) {
                if (Is[idx] < static_cast<int>(state.vsol.turb.size())) {
                    state.vsol.turb[Is[idx]] = false;
                }
            }
        }
    }
}

void StoreTransition(FoilState& state, int si, int i) {
    // Store xi and x transition locations using current vsol.xt
    // si: surface index (0=lower, 1=upper)
    // i: station index in Is array (transition is between i-1 and i)

    if (si < 2 && i > 0) {
        double xt = state.vsol.xt;  // Current transition xi location

        const auto& Is = state.vsol.Is[si];
        int i0 = Is[i - 1];  // Pre-transition node
        int i1 = Is[i];      // Post-transition node

        double xi0 = state.isol.xi(i0);
        double xi1 = state.isol.xi(i1);

        // Can only store transition on airfoil
        if (i0 < state.foil.N && i1 < state.foil.N) {
            double x0 = state.foil.x(0, i0);
            double x1 = state.foil.x(0, i1);

            state.vsol.Xt(si, 0) = xt;  // xi location
            // Interpolate x location
            if (std::abs(xi1 - xi0) > 1e-12) {
                state.vsol.Xt(si, 1) = x0 + (xt - xi0) / (xi1 - xi0) * (x1 - x0);
            } else {
                state.vsol.Xt(si, 1) = x0;
            }
        }
    }

    // Mark subsequent stations as turbulent
    const auto& Is = state.vsol.Is[si];
    bool found = false;
    for (size_t idx = 0; idx < Is.size(); ++idx) {
        if (static_cast<int>(idx) == i) found = true;
        if (found && Is[idx] < static_cast<int>(state.vsol.turb.size())) {
            state.vsol.turb[Is[idx]] = true;
        }
    }
}

}  // namespace cppfoil
