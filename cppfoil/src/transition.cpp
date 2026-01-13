// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/transition.h"

#include <cmath>
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/core/boundary_layer.h"

namespace cppfoil {

void MarchAmplification(FoilState& state, int si) {
    // March amplification factor along surface si
    // Using eN method

    const auto& Is = state.vsol.Is[si];
    if (Is.empty()) return;

    Param param = BuildParam(state, si);
    double ncrit = param.ncrit;

    for (size_t idx = 1; idx < Is.size(); ++idx) {
        int i = Is[idx];
        int ip = Is[idx - 1];

        StationParam(state, param, i);

        // Get amplification rate
        Vector4d Ui = state.glob.U.col(i);
        auto [damp, damp_U] = GetDamp(Ui, param);

        // Update amplification factor
        double ds = state.foil.s(i) - state.foil.s(ip);
        double dN = damp * ds;
        state.glob.U(2, i) = state.glob.U(2, ip) + dN;

        // Check for transition
        if (state.glob.U(2, i) > ncrit && !state.vsol.turb[i]) {
            StoreTransition(state, si, i);
        }
    }
}

void UpdateTransition(FoilState& state) {
    // Update transition on all surfaces
    for (int si = 0; si < 2; ++si) {  // Lower and upper only
        MarchAmplification(state, si);
    }
}

void StoreTransition(FoilState& state, int si, int i) {
    // Store transition location
    if (si < 2) {
        state.vsol.Xt(si, 0) = state.isol.xi(i);
        state.vsol.Xt(si, 1) = state.foil.x(0, i);
    }

    // Mark subsequent stations as turbulent
    const auto& Is = state.vsol.Is[si];
    bool found = false;
    for (int idx : Is) {
        if (idx == i) found = true;
        if (found && idx < static_cast<int>(state.vsol.turb.size())) {
            state.vsol.turb[idx] = true;
        }
    }
}

}  // namespace cppfoil
