// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/cppfoil.h"

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/paneling.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/core/forces.h"

namespace cppfoil {

Cppfoil::Cppfoil(const std::string& naca, int npanel) {
    NacaPoints(state_.geom, naca);
    MakePanels(state_.geom, npanel, state_.foil);
}

Cppfoil::Cppfoil(const MatrixXd& coords, int npanel) {
    SetCoords(state_.geom, coords);
    MakePanels(state_.geom, npanel, state_.foil);
}

void Cppfoil::SetOper(double alpha, double cl, double Re, double Ma) {
    state_.oper.alpha = alpha;
    state_.oper.Ma = Ma;

    if (cl >= 0) {
        state_.oper.cltgt = cl;
        state_.oper.givencl = true;
    } else {
        state_.oper.givencl = false;
    }

    if (Re > 0) {
        state_.oper.Re = Re;
        state_.oper.viscous = true;
    } else {
        state_.oper.viscous = false;
    }

    state_.oper.Vinf = 1.0;  // Normalized
}

void Cppfoil::Solve() {
    if (state_.oper.viscous) {
        // TODO: Implement viscous solver
        // For now, fall back to inviscid
        SolveInviscid(state_);
        CalcForce(state_);
    } else {
        SolveInviscid(state_);
        CalcForce(state_);
    }
}

}  // namespace cppfoil
