/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#ifndef CPPFOIL_CORE_VISCOUS_RESIDUAL_DERIV_H
#define CPPFOIL_CORE_VISCOUS_RESIDUAL_DERIV_H

#include "cppfoil/config.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// ============================================================================
// Phase 2 residual partials of R_V w.r.t. foil.x at fixed BL state U.
//
// R_BL_partial [3*Nsys x 2N] dependence on foil.x at fixed U is purely
// through xi (R_x chain) and wake gap (wgap chain, only at wake stations):
//
//   dR_BL/d(foil.x) = R_x · DXiDFoilX  +  R_BL_wgap · DWgapDFoilX
//
// where R_x is already assembled by BuildGlobSys (state.glob.R_x) and
// R_BL_wgap is derived from the existing R_U (because in ResidualStation:
//   U_modified = U - [0, Aux, 0, 0]   → ∂R/∂Aux = -∂R/∂(U[1])).
//
// PREREQUISITES: BuildGlobSys() must have been called (populates R_x and
// R_U), and SolveViscous() / StagPointFind() must have placed foil.s, xi,
// wgap, and gam consistent with mfoil's CS perturbation chain.
// ============================================================================
MatrixXd DRblDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

// R_ue [Nsys x 2N] dependence on foil.x at fixed U is through:
//   1. ueinv (foil-side: gamref + sgnue; wake-side: uewiref)
//   2. ue_sigma (mass-source influence at each station)
//   3. sigma_m (mass source assembly)
// All three depend on foil.x (panel positions) AND on wake.x (for
// wake-station rows). The wake.x chain uses DWakeXDFoilX from Phase 1.
//
// (Implementation pending — full chain through panel-influence derivatives.)
MatrixXd DRueDFoilX(const FoilState& state);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_VISCOUS_RESIDUAL_DERIV_H
