/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#ifndef SCIROCCO_CORE_VISCOUS_RESIDUAL_DERIV_H
#define SCIROCCO_CORE_VISCOUS_RESIDUAL_DERIV_H

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

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

// Reverse-mode vector-Jacobian product (VJP) of R_ue w.r.t. foil.x:
//   result[2N] = lambda_ue^T * dR_ue/dfoil_x
//
// where R_ue[i] = U[3, i] - ueinv[i] - sgue[i] * (ue_sigma * phi)[i],
//       phi = sigma_m * (U[1, :] .* U[3, :])      (mass-source contribution)
//       sgue = [sgnue (foil); ones (wake)]
//       ueinv = ue from inviscid solve (foil: sgnue .* gam; wake: uewi).
//
// This is what the viscous adjoint actually consumes (Phase 6's chain through
// the ue-coupling block); the full forward Jacobian would have been [Nsys x 2N]
// dense, much more expensive than this single VJP.
//
// Chain pieces (each contributing -lambda terms via reverse-mode):
//   A. -lambda^T * dueinv/dfoil_x
//        - foil-side (i ∈ [0, N)): -lambda[i] * sgnue[i] * dgam[i]/dfoil_x
//        - wake-side (i ∈ [N, Nsys)): -lambda[i] * duewi[i-N]/dfoil_x
//   B. -lambda_sg^T * dueisigma * phi (panel-velocity derivative chain)
//   C. -lambda_sg^T * ue_sigma * dphi
//      where dphi = dsigma_m * mass (sigma_m has 1/ds_panel entries → chain
//      via dfoil_s_dx and dwake_s_dx).
//
// Returns a [2N] vector in row-major foil.x layout (col = 2k + d).
VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const MatrixXd& dfoil_s_dx);

// Compute d(R_ue)/d(foil.x) [Nsys x 2N] forward Jacobian (FOR VALIDATION ONLY).
// Builds the full dense Jacobian by applying ComputeRueVJP with each row of
// the identity matrix as lambda. Slow (Nsys * Npanel work) — use VJP in
// production. Provided so the Phase 2 reference test can compare against
// comprehensive_reference.h5's dR_ue_dx.
MatrixXd DRueDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_VISCOUS_RESIDUAL_DERIV_H
