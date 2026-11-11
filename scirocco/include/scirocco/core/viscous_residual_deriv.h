/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#pragma once

#include "scirocco/config.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief Phase 2 residual partials of R_V w.r.t. foil.x at fixed BL state U.
 *
 * R_BL_partial [3*Nsys x 2N] dependence on foil.x at fixed U is purely
 * through xi (R_x chain) and wake gap (wgap chain, only at wake stations):
 *
 *   dR_BL/d(foil.x) = R_x · DXiDFoilX  +  R_BL_wgap · DWgapDFoilX
 *
 * where R_x is already assembled by BuildGlobSys (state.glob.R_x) and
 * R_BL_wgap is derived from the existing R_U (because in ResidualStation:
 *   U_modified = U - [0, Aux, 0, 0]   → ∂R/∂Aux = -∂R/∂(U[1])).
 *
 * PREREQUISITES: BuildGlobSys() must have been called (populates R_x and
 * R_U), and SolveViscous() / StagPointFind() must have placed foil.s, xi,
 * wgap, and gam consistent with mfoil's CS perturbation chain.
 *
 * @param state      Fully assembled foil state (R_x and R_U populated).
 * @param dfoil_s_dx Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 * @return dR_BL/d(foil.x) [3*Nsys x 2N].
 */
MatrixXd DRblDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

/**
 * @brief Reverse-mode vector-Jacobian product (VJP) of R_ue w.r.t. foil.x.
 *
 * result[2N] = lambda_ue^T * dR_ue/dfoil_x
 *
 * where R_ue[i] = U[3, i] - ueinv[i] - sgue[i] * (ue_sigma * phi)[i],
 *       phi = sigma_m * (U[1, :] .* U[3, :])      (mass-source contribution)
 *       sgue = [sgnue (foil); ones (wake)]
 *       ueinv = ue from inviscid solve (foil: sgnue .* gam; wake: uewi).
 *
 * This is what the viscous adjoint actually consumes (Phase 6's chain through
 * the ue-coupling block); the full forward Jacobian would have been [Nsys x 2N]
 * dense, much more expensive than this single VJP.
 *
 * Chain pieces (each contributing -lambda terms via reverse-mode):
 *   A. -lambda^T * dueinv/dfoil_x
 *        - foil-side (i ∈ [0, N)): -lambda[i] * sgnue[i] * dgam[i]/dfoil_x
 *        - wake-side (i ∈ [N, Nsys)): -lambda[i] * duewi[i-N]/dfoil_x
 *   B. -lambda_sg^T * dueisigma * phi (panel-velocity derivative chain)
 *   C. -lambda_sg^T * ue_sigma * dphi
 *      where dphi = dsigma_m * mass (sigma_m has 1/ds_panel entries → chain
 *      via dfoil_s_dx and dwake_s_dx).
 *
 * Returns a [2N] vector in row-major foil.x layout (col = 2k + d).
 *
 * This overload builds a RueVJPContext on every call. Callers that apply the
 * VJP to more than one lambda for the same state should build the context
 * once and use the overload below.
 *
 * @param state        Fully assembled foil state.
 * @param lambda_ue    Adjoint vector [Nsys].
 * @param dfoil_s_dx   Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 * @return VJP result [2N].
 */
VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const MatrixXd& dfoil_s_dx);

/**
 * @brief The parts of ComputeRueVJP that depend on the state but not on
 * lambda_ue.
 *
 * Building these (mostly d(gamma)/d(foil.x) and the wake Jacobians) costs far
 * more than applying the VJP, so a gradient builds them once and reuses them
 * for every output. Only valid for the state it was built from.
 */
struct RueVJPContext {
    MatrixXd dfoil_s_dx; /**< d(foil.s)/d(foil.x) [N x 2N] */
    MatrixXd dgam_dx;    /**< d(gamma)/d(foil.x) at oper.alpha [N+1 x 2N] */
    MatrixXd dwake_x_dx; /**< d(wake.x)/d(foil.x) [2*Nw x 2N] */
    MatrixXd dwake_t_dx; /**< d(wake.t)/d(foil.x) [2*Nw x 2N] */
    MatrixXd dwake_s_dx; /**< d(wake.s)/d(foil.x) [Nw x 2N] */
    Eigen::PartialPivLU<MatrixXd> aic_t_lu; /**< LU of AIC^T */
};

/**
 * @brief Build the lambda-independent context for ComputeRueVJP.
 *
 * @param state        Fully assembled foil state.
 * @param dfoil_s_dx   Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 */
RueVJPContext BuildRueVJPContext(const FoilState& state,
                                 const MatrixXd& dfoil_s_dx);

/**
 * @brief ComputeRueVJP with a prebuilt context. Same result as the overload
 * above, without rebuilding the context.
 *
 * @param state     The state ctx was built from.
 * @param lambda_ue Adjoint vector [Nsys].
 * @param ctx       BuildRueVJPContext(state, dfoil_s_dx).
 * @return VJP result [2N].
 */
VectorXd ComputeRueVJP(const FoilState& state, const VectorXd& lambda_ue,
                       const RueVJPContext& ctx);

/**
 * @brief Compute d(R_ue)/d(foil.x) [Nsys x 2N] forward Jacobian (FOR VALIDATION
 * ONLY).
 *
 * Builds the full dense Jacobian by applying ComputeRueVJP with each row of
 * the identity matrix as lambda. Slow (Nsys * Npanel work) — use VJP in
 * production. Provided so the Phase 2 reference test can compare against
 * comprehensive_reference.h5's dR_ue_dx.
 *
 * @param state        Fully assembled foil state.
 * @param dfoil_s_dx   Sensitivity of foil arc-lengths w.r.t. foil.x [N x 2N].
 * @return dR_ue/d(foil.x) [Nsys x 2N].
 */
MatrixXd DRueDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx);

} /* namespace scirocco */
