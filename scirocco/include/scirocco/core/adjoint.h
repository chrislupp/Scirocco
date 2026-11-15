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
#pragma once

#include "scirocco/data/sensitivity.h"
#include "scirocco/data/state.h"

namespace scirocco {

/**
 * @brief Compute sensitivities for the inviscid problem using the adjoint
 * method. Must be called after SolveInviscid() and CalcForce().
 *
 * Objectives: cl, cm, cdpi
 * Design variables: foil.x [2 x N] panel node coordinates
 *
 * The adjoint equation is: AIC^T * lambda = (dJ/dgamma)^T
 * Total sensitivity: dJ/dx = dJ/dx_partial - lambda^T * d(A*gamma - rhs)/dx
 * @param state foil state containing the converged inviscid solution
 * Returns shape and scalar derivatives, with the lift constraint applied
 * when givencl is set.
 */
Sensitivity InviscidSensitivities(FoilState& state);

/**
 * @brief Compute force partial derivatives w.r.t. panel coordinates (holding
 * gamma fixed) Populates: dcl_dx_partial, dcm_dx_partial, dcdpi_dx_partial in
 * state.post
 * @param state foil state to update with partial derivative fields
 */
void CalcForceLin(FoilState& state);

/**
 * @brief Compute sensitivities of cl, cd, cm w.r.t. angle of attack (alpha).
 * Inviscid fixed-angle helper: uses stored gamref and cl_ue/cm_ue/cdpi_ue.
 * Use InviscidSensitivities for complete lift-constrained derivatives.
 * @param state foil state containing the converged inviscid solution
 * @param sens Sensitivity struct to populate with dcl_dalpha, dcm_dalpha,
 * dcd_dalpha
 */
void AlphaSensitivities(FoilState& state, Sensitivity& sens);

/**
 * @brief Compute sensitivities of cl, cd, cm w.r.t. Mach number.
 * Inviscid fixed-angle helper: differentiates KT compressibility analytically.
 * Use InviscidSensitivities for complete lift-constrained derivatives.
 * @param state foil state containing the converged inviscid solution
 * @param sens Sensitivity struct to populate with dcl_dMa, dcm_dMa, dcd_dMa
 */
void MachSensitivities(FoilState& state, Sensitivity& sens);

/**
 * @brief Solve the viscous adjoint transpose system R_V^T · λ = rhs.
 *
 * Requires that the viscous primal has been converged and that
 * StoreAdjointFactorization has populated state.glob.coupled_factorization.
 * If it is null this call rebuilds it (idempotent). Returns λ ∈ ℝ^{4·Nsys}.
 * @param state foil state with a valid adjoint factorization
 * @param rhs right-hand side vector of size 4*Nsys
 * @return adjoint solution vector λ of size 4*Nsys
 */
VectorXd SolveViscousAdjoint(FoilState& state, const VectorXd& rhs);

/**
 * @brief Compute sensitivities for the viscous problem using the adjoint
 * method. Must be called after SolveViscous() and CalcForce().
 *
 * Objectives: cl, cm, cd (Squire-Young total drag)
 * Design variables: foil.x [2 x N] panel node coordinates
 *
 * Adjoint formulation:
 *   R_V^T λ_J = ∂J/∂V
 *   dJ/dx_foil = ∂J/∂x_partial - λ_J^T · ∂R_V/∂x_foil
 * where the VJP λ_J^T · ∂R_V/∂x_foil = λ_bl^T · DRbl/∂x + ComputeRueVJP(λ_ue).
 * @param state foil state containing the converged viscous solution
 * @return Sensitivity struct populated with dcl_dx, dcm_dx, dcd_dx and scalar
 * sensitivities
 */
Sensitivity ViscousSensitivities(FoilState& state);

/** Schur complement of Cl(V, alpha, p) = Cl_target.
 * wake_seed contains fixed-final-alpha output derivatives with respect to the
 * input angle that generated a frozen wake, in cl, cd, cm, cdpi order. */
void ApplyLiftConstraint(Sensitivity& sens,
                         const Vector4d& wake_seed = Vector4d::Zero());

/** @brief Total residual derivative with respect to prescribed alpha.
 * Includes circulation, wake tracing, and source-coupling chains, at fixed U.
 * Returns a vector of size 4*Nsys. */
VectorXd DRVDAlpha(FoilState& state);

/** @brief Analytic residual Mach derivative through BL closures and
 * thermodynamics. The h argument is retained for source compatibility and is
 * unused. */
VectorXd DRVDMa(FoilState& state, double h = 1e-7);

/** @brief Analytic residual Reynolds derivative through viscosity and Ret.
 * The h argument is retained for source compatibility and is unused. */
VectorXd DRVDRe(FoilState& state, double h = 1e-7);

}  // namespace scirocco
