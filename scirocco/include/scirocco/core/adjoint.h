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
 * @brief Compute sensitivities for the inviscid problem using the adjoint method.
 * Must be called after SolveInviscid() and CalcForce().
 *
 * Objectives: cl, cm, cdpi
 * Design variables: foil.x [2 x N] panel node coordinates
 *
 * The adjoint equation is: AIC^T * lambda = (dJ/dgamma)^T
 * Total sensitivity: dJ/dx = dJ/dx_partial + lambda^T * d(A*gamma - rhs)/dx
 * @param state foil state containing the converged inviscid solution
 * @return Sensitivity struct populated with dcl_dx, dcm_dx, dcdpi_dx
 */
Sensitivity InviscidSensitivities(FoilState& state);

/**
 * @brief Compute force partial derivatives w.r.t. panel coordinates (holding gamma
 * fixed) Populates: dcl_dx_partial, dcm_dx_partial, dcdpi_dx_partial in
 * state.post
 * @param state foil state to update with partial derivative fields
 */
void CalcForceLin(FoilState& state);

/**
 * @brief Compute sensitivities of cl, cd, cm w.r.t. angle of attack (alpha).
 * Inviscid only: uses stored gamref and cl_ue/cm_ue/cdpi_ue.
 * @param state foil state containing the converged inviscid solution
 * @param sens Sensitivity struct to populate with dcl_dalpha, dcm_dalpha, dcd_dalpha
 */
void AlphaSensitivities(FoilState& state, Sensitivity& sens);

/**
 * @brief Compute sensitivities of cl, cd, cm w.r.t. Mach number.
 * Inviscid only: differentiates KT compressibility correction analytically.
 * @param state foil state containing the converged inviscid solution
 * @param sens Sensitivity struct to populate with dcl_dMa, dcm_dMa, dcd_dMa
 */
void MachSensitivities(FoilState& state, Sensitivity& sens);

/**
 * @brief Solve the viscous adjoint transpose system R_V^T · λ = rhs.
 *
 * Requires that the viscous primal has been converged and that
 * StoreAdjointFactorization has populated state.glob.R_V and
 * state.glob.lu_factorization. If lu_factorization is null this call
 * rebuilds it (idempotent). Returns λ ∈ ℝ^{4·Nsys}.
 * @param state foil state with a valid adjoint LU factorization
 * @param rhs right-hand side vector of size 4*Nsys
 * @return adjoint solution vector λ of size 4*Nsys
 */
VectorXd SolveViscousAdjoint(FoilState& state, const VectorXd& rhs);

/**
 * @brief Compute sensitivities for the viscous problem using the adjoint method.
 * Must be called after SolveViscous() and CalcForce().
 *
 * Objectives: cl, cm, cd (Squire-Young total drag)
 * Design variables: foil.x [2 x N] panel node coordinates
 *
 * Adjoint formulation:
 *   R_V^T λ_J = ∂J/∂V
 *   dJ/dx_foil = ∂J/∂x_partial - λ_J^T · ∂R_V/∂x_foil
 * where the VJP λ_J^T · ∂R_V/∂x_foil = λ_bl^T · DRbl/∂x + ComputeRueVJP(λ_ue).
 * @param state foil state containing the converged viscous solution
 * @return Sensitivity struct populated with dcl_dx, dcm_dx, dcd_dx and scalar sensitivities
 */
Sensitivity ViscousSensitivities(FoilState& state);

/**
 * @brief Compute ∂R_V/∂alpha as a [4*Nsys] column vector.
 *
 * R_BL rows: zero (BL residual depends on U and xi only, not alpha directly).
 * R_ue rows: -∂ueinv/∂alpha, where
 *   foil i (0..N-1):  ueinv[i] = sgnue[i] * gamma[i],
 *                     gamma = gamref[:,0]*cos α + gamref[:,1]*sin α
 *   wake j=0 (continuity): same chain through gamma[N-1]
 *   wake j>0:  ueinv[N+j] = uewiref[j,0]*cos α + uewiref[j,1]*sin α
 * @param state foil state with a converged viscous solution
 * @return residual sensitivity vector ∂R_V/∂alpha of size 4*Nsys
 */
VectorXd DRVDAlpha(FoilState& state);

/**
 * @brief Compute ∂R_V/∂Ma as a [4*Nsys] column vector via central FD on
 * state.param.Minf with the BL closures' Mach-dependence integrated through
 * the residual. The Mach chain enters via GetUk (KT correction), GetMach2,
 * GetHss, GetRet (compressible density and Sutherland viscosity). A clean
 * analytic chain would require differentiating each closure w.r.t. Ma which
 * is invasive; the scalar FD here is well-conditioned (single direction)
 * and limited to ~1e-7 abs noise floor at h=1e-7.
 * @param state foil state with a converged viscous solution
 * @param h finite-difference step size (default 1e-7)
 * @return residual sensitivity vector ∂R_V/∂Ma of size 4*Nsys
 */
VectorXd DRVDMa(FoilState& state, double h = 1e-7);

/**
 * @brief Compute ∂R_V/∂Re as a [4*Nsys] column vector via central FD on
 * state.param.mu0. Re enters only through Ret in BL closures.
 * dJ/dRe = ∂J/∂Re - λ^T · ∂R_V/∂Re. Note ∂param.mu0/∂Re = -mu0/Re.
 * @param state foil state with a converged viscous solution
 * @param h finite-difference step size (default 1e-7)
 * @return residual sensitivity vector ∂R_V/∂Re of size 4*Nsys
 */
VectorXd DRVDRe(FoilState& state, double h = 1e-7);

}  // namespace scirocco
