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
#ifndef CPPFOIL_CORE_ADJOINT_H
#define CPPFOIL_CORE_ADJOINT_H

#include "cppfoil/data/sensitivity.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

// Compute sensitivities for the inviscid problem using the adjoint method.
// Must be called after SolveInviscid() and CalcForce().
//
// Objectives: cl, cm, cdpi
// Design variables: foil.x [2 x N] panel node coordinates
//
// The adjoint equation is: AIC^T * lambda = (dJ/dgamma)^T
// Total sensitivity: dJ/dx = dJ/dx_partial + lambda^T * d(A*gamma - rhs)/dx
Sensitivity InviscidSensitivities(FoilState& state);

// Compute force partial derivatives w.r.t. panel coordinates (holding gamma fixed)
// Populates: dcl_dx_partial, dcm_dx_partial, dcdpi_dx_partial in state.post
void CalcForceLin(FoilState& state);

// Compute sensitivities of cl, cd, cm w.r.t. angle of attack (alpha).
// Inviscid only: uses stored gamref and cl_ue/cm_ue/cdpi_ue.
void AlphaSensitivities(FoilState& state, Sensitivity& sens);

// Compute sensitivities of cl, cd, cm w.r.t. Mach number.
// Inviscid only: differentiates KT compressibility correction analytically.
void MachSensitivities(FoilState& state, Sensitivity& sens);

// Solve the viscous adjoint transpose system R_V^T · λ = rhs.
//
// Requires that the viscous primal has been converged and that
// StoreAdjointFactorization has populated state.glob.R_V and
// state.glob.lu_factorization. If lu_factorization is null this call
// rebuilds it (idempotent). Returns λ ∈ ℝ^{4·Nsys}.
VectorXd SolveViscousAdjoint(FoilState& state, const VectorXd& rhs);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_ADJOINT_H
