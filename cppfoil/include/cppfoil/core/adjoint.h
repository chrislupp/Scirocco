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

// Compute sensitivities for the viscous problem using the adjoint method.
// Must be called after SolveViscous() and CalcForce().
// Requires stored LU factorization from StoreAdjointFactorization().
//
// Objectives: cl, cd, cm
// Design variables: foil.x [2 x N] panel node coordinates
Sensitivity ViscousSensitivities(FoilState& state);

// Compute force partial derivatives w.r.t. panel coordinates (holding gamma fixed)
// Populates: dcl_dx_partial, dcm_dx_partial, dcdpi_dx_partial in state.post
void CalcForceLin(FoilState& state);

// Compute full dcd/dV [4*Nsys] linearization for Squire-Young + skin friction
void CalcCdLin(FoilState& state);

// Compute sensitivities of cl, cd, cm w.r.t. angle of attack (alpha).
// Must be called after CalcForce() and CalcForceLin().
// For inviscid: uses stored gamref and cl_ue/cm_ue/cdpi_ue.
// For viscous: also uses stored LU factorization for adjoint solve.
void AlphaSensitivities(FoilState& state, Sensitivity& sens);

// Compute dxi/dx_foil: sensitivity of arc-length coordinate xi to panel coords
// Returns sparse matrix [Nsys x 2*N]
SparseMatrixd ComputeDxiDxFoil(const FoilState& state);

// Compute dR_ue/dx_foil: sensitivity of 4th equation to panel coords
// Returns dense matrix [Nsys x 2*N]
MatrixXd ComputeDRueDxFoil(FoilState& state);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_ADJOINT_H
