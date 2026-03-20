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

// Compute sensitivities of cl, cd, cm w.r.t. Mach number.
// Inviscid: differentiates KT compressibility correction analytically.
// Viscous: also computes dR_BL/dMa and uses adjoint solve.
void MachSensitivities(FoilState& state, Sensitivity& sens);

// Compute sensitivities of cl, cd, cm w.r.t. Reynolds number.
// Inviscid: zero (Re doesn't affect inviscid solution).
// Viscous: Re enters through mu0; reuses Mach companion functions.
void ReSensitivities(FoilState& state, Sensitivity& sens);

// Compute dxi/dx_foil: sensitivity of arc-length coordinate xi to panel coords
// Returns sparse matrix [Nsys x 2*N]
SparseMatrixd ComputeDxiDxFoil(const FoilState& state);

// Compute dR_ue/dx_foil: sensitivity of 4th equation to panel coords
// Returns dense matrix [Nsys x 2*N]
MatrixXd ComputeDRueDxFoil(FoilState& state);

// Compute VJP of wake geometry through R_ue equation.
// Returns lambda_ue^T * (∂R_ue/∂wake.x · dwake.x/dx_foil + ∂R_ue/∂wake.t · dwake.t/dx_foil)
// This captures the missing chain rule terms where wake positions and tangents
// depend on foil coordinates through streamline tracing.
// INPUT:
//   state : FoilState with converged viscous solution
//   lambda_ue : Adjoint vector for the ue equation [Nsys]
// OUTPUT:
//   VJP vector [2*N] in foil.x space
VectorXd ComputeWakeGeomVJP(FoilState& state, const VectorXd& lambda_ue);

// Compute VJP of stagnation point location through xi coordinate.
// Returns lambda_bl^T * R_x * (∂xi/∂sstag · dsstag/dx_foil)
// This captures the sensitivity of BL arc-length coordinate xi to the
// stagnation point location, which depends on gamma which depends on foil.x.
// INPUT:
//   state : FoilState with converged viscous solution
//   lambda_bl : Adjoint vector for the BL equations [3*Nsys]
// OUTPUT:
//   VJP vector [2*N] in foil.x space
VectorXd ComputeXiVJP(FoilState& state, const VectorXd& lambda_bl);

// Compute VJP of wake gap (wgap) through BL residual.
// Returns lambda_bl^T * dR_BL/d(wgap) * d(wgap)/d(foil.x)
// The BL residual subtracts wgap from delta* at wake stations;
// wgap depends on foil.x through hTE and dtdx at TE nodes {0,1,N-2,N-1}.
// INPUT:
//   state : FoilState with converged viscous solution
//   lambda_bl : Adjoint vector for the BL equations [3*Nsys]
// OUTPUT:
//   VJP vector [2*N] in foil.x space
VectorXd ComputeWgapVJP(FoilState& state, const VectorXd& lambda_bl);

// Compute d(tcp)/d(x[k]_dim) and d(tdp)/d(x[k]_dim) analytically
void GetTEInfoDerivPublic(const MatrixXd& x, int k, int dim,
                           double& dtcp, double& dtdp);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_ADJOINT_H
