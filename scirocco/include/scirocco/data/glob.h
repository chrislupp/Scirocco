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

#include <memory>

#include "scirocco/config.h"

namespace scirocco {

/* Global system variables */
struct Glob {
    int Nsys = 0;        /**< Number of equations and states */
    MatrixXd U;          /**< Primary states [th, ds, sa, ue] [4 x Nsys] */
    MatrixXd dU;         /**< Primary state update [4 x Nsys] */
    double dalpha = 0.0; /**< Angle of attack update */
    bool conv = false;   /**< Converged flag */
    VectorXd R;          /**< Residuals [3*Nsys] */
    SparseMatrixd R_U;   /**< Residual Jacobian w.r.t. primary states */
    SparseMatrixd R_x;   /**< Residual Jacobian w.r.t. xi [3*Nsys x Nsys] */
    SparseMatrixd R_V;   /**< Global Jacobian [4*Nsys x 4*Nsys] */
    /* ∂R_BL/∂Re·Re for the analytic Re chain (accumulated in BuildGlobSys
       from each station's ResidualLin.lin_Ret). Re enters R only through
       Ret = ρ·ue·θ/μ_0 with ∂Ret/∂Re = Ret/Re, so ∂R/∂Re = dR_dRe_times_Re /
       Re. Size 3*Nsys (BL rows); ue-equation rows (3*Nsys..4*Nsys) are
       Re-indep. */
    VectorXd dR_dRe_times_Re;
    /* ∂R_BL/∂Ma·Ma for the analytic Ma chain (accumulated from each station's
       ResidualLin.lin_Ma). Ma enters R via KTb/KTl/H0/cps/mu0/rho0 and Minf,
       and pre-multiplying by Ma here keeps the BuildGlobSys aggregation
       uniform: ∂R/∂Ma = dR_dMa_times_Ma / Ma. Size 3*Nsys (BL rows). */
    VectorXd dR_dMa_times_Ma;
    bool realloc = false; /**< If true, system Jacobians will be reallocated */

    /* Stored SparseLU factorization of R_V for adjoint solves */
    std::shared_ptr<Eigen::SparseLU<SparseMatrixd>> lu_factorization;

    /* PTC state — current pseudo-time step and the previous-iter residual
       norm, used by the SER tau update in SolveCoupled. */
    double tau = 0.0;
    double R_norm_prev = 0.0;
};

}  // namespace scirocco
