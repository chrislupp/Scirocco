// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/coupled_solver.h"

#include <iostream>
#include <cmath>
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/core/boundary_layer.h"
#include "cppfoil/core/transition.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

void SolveViscous(FoilState& state) {
    // Top-level viscous solution

    // Build inviscid solution first
    BuildGamma(state, state.oper.alpha);
    StagPointFind(state);

    // Build wake
    BuildWake(state);

    // Identify surfaces
    IdentifySurfaces(state);

    // Set wake gap
    SetWakeGap(state);

    // Calculate ue sensitivities
    CalcUeM(state);

    // Initialize boundary layer
    if (state.oper.initbl) {
        InitBoundaryLayer(state);
    }

    // Solve coupled system
    SolveCoupled(state);

    // Calculate forces
    CalcForce(state);
    GetDistributions(state);
}

void SolveCoupled(FoilState& state) {
    // Coupled Newton iteration

    int niglob = state.param.niglob;
    double rtol = state.param.rtol;

    for (int iter = 0; iter < niglob; ++iter) {
        // Build global system
        BuildGlobSys(state);

        // Check convergence
        double rnorm = state.glob.R.norm();
        if (state.param.verb >= 1) {
            std::cout << "  Iteration " << iter << ": |R| = " << rnorm << std::endl;
        }

        if (rnorm < rtol) {
            state.glob.conv = true;
            return;
        }

        // Solve and update
        SolveGlob(state);
        UpdateState(state);

        // Update stagnation point
        StagPointMove(state);

        // Update transition
        UpdateTransition(state);
    }

    state.glob.conv = false;
    std::cout << "  Warning: viscous solver did not converge" << std::endl;
}

void BuildGlobSys(FoilState& state) {
    // Build global residual and Jacobian

    int Nsys = state.glob.Nsys;

    // Allocate residual and Jacobians if needed
    bool alloc = state.glob.realloc || state.glob.R.size() != 3 * Nsys;
    if (alloc) {
        state.glob.R.resize(3 * Nsys);
        state.glob.R_U.resize(3 * Nsys, 4 * Nsys);
        state.glob.R_x.resize(3 * Nsys, Nsys);
        state.glob.realloc = false;
    }

    state.glob.R.setZero();

    // Build using triplets for sparse matrix
    std::vector<Triplet> triplets_U, triplets_x;
    triplets_U.reserve(Nsys * 24);  // Rough estimate
    triplets_x.reserve(Nsys * 6);

    // Loop over surfaces (lower, upper, wake)
    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        if (Is.empty()) continue;

        int N = static_cast<int>(Is.size());

        // Get parameter structure for this surface
        Param param = BuildParam(state, si);

        // Get xi values for this surface
        VectorXd xi(N);
        for (int idx = 0; idx < N; ++idx) {
            xi(idx) = state.isol.xi(Is[idx]);
        }

        // Get auxiliary data (wake gap)
        VectorXd Aux = VectorXd::Zero(N);
        if (si == 2 && state.vsol.wgap.size() > 0) {
            for (int idx = 0; idx < N && idx < static_cast<int>(state.vsol.wgap.size()); ++idx) {
                Aux(idx) = state.vsol.wgap(idx);
            }
        }

        // Determine starting index (skip if first xi is tiny - stagnation)
        int i0 = (si < 2 && xi(0) < 1e-8 * xi(N - 1)) ? 1 : 0;

        // First point system
        Vector3d R1;
        Matrix38d R1_U;
        Matrix32d R1_x;  // 3x2 for x linearization
        std::vector<int> J;

        if (si < 2) {
            // Stagnation state for lower/upper surfaces
            int ip0 = i0, ip1 = i0 + 1;
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip0]);
            Uin.col(1) = state.glob.U.col(Is[ip1]);
            Vector2d xin(xi(ip0), xi(ip1));

            auto stag = StagnationState(Uin, xin);

            // Similarity station residual
            param.turb = false;
            param.simi = true;
            Matrix42d Ust_pair;
            Ust_pair.col(0) = stag.Ust;
            Ust_pair.col(1) = stag.Ust;
            Vector2d xst_pair(stag.xst, stag.xst);
            Vector2d Aux_pair(Aux(ip0), Aux(ip0));

            auto res = ResidualStation(param, xst_pair, Ust_pair, Aux_pair);
            param.simi = false;

            // Chain rule: R1_U = R1_Ut * [Ust_U]
            // R1_Ust is 3x4: sum of both 3x4 blocks from res.lin_U
            Matrix34d R1_Ust = res.lin_U.block<3, 4>(0, 0) + res.lin_U.block<3, 4>(0, 4);
            R1 = res.value;
            R1_U.setZero();
            R1_U.block<3, 4>(0, 0) = R1_Ust * stag.Ust_U.block<4, 4>(0, 0);
            R1_U.block<3, 4>(0, 4) = R1_Ust * stag.Ust_U.block<4, 4>(0, 4);
            R1_x = R1_Ust * stag.Ust_x;

            J = {Is[ip0], Is[ip1]};

            // Handle i0=1 case (stagnation on node)
            if (i0 == 1) {
                int Ig = 3 * Is[0];
                state.glob.R.segment<3>(Ig) = state.glob.U.col(Is[0]).head<3>() - stag.Ust.head<3>();
                // d/dU[Is[0]]
                for (int k = 0; k < 3; ++k) {
                    triplets_U.emplace_back(Ig + k, 4 * Is[0] + k, 1.0);
                }
                // d/dU[J]
                for (int k = 0; k < 3; ++k) {
                    for (int c = 0; c < 4; ++c) {
                        if (std::abs(stag.Ust_U(k, c)) > 1e-15)
                            triplets_U.emplace_back(Ig + k, 4 * J[0] + c, -stag.Ust_U(k, c));
                        if (std::abs(stag.Ust_U(k, c + 4)) > 1e-15)
                            triplets_U.emplace_back(Ig + k, 4 * J[1] + c, -stag.Ust_U(k, c + 4));
                    }
                }
            }
        } else {
            // Wake initialization
            auto wake = WakeSys(state, param);
            R1 = wake.R;
            R1_U.setZero();
            // Map 3x12 to 3x8 format (we only store the needed blocks)
            R1_U.block<3, 4>(0, 0) = wake.R_U.block<3, 4>(0, 0);  // Lower
            R1_U.block<3, 4>(0, 4) = wake.R_U.block<3, 4>(0, 4);  // Upper
            // Wake contribution (indices 8-11 in R_U) stored separately
            R1_x.setZero();
            J = {wake.J[0], wake.J[1]};
            param.turb = true;
            param.wake = true;

            // Store first point system
            int Ig = 3 * Is[i0];
            state.glob.R.segment<3>(Ig) = R1;
            // Add wake state contribution
            for (int k = 0; k < 3; ++k) {
                for (int c = 0; c < 4; ++c) {
                    if (std::abs(wake.R_U(k, c)) > 1e-15)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[0] + c, wake.R_U(k, c));
                    if (std::abs(wake.R_U(k, c + 4)) > 1e-15)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[1] + c, wake.R_U(k, c + 4));
                    if (std::abs(wake.R_U(k, c + 8)) > 1e-15)
                        triplets_U.emplace_back(Ig + k, 4 * wake.J[2] + c, wake.R_U(k, c + 8));
                }
            }
        }

        // Store first point system for non-wake surfaces
        if (si < 2) {
            int Ig = 3 * Is[i0];
            state.glob.R.segment<3>(Ig) = R1;
            for (size_t j = 0; j < J.size(); ++j) {
                for (int k = 0; k < 3; ++k) {
                    for (int c = 0; c < 4; ++c) {
                        double val = R1_U(k, 4 * j + c);
                        if (std::abs(val) > 1e-15)
                            triplets_U.emplace_back(Ig + k, 4 * J[j] + c, val);
                    }
                }
            }
            // R1_x linearization w.r.t. xi is handled through stagnation point sens
            for (size_t j = 0; j < J.size(); ++j) {
                for (int k = 0; k < 3; ++k) {
                    double val = R1_x(k, j);
                    if (std::abs(val) > 1e-15)
                        triplets_x.emplace_back(Ig + k, J[j], val);
                }
            }
        }

        // March over rest of points
        for (int i = i0 + 1; i < N; ++i) {
            int ip = i - 1;
            bool tran = (Is[ip] < static_cast<int>(state.vsol.turb.size()) &&
                        Is[i] < static_cast<int>(state.vsol.turb.size())) &&
                        (state.vsol.turb[Is[ip]] != state.vsol.turb[Is[i]]);

            // Update param for this station
            StationParam(state, param, Is[i]);

            // Get states
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip]);
            Uin.col(1) = state.glob.U.col(Is[i]);
            Vector2d xin(xi(ip), xi(i));
            Vector2d Aux_in(Aux(ip), Aux(i));

            // Compute residual
            auto res = ResidualStation(param, xin, Uin, Aux_in);

            // Store in global system
            int Ig = 3 * Is[i];
            state.glob.R.segment<3>(Ig) += res.value;

            // Store Jacobian entries
            for (int k = 0; k < 3; ++k) {
                for (int c = 0; c < 4; ++c) {
                    double val1 = res.lin_U(k, c);
                    double val2 = res.lin_U(k, c + 4);
                    if (std::abs(val1) > 1e-15)
                        triplets_U.emplace_back(Ig + k, 4 * Is[ip] + c, val1);
                    if (std::abs(val2) > 1e-15)
                        triplets_U.emplace_back(Ig + k, 4 * Is[i] + c, val2);
                }
                double val_x1 = res.lin_x(k, 0);
                double val_x2 = res.lin_x(k, 1);
                if (std::abs(val_x1) > 1e-15)
                    triplets_x.emplace_back(Ig + k, Is[ip], val_x1);
                if (std::abs(val_x2) > 1e-15)
                    triplets_x.emplace_back(Ig + k, Is[i], val_x2);
            }

            // If transition occurred, subsequent stations are turbulent
            if (tran) {
                param.turb = true;
            }
        }
    }

    // Build sparse matrices from triplets
    state.glob.R_U.setFromTriplets(triplets_U.begin(), triplets_U.end());
    state.glob.R_x.setFromTriplets(triplets_x.begin(), triplets_x.end());

    // Include effects of R_x into R_U: R_ue += R_x * x_st * st_ue
    // x_st = -sgnue (sensitivity of xi to stagnation point location)
    VectorXd x_st(Nsys);
    int N = state.foil.N;
    for (int i = 0; i < N; ++i) {
        x_st(i) = -state.isol.sgnue(i);
    }
    for (int i = N; i < Nsys; ++i) {
        x_st(i) = -1.0;  // Wake has same sens as upper surface
    }

    VectorXd R_st = state.glob.R_x * x_st;
    const auto& Istag = state.isol.Istag;
    const auto& st_ue = state.isol.sstag_ue;

    // Add contributions at stagnation point ue indices
    for (int row = 0; row < 3 * Nsys; ++row) {
        if (std::abs(R_st(row)) > 1e-15) {
            // ue index for stagnation points
            int ue_idx0 = 4 * Istag[0] + 3;
            int ue_idx1 = 4 * Istag[1] + 3;
            state.glob.R_U.coeffRef(row, ue_idx0) += R_st(row) * st_ue(0);
            state.glob.R_U.coeffRef(row, ue_idx1) += R_st(row) * st_ue(1);
        }
    }
}

void SolveGlob(FoilState& state) {
    // Solve the global linear system using augmented system
    // Variables = [th, ds, sa, ue] at each node
    // Fourth equation: ue = ueinv + ue_m * (ds * ue)

    int Nsys = state.glob.Nsys;
    bool docl = state.oper.givencl;

    // Get edge velocity and displacement thickness
    VectorXd ue = state.glob.U.row(3).transpose();
    VectorXd ds = state.glob.U.row(1).transpose();
    double uemax = ue.cwiseAbs().maxCoeff();
    for (int i = 0; i < Nsys; ++i) {
        ue(i) = std::max(ue(i), 1e-10 * uemax);
    }

    // Get inviscid edge velocity
    VectorXd ueinv = GetUeInv(state);

    // Build augmented system
    int NN = 4 * Nsys + (docl ? 1 : 0);
    bool alloc = state.glob.realloc ||
                 state.glob.R_V.rows() != NN || state.glob.R_V.cols() != NN;

    if (alloc) {
        state.glob.R_V.resize(NN, NN);
    }

    // Build residual vector
    VectorXd R(NN);
    R.head(3 * Nsys) = state.glob.R;

    // Fourth equation residual: ue - (ueinv + ue_m * (ds .* ue))
    VectorXd mass = ds.cwiseProduct(ue);
    VectorXd ue_correction = state.vsol.ue_m * mass;
    R.segment(3 * Nsys, Nsys) = ue - (ueinv + ue_correction);

    // Build Jacobian using triplets
    std::vector<Triplet> triplets;
    triplets.reserve(state.glob.R_U.nonZeros() + 4 * Nsys);

    // Copy existing R_U entries
    for (int k = 0; k < state.glob.R_U.outerSize(); ++k) {
        for (SparseMatrixd::InnerIterator it(state.glob.R_U, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), it.value());
        }
    }

    // Fourth equation Jacobian
    // d(R_ue)/d(ue) = I - ue_m * diag(ds)
    // d(R_ue)/d(ds) = -ue_m * diag(ue)
    for (int i = 0; i < Nsys; ++i) {
        int row = 3 * Nsys + i;
        int ue_idx = 4 * i + 3;  // ue index in state vector

        // Identity part for ue
        triplets.emplace_back(row, ue_idx, 1.0);

        // ue_m contributions
        for (int j = 0; j < Nsys; ++j) {
            double ue_m_ij = state.vsol.ue_m(i, j);
            if (std::abs(ue_m_ij) > 1e-15) {
                // d/d(ue_j): -ue_m(i,j) * ds(j)
                triplets.emplace_back(row, 4 * j + 3, -ue_m_ij * ds(j));
                // d/d(ds_j): -ue_m(i,j) * ue(j)
                triplets.emplace_back(row, 4 * j + 1, -ue_m_ij * ue(j));
            }
        }
    }

    // Cl-alpha constraint if requested
    if (docl) {
        double cl = state.post.cl;
        double cltgt = state.oper.cltgt;
        R(NN - 1) = cl - cltgt;

        // d(cl)/d(alpha)
        triplets.emplace_back(NN - 1, NN - 1, state.post.cl_alpha);

        // d(cl)/d(ue) at airfoil nodes
        int N = state.foil.N;
        for (int i = 0; i < N; ++i) {
            if (i < static_cast<int>(state.post.cl_ue.size())) {
                triplets.emplace_back(NN - 1, 4 * i + 3, state.post.cl_ue(i));
            }
        }

        // d(ue residual)/d(alpha)
        // uinv = uinvref * [cos(alpha), sin(alpha)]
        // d(uinv)/d(alpha) = uinvref * [-sin(alpha), cos(alpha)] * pi/180
        double alpha = state.oper.alpha;
        for (int i = 0; i < Nsys; ++i) {
            double duinv_dalpha = 0.0;
            if (i < state.foil.N) {
                duinv_dalpha = -state.isol.gamref(i, 0) * Sind(alpha) +
                               state.isol.gamref(i, 1) * Cosd(alpha);
            } else if (i - state.foil.N < static_cast<int>(state.isol.uewiref.rows())) {
                int iw = i - state.foil.N;
                duinv_dalpha = -state.isol.uewiref(iw, 0) * Sind(alpha) +
                               state.isol.uewiref(iw, 1) * Cosd(alpha);
            }
            duinv_dalpha *= kDegToRad;
            triplets.emplace_back(3 * Nsys + i, NN - 1, -duinv_dalpha);
        }
    }

    // Build sparse matrix
    state.glob.R_V.setFromTriplets(triplets.begin(), triplets.end());

    // Solve sparse system
    Eigen::SparseLU<SparseMatrixd> solver;
    solver.compute(state.glob.R_V);
    if (solver.info() != Eigen::Success) {
        std::cerr << "  Warning: SparseLU decomposition failed" << std::endl;
        state.glob.dU.setZero();
        state.glob.dalpha = 0.0;
        return;
    }

    VectorXd dV = -solver.solve(R);
    if (solver.info() != Eigen::Success) {
        std::cerr << "  Warning: SparseLU solve failed" << std::endl;
        state.glob.dU.setZero();
        state.glob.dalpha = 0.0;
        return;
    }

    // Store dU (reshaped from 4*Nsys vector to 4xNsys matrix)
    if (state.glob.dU.cols() != Nsys) {
        state.glob.dU.resize(4, Nsys);
    }
    for (int i = 0; i < Nsys; ++i) {
        state.glob.dU.col(i) = dV.segment<4>(4 * i);
    }

    if (docl) {
        state.glob.dalpha = dV(NN - 1);
    } else {
        state.glob.dalpha = 0.0;
    }
}

void UpdateState(FoilState& state) {
    // Update state with under-relaxation

    double omega = 0.5;  // Under-relaxation factor

    state.glob.U += omega * state.glob.dU;
    state.oper.alpha += omega * state.glob.dalpha;

    // Enforce physical constraints
    // - theta > 0
    // - delta* > 0
    // - bounded shape parameters
    for (int i = 0; i < state.glob.Nsys; ++i) {
        state.glob.U(0, i) = std::max(state.glob.U(0, i), 1e-10);
        state.glob.U(1, i) = std::max(state.glob.U(1, i), 1e-10);
    }
}

void CalcUeM(FoilState& state) {
    // Calculate d(ue)/d(mass) sensitivities
    // This is a placeholder - full implementation would compute
    // influence matrices for transpiration BC

    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    state.vsol.ue_m.resize(Nsys, Nsys);
    state.vsol.ue_m.setIdentity();
}

void RebuildUeM(FoilState& state) {
    // Rebuild sensitivities after stagnation point move
    CalcUeM(state);
}

}  // namespace cppfoil
