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
#include "cppfoil/core/panel_influence.h"
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

void SolveViscous(FoilState& state) {
    // Top-level viscous solution

    if (state.param.verb >= 1) {
        std::cout << "SolveViscous: starting" << std::endl;
    }

    // Initialize thermodynamic parameters
    state.oper.viscous = true;
    InitThermo(state);

    // Build inviscid solution first
    BuildGamma(state, state.oper.alpha);
    StagPointFind(state);

    if (state.param.verb >= 1) {
        std::cout << "  Inviscid done, N=" << state.foil.N
                  << ", Istag=" << state.isol.Istag[0] << "," << state.isol.Istag[1]
                  << std::endl;
    }

    // Build wake
    BuildWake(state);

    // Extend xi to include wake points
    ExtendXiToWake(state);

    if (state.param.verb >= 1) {
        std::cout << "  Wake done, Nw=" << state.wake.N << std::endl;
    }

    // Identify surfaces
    IdentifySurfaces(state);

    // Set wake gap
    SetWakeGap(state);

    if (state.param.verb >= 1) {
        std::cout << "  Surfaces: lower=" << state.vsol.Is[0].size()
                  << ", upper=" << state.vsol.Is[1].size()
                  << ", wake=" << state.vsol.Is[2].size() << std::endl;
    }

    // Calculate ue sensitivities
    CalcUeM(state);

    // Initialize boundary layer
    if (state.oper.initbl) {
        InitBoundaryLayer(state);
        if (state.param.verb >= 1) {
            std::cout << "  BL init done, Nsys=" << state.glob.Nsys << std::endl;
        }
    }

    // Solve coupled system
    SolveCoupled(state);

    // Calculate forces
    CalcForce(state);
    GetDistributions(state);

    if (state.param.verb >= 1) {
        std::cout << "SolveViscous: done, cl=" << state.post.cl
                  << ", cd=" << state.post.cd << std::endl;
    }
}

void CheckJacobian(FoilState& state, int station) {
    // Finite-difference check of BL residual Jacobian at one station

    const auto& Is = state.vsol.Is[1];  // Upper surface
    if (station < 1 || station >= static_cast<int>(Is.size())) return;

    int ip = station - 1;
    int i = station;

    Param param = BuildParam(state, 1);  // Upper surface
    StationParam(state, param, Is[i]);

    VectorXd xi_surf(Is.size());
    for (size_t j = 0; j < Is.size(); ++j) {
        xi_surf(j) = state.isol.xi(Is[j]);
    }

    Matrix42d U0;
    U0.col(0) = state.glob.U.col(Is[ip]);
    U0.col(1) = state.glob.U.col(Is[i]);
    Vector2d x0(xi_surf(ip), xi_surf(i));
    Vector2d Aux0 = Vector2d::Zero();

    auto res0 = ResidualStation(param, x0, U0, Aux0);

    std::cout << "Jacobian check at station " << station << " (idx=" << Is[i] << "):" << std::endl;
    std::cout << "  U1 = (th=" << U0(0,0) << ", ds=" << U0(1,0) << ", sa=" << U0(2,0) << ", ue=" << U0(3,0) << ")" << std::endl;
    std::cout << "  U2 = (th=" << U0(0,1) << ", ds=" << U0(1,1) << ", sa=" << U0(2,1) << ", ue=" << U0(3,1) << ")" << std::endl;
    std::cout << "  x = (" << x0(0) << ", " << x0(1) << ")" << std::endl;
    std::cout << "  R = (" << res0.value(0) << ", " << res0.value(1) << ", " << res0.value(2) << ")" << std::endl;

    double eps = 1e-5;  // Larger step for better numerical accuracy
    for (int var = 0; var < 8; ++var) {
        Matrix42d Up = U0;
        int col = var / 4;
        int row = var % 4;
        double base = std::abs(Up(row, col));
        double h = eps * std::max(1e-6, base);  // Ensure h isn't too tiny
        Up(row, col) += h;

        auto resp = ResidualStation(param, x0, Up, Aux0);
        Vector3d dRdU_fd = (resp.value - res0.value) / h;
        Vector3d dRdU_an = res0.lin_U.col(var);
        double err = (dRdU_fd - dRdU_an).norm() / std::max(1.0, dRdU_an.norm());

        std::cout << "  dR/dU[" << var << "]: fd=(" << dRdU_fd(0) << "," << dRdU_fd(1) << "," << dRdU_fd(2)
                  << ") an=(" << dRdU_an(0) << "," << dRdU_an(1) << "," << dRdU_an(2)
                  << ") err=" << err << std::endl;
    }
}

void SolveCoupled(FoilState& state) {
    // Coupled Newton iteration

    int niglob = state.param.niglob;
    double rtol = state.param.rtol;

    // Debug: check Jacobian at a station
    if (state.param.verb >= 2) {
        CheckJacobian(state, 10);
    }

    for (int iter = 0; iter < niglob; ++iter) {
        // Build global system
        BuildGlobSys(state);

        // Check convergence
        double rnorm = state.glob.R.norm();
        if (state.param.verb >= 1) {
            std::cout << "  Iteration " << iter << ": |R| = " << rnorm << std::endl;

            // Debug: print max residual components
            if (iter < 3) {
                int Nsys = state.glob.Nsys;
                double maxRmom = 0, maxRshp = 0, maxRlag = 0;
                int imaxRmom = 0, imaxRshp = 0, imaxRlag = 0;
                for (int i = 0; i < Nsys; ++i) {
                    if (std::abs(state.glob.R(3*i)) > maxRmom) {
                        maxRmom = std::abs(state.glob.R(3*i));
                        imaxRmom = i;
                    }
                    if (std::abs(state.glob.R(3*i+1)) > maxRshp) {
                        maxRshp = std::abs(state.glob.R(3*i+1));
                        imaxRshp = i;
                    }
                    if (std::abs(state.glob.R(3*i+2)) > maxRlag) {
                        maxRlag = std::abs(state.glob.R(3*i+2));
                        imaxRlag = i;
                    }
                }
                std::cout << "    maxRmom=" << maxRmom << " @" << imaxRmom
                          << ", maxRshp=" << maxRshp << " @" << imaxRshp
                          << ", maxRlag=" << maxRlag << " @" << imaxRlag << std::endl;

                // Print state at max shape residual location
                std::cout << "    U@" << imaxRshp << ": th=" << state.glob.U(0, imaxRshp)
                          << ", ds=" << state.glob.U(1, imaxRshp)
                          << ", sa=" << state.glob.U(2, imaxRshp)
                          << ", ue=" << state.glob.U(3, imaxRshp)
                          << ", xi=" << state.isol.xi(imaxRshp) << std::endl;
                // Previous point
                if (imaxRshp > 0) {
                    int ip = imaxRshp - 1;
                    std::cout << "    U@" << ip << ": th=" << state.glob.U(0, ip)
                              << ", ds=" << state.glob.U(1, ip)
                              << ", sa=" << state.glob.U(2, ip)
                              << ", ue=" << state.glob.U(3, ip)
                              << ", xi=" << state.isol.xi(ip) << std::endl;
                }
            }
        }

        if (rnorm < rtol) {
            state.glob.conv = true;
            return;
        }

        if (!std::isfinite(rnorm)) {
            state.glob.conv = false;
            std::cout << "  Warning: NaN in residual, stopping" << std::endl;
            return;
        }

        // Solve and update
        SolveGlob(state);

        // Debug: print max update
        if (state.param.verb >= 1 && iter < 3) {
            double maxdth = state.glob.dU.row(0).cwiseAbs().maxCoeff();
            double maxdds = state.glob.dU.row(1).cwiseAbs().maxCoeff();
            double maxdsa = state.glob.dU.row(2).cwiseAbs().maxCoeff();
            double maxdue = state.glob.dU.row(3).cwiseAbs().maxCoeff();
            std::cout << "    max|dU|: dth=" << maxdth << ", dds=" << maxdds
                      << ", dsa=" << maxdsa << ", due=" << maxdue << std::endl;
        }

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

            // Get states
            Matrix42d Uin;
            Uin.col(0) = state.glob.U.col(Is[ip]);
            Uin.col(1) = state.glob.U.col(Is[i]);
            Vector2d xin(xi(ip), xi(i));
            Vector2d Aux_in(Aux(ip), Aux(i));

            // Compute residual - use transition residual at transition stations
            ResidualLin res;
            if (tran) {
                // Transition station: combined laminar + turbulent residual
                res = ResidualTransition(state, param, xin, Uin, Aux_in);
                // After transition, subsequent stations are turbulent
                param.turb = true;
            } else {
                // Update param for this station
                StationParam(state, param, Is[i]);
                // Regular residual
                res = ResidualStation(param, xin, Uin, Aux_in);
            }

            // Store in global system
            int Ig = 3 * Is[i];
            state.glob.R.segment<3>(Ig) += res.value;

            // Check for NaN in residual
            if (!res.value.allFinite()) {
                std::cout << "  BuildGlobSys NaN detected: si=" << si << " i=" << i
                          << " Is[i]=" << Is[i] << " tran=" << tran << " turb=" << param.turb << std::endl;
                std::cout << "    R=(" << res.value(0) << "," << res.value(1) << "," << res.value(2) << ")" << std::endl;
                std::cout << "    U[" << Is[ip] << "]: th=" << Uin(0,0) << " ds=" << Uin(1,0)
                          << " sa=" << Uin(2,0) << " ue=" << Uin(3,0) << std::endl;
                std::cout << "    U[" << Is[i] << "]: th=" << Uin(0,1) << " ds=" << Uin(1,1)
                          << " sa=" << Uin(2,1) << " ue=" << Uin(3,1) << std::endl;
                std::cout << "    xi=(" << xin(0) << "," << xin(1) << ") Aux=(" << Aux_in(0) << "," << Aux_in(1) << ")" << std::endl;
            }

            // Debug for node 0
            if (Is[i] == 0) {
                std::cout << "  BuildGlobSys node0: si=" << si << " i=" << i
                          << " Is[i-1]=" << Is[ip] << " Is[i]=" << Is[i]
                          << " R=(" << res.value(0) << "," << res.value(1) << "," << res.value(2) << ")"
                          << " turb=" << param.turb << " tran=" << tran << std::endl;
                std::cout << "    U[" << Is[ip] << "]: th=" << Uin(0,0) << " ds=" << Uin(1,0)
                          << " sa=" << Uin(2,0) << " ue=" << Uin(3,0) << std::endl;
                std::cout << "    U[" << Is[i] << "]: th=" << Uin(0,1) << " ds=" << Uin(1,1)
                          << " sa=" << Uin(2,1) << " ue=" << Uin(3,1) << std::endl;
                std::cout << "    xi=(" << xin(0) << "," << xin(1) << ")" << std::endl;
            }

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

    // Debug: check ue_m magnitude and NaN
    double max_mass = mass.cwiseAbs().maxCoeff();
    double max_ue_m = state.vsol.ue_m.cwiseAbs().maxCoeff();
    if (state.param.verb >= 1) {
        std::cout << "    Debug: max|mass|=" << max_mass << ", max|ue_m|=" << max_ue_m
                  << ", ue_m size=" << state.vsol.ue_m.rows() << "x" << state.vsol.ue_m.cols() << std::endl;

        // Check for NaN in BL residuals
        for (int i = 0; i < 3 * Nsys; ++i) {
            if (!std::isfinite(state.glob.R(i))) {
                std::cout << "    NaN in state.glob.R at i=" << i << " (node " << i/3 << ", comp " << i%3 << ")" << std::endl;
                // Print state at this node
                int node = i / 3;
                std::cout << "      U[" << node << "]: th=" << state.glob.U(0, node)
                          << ", ds=" << state.glob.U(1, node)
                          << ", sa=" << state.glob.U(2, node)
                          << ", ue=" << state.glob.U(3, node) << std::endl;
                break;
            }
        }
        // Check for NaN in ue residuals
        for (int i = 0; i < Nsys; ++i) {
            if (!std::isfinite(R(3 * Nsys + i))) {
                std::cout << "    NaN in ue residual at i=" << i
                          << " (ue=" << ue(i) << ", ueinv=" << ueinv(i)
                          << ", corr=" << ue_correction(i) << ")" << std::endl;
                break;
            }
        }
    }

    // Debug: check ue residual magnitude
    VectorXd R_ue = ue - (ueinv + ue_correction);
    double R_bl_norm = state.glob.R.norm();
    double R_ue_norm = R_ue.norm();
    if (state.param.verb >= 1) {
        // Find max ue residual
        int imax_ue = 0;
        double max_ue_res = 0;
        for (int i = 0; i < Nsys; ++i) {
            if (std::abs(R_ue(i)) > max_ue_res) {
                max_ue_res = std::abs(R_ue(i));
                imax_ue = i;
            }
        }
        std::cout << "    SolveGlob: |R_bl|=" << R_bl_norm << ", |R_ue|=" << R_ue_norm
                  << ", maxR_ue=" << max_ue_res << " @" << imax_ue
                  << " (ue=" << ue(imax_ue) << ", ueinv=" << ueinv(imax_ue)
                  << ", corr=" << ue_correction(imax_ue) << ")" << std::endl;
    }

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
    // Update state with adaptive under-relaxation (matching mfoil approach)

    int Nsys = state.glob.Nsys;
    double omega = 1.0;

    // Get max ctau for threshold calculations
    double ctmax = 0.0;
    std::vector<int> turbIdx;
    for (int i = 0; i < Nsys; ++i) {
        if (i < static_cast<int>(state.vsol.turb.size()) && state.vsol.turb[i]) {
            ctmax = std::max(ctmax, state.glob.U(2, i));
            turbIdx.push_back(i);
        }
    }

    // 1. Limit theta and delta* decreases to 50%
    for (int k = 0; k < 2; ++k) {
        double fmin = 0.0;
        for (int i = 0; i < Nsys; ++i) {
            double Uk = state.glob.U(k, i);
            double dUk = state.glob.dU(k, i);
            if (Uk > 1e-12) {
                fmin = std::min(fmin, dUk / Uk);
            }
        }
        if (fmin < -0.5) {
            double om = std::abs(0.5 / fmin);
            if (om < omega) {
                omega = om;
                if (state.param.verb >= 3) {
                    std::cout << "    th/ds decrease: omega = " << omega << std::endl;
                }
            }
        }
    }

    // 2. Prevent negative amplification/ctau (mfoil approach)
    for (int i = 0; i < Nsys; ++i) {
        if (i >= static_cast<int>(state.vsol.turb.size())) continue;
        double Uk = state.glob.U(2, i);
        double dUk = state.glob.dU(2, i);
        // Skip small values (don't limit small amp/ctau)
        if (!state.vsol.turb[i] && Uk < 0.2) continue;   // small amp
        if (state.vsol.turb[i] && Uk < 0.1 * ctmax) continue;  // small ctau
        if (Uk == 0.0 || dUk == 0.0) continue;
        // Check if update would make value negative
        if (Uk + dUk < 0.0) {
            double om = 0.8 * std::abs(Uk / dUk);
            if (om < omega) {
                omega = om;
                if (state.param.verb >= 3) {
                    std::cout << "    neg sa: omega = " << omega << std::endl;
                }
            }
        }
    }

    // 3. Limit amplification changes (laminar regions) to 2.0
    // Include all laminar points (matching mfoil)
    double dumax_amp = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        if (i < static_cast<int>(state.vsol.turb.size()) && !state.vsol.turb[i]) {
            dumax_amp = std::max(dumax_amp, std::abs(state.glob.dU(2, i)));
        }
    }
    if (dumax_amp > 0) {
        double om = 2.0 / dumax_amp;
        if (om < omega) {
            omega = om;
            if (state.param.verb >= 3) {
                std::cout << "    amp: omega = " << omega << std::endl;
            }
        }
    }

    // 4. Limit ctau changes (turbulent regions) to 0.05
    // Include all turbulent points (matching mfoil)
    double dumax_ct = 0.0;
    int imax_ct = -1;
    for (int i : turbIdx) {
        double du = std::abs(state.glob.dU(2, i));
        if (du > dumax_ct) {
            dumax_ct = du;
            imax_ct = i;
        }
    }
    if (state.param.verb >= 1) {
        std::cout << "    ctau limit check: turbIdx.size=" << turbIdx.size()
                  << " dumax_ct=" << dumax_ct << " imax_ct=" << imax_ct
                  << " turb[21]=" << (21 < static_cast<int>(state.vsol.turb.size()) ? state.vsol.turb[21] : -1)
                  << " dU[2,21]=" << state.glob.dU(2, 21)
                  << " U[2,21]=" << state.glob.U(2, 21) << std::endl;
    }
    if (dumax_ct > 0) {
        double om = 0.05 / dumax_ct;
        if (om < omega) {
            omega = om;
            if (state.param.verb >= 1) {
                std::cout << "    ctau limit applied: omega = " << omega
                          << " @" << imax_ct << " ct=" << state.glob.U(2, imax_ct)
                          << " dct=" << state.glob.dU(2, imax_ct) << std::endl;
            }
        }
    }

    // 5. Limit ue changes to 0.2 * Vinf
    double fmax_ue = 0.0;
    for (int i = 0; i < Nsys; ++i) {
        fmax_ue = std::max(fmax_ue, std::abs(state.glob.dU(3, i)) / state.oper.Vinf);
    }
    if (fmax_ue > 0.2) {
        double om = 0.2 / fmax_ue;
        if (om < omega) {
            omega = om;
            if (state.param.verb >= 3) {
                std::cout << "    ue: omega = " << omega << std::endl;
            }
        }
    }

    // 6. Limit alpha changes to 2 degrees
    if (std::abs(state.glob.dalpha) > 2.0) {
        omega = std::min(omega, 2.0 / std::abs(state.glob.dalpha));
    }

    if (state.param.verb >= 1) {
        std::cout << "    state update: omega = " << omega << std::endl;
    }

    // Store omega for convergence check
    double alpha_change = omega * state.glob.dalpha;

    // Apply update
    double U21_before = state.glob.U(2, 21);
    double dU21 = state.glob.dU(2, 21);
    state.glob.U += omega * state.glob.dU;
    state.oper.alpha += alpha_change;

    if (state.param.verb >= 1) {
        std::cout << "    Update check: U[2,21] " << U21_before << " + " << omega << " * " << dU21
                  << " = " << U21_before + omega * dU21 << " (actual=" << state.glob.U(2, 21) << ")" << std::endl;
    }

    // Fix bad Hk after update (using proper Hk calculation)
    for (int si = 0; si < 3; ++si) {
        double Hkmin = (si == 2) ? 1.00005 : 1.02;
        const auto& Is = state.vsol.Is[si];
        Param param = BuildParam(state, si);
        for (int idx : Is) {
            if (idx >= Nsys) continue;
            StationParam(state, param, idx);
            Vector4d Uj = state.glob.U.col(idx);
            auto [Hk, Hk_U] = GetHk(Uj, param);
            if (Hk < Hkmin) {
                state.glob.U(1, idx) += 2.0 * (Hkmin - Hk) * state.glob.U(1, idx);
            }
        }
    }

    // Fix negative ctau after update
    // Recompute ctmax after update
    ctmax = 0.0;
    for (int i : turbIdx) {
        ctmax = std::max(ctmax, state.glob.U(2, i));
    }
    for (int i : turbIdx) {
        if (state.glob.U(2, i) < 0) {
            state.glob.U(2, i) = 0.1 * ctmax;
        }
    }

    // Ensure positivity constraints
    for (int i = 0; i < Nsys; ++i) {
        state.glob.U(0, i) = std::max(state.glob.U(0, i), 1e-10);
        state.glob.U(1, i) = std::max(state.glob.U(1, i), 1e-10);
    }

    // Rebuild inviscid solution if alpha changed
    if (std::abs(alpha_change) > 1e-10) {
        RebuildIsol(state);
    }
}

void CalcUeM(FoilState& state) {
    // Calculate d(ue)/d(mass) sensitivities
    // This implements the transpiration BC coupling

    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    if (Nw == 0) {
        // No wake - just use identity
        state.vsol.ue_m.resize(Nsys, Nsys);
        state.vsol.ue_m.setIdentity();
        return;
    }

    // 1. Cgam = d(wake uei)/d(gamma) [Nw x N]
    // Use full InviscidVelocity to get linearization w.r.t. gamma, which includes
    // TE source and TE vortex contributions (matching mfoil's inviscid_velocity call)
    MatrixXd Cgam = MatrixXd::Zero(Nw, N);
    for (int i = 0; i < Nw; ++i) {
        Vector2d xi = state.wake.x.col(i);
        Vector2d ti = state.wake.t.col(i);

        MatrixXd V_G;
        InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi, true, &V_G);

        // Cgam(i,:) = V_G[0,:]*ti[0] + V_G[1,:]*ti[1] (tangent component)
        for (int j = 0; j < N; ++j) {
            Cgam(i, j) = V_G(0, j) * ti(0) + V_G(1, j) * ti(1);
        }
    }

    // 2. B = d(airfoil streamfunction)/d(source) [(N+1) x (N+Nw-2)]
    int Npanel = N + Nw - 2;  // N-1 airfoil panels + Nw-1 wake panels
    MatrixXd B = MatrixXd::Zero(N + 1, Npanel);

    for (int i = 0; i < N; ++i) {
        Vector2d xi = state.foil.x.col(i);

        // Constant source panels on airfoil
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);
            B(i, j) = PanelConstSourceStream(Xj, xi);
        }

        // Linear source panels on wake
        for (int j = 0; j < Nw - 1; ++j) {
            Matrix2d Xj;
            Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
            Xj.col(0) = state.wake.x.col(j);
            Xj.col(1) = Xm;

            // Ghost extension at last point
            Vector2d X_right;
            if (j == Nw - 2) {
                X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
            } else {
                X_right = state.wake.x.col(j + 1);
            }

            double a_left, b_left, a_right, b_right;

            // Left half panel
            PanelLinSourceStream(Xj, xi, a_left, b_left);
            if (j > 0) {
                B(i, N - 1 + j) += 0.5 * a_left + b_left;
                B(i, N - 1 + j - 1) += 0.5 * a_left;
            } else {
                B(i, N - 1 + j) += b_left;
            }

            // Right half panel
            Xj.col(0) = Xm;
            Xj.col(1) = X_right;
            PanelLinSourceStream(Xj, xi, a_right, b_right);
            B(i, N - 1 + j) += a_right + 0.5 * b_right;
            if (j < Nw - 2) {
                B(i, N - 1 + j + 1) += 0.5 * b_right;
            } else {
                B(i, N - 1 + j) += 0.5 * b_right;
            }
        }
    }

    // 3. Bp = -inv(AIC) * B [N x (N+Nw-2)]
    // Note: AIC is (N+1) x (N+1), but last row is Kutta condition
    MatrixXd Bp = -state.isol.AIC.lu().solve(B);
    Bp.conservativeResize(N, Npanel);  // Remove last row

    // Debug: check intermediate matrix magnitudes
    std::cout << "  CalcUeM debug: max|Cgam|=" << Cgam.cwiseAbs().maxCoeff()
              << ", max|B|=" << B.cwiseAbs().maxCoeff()
              << ", max|Bp|=" << Bp.cwiseAbs().maxCoeff() << std::endl;

    // 4. Csig = d(wake uei)/d(source) [Nw x (N+Nw-2)]
    MatrixXd Csig = MatrixXd::Zero(Nw, Npanel);

    // Debug: track contributions to Csig(1, 0)
    double csig10_airfoil = 0.0;
    double csig10_wake = 0.0;

    for (int i = 0; i < Nw; ++i) {
        Vector2d xi = state.wake.x.col(i);
        Vector2d ti = state.wake.t.col(i);

        // Constant source panels on airfoil
        int jstart = (i == 0) ? 1 : 0;
        int jend = (i == 0) ? N - 2 : N - 1;
        for (int j = jstart; j < jend; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);
            Vector2d v = PanelConstSourceVelocity(Xj, xi, &ti);
            Csig(i, j) = v(0);
            if (i == 1 && j == 0) {
                csig10_airfoil = v(0);
            }
        }

        // Piecewise linear source on wake
        for (int j = 0; j < Nw; ++j) {
            int il = std::max(j - 1, 0);
            int im = j;
            int ir = std::min(j + 1, Nw - 1);

            Vector2d Xl = 0.5 * (state.wake.x.col(il) + state.wake.x.col(im));
            Vector2d Xm = state.wake.x.col(im);
            Vector2d Xr = 0.5 * (state.wake.x.col(im) + state.wake.x.col(ir));

            if (j == Nw - 1) {
                Xr = 2.0 * Xm - Xl;  // Ghost extension
            }

            double d1 = (Xm - Xl).norm();
            double d2 = (Xr - Xm).norm();

            if (i == j) {
                // Self-influence (singular)
                if (j == 0) {
                    // First point: TE special case
                    double dl = (state.foil.x.col(1) - state.foil.x.col(0)).norm();
                    double du = (state.foil.x.col(N - 1) - state.foil.x.col(N - 2)).norm();
                    Csig(i, 0) += (0.5 / kPi) * (SafeLog(dl / d2) + 1.0);
                    Csig(i, N - 2) += (0.5 / kPi) * (SafeLog(du / d2) + 1.0);
                    Csig(i, N - 1) += -0.5 / kPi;
                } else if (j == Nw - 1) {
                    // Last point: no self effect
                } else {
                    double aa = (0.25 / kPi) * SafeLog(d1 / d2);
                    Csig(i, N - 1 + j - 1) += aa + 0.5 / kPi;
                    Csig(i, N - 1 + j) += aa - 0.5 / kPi;
                }
            } else {
                // Non-self influence
                if (j == 0) {
                    Matrix2d Xj;
                    Xj.col(0) = Xm;
                    Xj.col(1) = Xr;
                    Vector2d a, b;
                    PanelLinSourceVelocity(Xj, xi, a, b);
                    Csig(i, N - 1) += b.dot(ti);
                    Csig(i, 0) += a.dot(ti);
                    Csig(i, N - 2) += a.dot(ti);
                    if (i == 1) {
                        csig10_wake = a.dot(ti);
                    }
                } else if (j == Nw - 1) {
                    Matrix2d Xj;
                    Xj.col(0) = Xl;
                    Xj.col(1) = Xr;
                    Vector2d v = PanelConstSourceVelocity(Xj, xi, &ti);
                    Csig(i, N + Nw - 3) += v(0);
                } else {
                    Matrix2d Xj;
                    Vector2d a1, b1, a2, b2;

                    // Left half-panel
                    Xj.col(0) = Xl;
                    Xj.col(1) = Xm;
                    PanelLinSourceVelocity(Xj, xi, a1, b1);

                    // Right half-panel
                    Xj.col(0) = Xm;
                    Xj.col(1) = Xr;
                    PanelLinSourceVelocity(Xj, xi, a2, b2);

                    Csig(i, N - 1 + j - 1) += a1.dot(ti) + 0.5 * b1.dot(ti);
                    Csig(i, N - 1 + j) += 0.5 * a2.dot(ti) + b2.dot(ti);
                }
            }
        }
    }

    // Debug: Csig(1,0) breakdown
    std::cout << "  CalcUeM Csig(1,0) breakdown: airfoil=" << csig10_airfoil
              << ", wake=" << csig10_wake
              << ", total=" << Csig(1, 0) << std::endl;

    // 5. ue_sigma = d(unsigned ue)/d(source) [(N+Nw) x (N+Nw-2)]
    // Df = Bp for foil, Dw = Cgam*Bp + Csig for wake
    MatrixXd CgamBp = Cgam * Bp;
    MatrixXd Dw = CgamBp + Csig;
    Dw.row(0) = Bp.row(N - 1);  // First wake point same as TE

    // Debug: check components at row 132 = ue_sigma(132, 0)
    int row132_wake_idx = 132 - N;
    std::cout << "  CalcUeM: N=" << N << ", Nw=" << Nw
              << ", row 132 maps to wake idx " << row132_wake_idx << std::endl;
    if (row132_wake_idx >= 0 && row132_wake_idx < Nw) {
        std::cout << "  CalcUeM debug at wake(" << row132_wake_idx << ",0): CgamBp="
                  << CgamBp(row132_wake_idx, 0)
                  << ", Csig=" << Csig(row132_wake_idx, 0)
                  << ", Dw=" << Dw(row132_wake_idx, 0) << std::endl;
    } else if (132 < N) {
        std::cout << "  CalcUeM debug: row 132 is in Bp: Bp(132,0)=" << Bp(132, 0) << std::endl;
    }

    state.vsol.ue_sigma.resize(Nsys, Npanel);
    state.vsol.ue_sigma.topRows(N) = Bp;
    state.vsol.ue_sigma.bottomRows(Nw) = Dw;

    // 6. sigma_m = d(source)/d(mass) [(N+Nw-2) x (N+Nw)]
    std::vector<Triplet> triplets;
    for (int i = 0; i < N - 1; ++i) {
        double ds = state.foil.s(i + 1) - state.foil.s(i);
        if (ds > 1e-12) {
            triplets.emplace_back(i, i, -state.isol.sgnue(i) / ds);
            triplets.emplace_back(i, i + 1, state.isol.sgnue(i + 1) / ds);
        }
    }
    for (int i = 0; i < Nw - 1; ++i) {
        double ds = state.wake.s(i + 1) - state.wake.s(i);
        if (ds > 1e-12) {
            triplets.emplace_back(N - 1 + i, N + i, -1.0 / ds);
            triplets.emplace_back(N - 1 + i, N + i + 1, 1.0 / ds);
        }
    }
    state.vsol.sigma_m.resize(Npanel, Nsys);
    state.vsol.sigma_m.setFromTriplets(triplets.begin(), triplets.end());

    // 7. ue_m = sgnue * ue_sigma * sigma_m
    VectorXd sgue(Nsys);
    sgue.head(N) = state.isol.sgnue;
    sgue.tail(Nw).setOnes();

    state.vsol.ue_m = sgue.asDiagonal() * state.vsol.ue_sigma * state.vsol.sigma_m;

    // Debug: check matrix magnitudes
    double max_ue_sigma = state.vsol.ue_sigma.cwiseAbs().maxCoeff();
    double max_sigma_m = 0.0;
    for (int k = 0; k < state.vsol.sigma_m.outerSize(); ++k) {
        for (SparseMatrixd::InnerIterator it(state.vsol.sigma_m, k); it; ++it) {
            max_sigma_m = std::max(max_sigma_m, std::abs(it.value()));
        }
    }
    std::cout << "  CalcUeM: max|ue_sigma|=" << max_ue_sigma
              << ", max|sigma_m|=" << max_sigma_m << std::endl;
}

void RebuildUeM(FoilState& state) {
    // Rebuild sensitivities after stagnation point move
    CalcUeM(state);
}

void RebuildIsol(FoilState& state) {
    // Rebuilds inviscid solution after angle of attack change
    // Recompute gamma from reference solutions

    double alpha = state.oper.alpha;
    double cosa = Cosd(alpha);
    double sina = Sind(alpha);

    // Rebuild gamma on airfoil: gam = gamref[:,0]*cos(alpha) + gamref[:,1]*sin(alpha)
    int N = state.foil.N;
    for (int i = 0; i < N; ++i) {
        state.isol.gam(i) = state.isol.gamref(i, 0) * cosa +
                            state.isol.gamref(i, 1) * sina;
    }

    // For viscous flow, may need to rebuild wake and ue_m
    if (state.oper.viscous && state.oper.redowake) {
        BuildWake(state);
        IdentifySurfaces(state);
        CalcUeM(state);
    }
}

}  // namespace cppfoil
