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

#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/data/state.h"

using namespace cppfoil;
using namespace Eigen;

// ---------------------------------------------------------------------------
// Tolerances
// ---------------------------------------------------------------------------
// Total VJP check: C++ sum of all VJP chains vs lambda^T @ dR/dx_CS
constexpr double TOL_VJP = 1e-4;
// Total sensitivity: end-to-end check
constexpr double TOL_TOTAL = 1e-4;

// ---------------------------------------------------------------------------
// HDF5 helpers
// ---------------------------------------------------------------------------
template <typename T>
static T readScalar(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    T value;
    dataset.read(&value, dataset.getDataType());
    return value;
}

static VectorXd readVector(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(buffer.data(), dims[0]);
}

static MatrixXd readMatrix(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd result(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j)
            result(i, j) = buffer[i * dims[1] + j];
    return result;
}

static std::vector<double> read3D(H5::H5File& file, const std::string& path,
                                   hsize_t& d0, hsize_t& d1, hsize_t& d2) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[3];
    dataspace.getSimpleExtentDims(dims);
    d0 = dims[0]; d1 = dims[1]; d2 = dims[2];
    std::vector<double> buffer(d0 * d1 * d2);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    return buffer;
}

/// Row-major 3D indexing: [i][j][k] in shape [d0 x d1 x d2]
static inline size_t idx3(int i, int j, int k, int d1, int d2) {
    return static_cast<size_t>(i) * d1 * d2 + j * d2 + k;
}

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------
class ComprehensiveReferenceTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool file_loaded_;

    static void SetUpTestSuite() {
        try {
            file_ = new H5::H5File(COMPREHENSIVE_REFERENCE_FILE, H5F_ACC_RDONLY);
            file_loaded_ = true;
            std::cout << "Loaded comprehensive reference: "
                      << COMPREHENSIVE_REFERENCE_FILE << std::endl;
        } catch (H5::Exception& e) {
            std::cerr << "Failed to open HDF5 file: " << e.getCDetailMsg()
                      << std::endl;
            file_loaded_ = false;
        }
    }

    static void TearDownTestSuite() {
        if (file_) { file_->close(); delete file_; file_ = nullptr; }
    }

    void SetUp() override {
        if (!file_loaded_)
            GTEST_SKIP() << "Comprehensive reference file not loaded";
    }

    H5::H5File& getFile() { return *file_; }
};

H5::H5File* ComprehensiveReferenceTest::file_ = nullptr;
bool ComprehensiveReferenceTest::file_loaded_ = false;

// ===========================================================================
// ComprehensiveVJPValidation
//
// The definitive test: for EVERY design variable (k=0..N-1, dim=0,1):
//
// 1. Total VJP check:
//    C++ vjp_total = vjp_bl + vjp_ue + vjp_wake + vjp_wgap
//    should equal lambda_cpp^T @ dR/dx_CS (using C++ lambda, CS dR/dx)
//
// 2. Total sensitivity check:
//    C++ total = partial_cpp - vjp_total_cpp
//    should match: partial_CS - lambda_cpp^T @ dR/dx_CS
//
// 3. Equation-level diagnostic (informational):
//    Split the CS total into BL-row and ue-row contributions to localize bugs.
//    Also split C++ VJP into xi-only, sigma-only, ue, wake, wgap contributions.
// ===========================================================================
TEST_F(ComprehensiveReferenceTest, ComprehensiveVJPValidation) {
    // ---- Load reference data ----
    MatrixXd x_base = readMatrix(getFile(), "comprehensive/foil_x");
    double cl_ref = readScalar<double>(getFile(), "comprehensive/cl_base");
    double cd_ref = readScalar<double>(getFile(), "comprehensive/cd_base");
    MatrixXd dcl_partial_ref = readMatrix(getFile(), "comprehensive/dcl_partial");
    MatrixXd dcd_partial_ref = readMatrix(getFile(), "comprehensive/dcd_partial");

    hsize_t d0, d1, d2;
    auto dR_dx_buf = read3D(getFile(), "comprehensive/dR_dx", d0, d1, d2);
    hsize_t bd0, bd1, bd2;
    auto dR_BL_dx_buf = read3D(getFile(), "comprehensive/dR_BL_dx", bd0, bd1, bd2);
    hsize_t ud0, ud1, ud2;
    auto dR_ue_dx_buf = read3D(getFile(), "comprehensive/dR_ue_dx", ud0, ud1, ud2);

    int N_ref = static_cast<int>(d2);

    // ---- Solve same case in C++ ----
    Cppfoil foil("0012", N_ref - 1);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& state = foil.GetState();
    state.foil.x = x_base;
    state.foil.N = N_ref;
    RecomputeArcLength(state);
    state.oper.viscous = true;
    state.oper.initbl = true;
    state.param.verb = 0;
    SolveViscous(state);
    if (!state.glob.conv) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }
    CalcForce(state);

    int N = state.foil.N;
    int Nsys = state.glob.Nsys;

    double cl_err = std::abs(state.post.cl - cl_ref) / std::max(std::abs(cl_ref), 1e-10);
    double cd_err = std::abs(state.post.cd - cd_ref) / std::max(std::abs(cd_ref), 1e-10);
    EXPECT_LT(cl_err, 1e-3) << "Base cl mismatch";
    EXPECT_LT(cd_err, 1e-2) << "Base cd mismatch";
    std::cout << "Base: cl=" << state.post.cl << " (ref=" << cl_ref
              << ") cd=" << state.post.cd << " (ref=" << cd_ref << ")" << std::endl;

    // ---- Compute C++ adjoint data ----
    CalcForceLin(state);
    CalcCdLin(state);

    VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
    VectorXd dcd_dV = state.post.cd_V;
    for (int i = 0; i < N; ++i)
        dcl_dV(4 * i + 3) = state.post.cl_ue(i);

    auto& lu = *state.glob.lu_factorization;
    VectorXd lambda_cl = lu.transpose().solve(dcl_dV);
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);
    VectorXd lambda_bl_cl = lambda_cl.head(3 * Nsys);
    VectorXd lambda_ue_cl = lambda_cl.tail(Nsys);
    VectorXd lambda_bl_cd = lambda_cd.head(3 * Nsys);
    VectorXd lambda_ue_cd = lambda_cd.tail(Nsys);

    // Compute each VJP chain
    MatrixXd dRue_dx = ComputeDRueDxFoil(state);
    VectorXd vjp_ue_cl = dRue_dx.transpose() * lambda_ue_cl;
    VectorXd vjp_ue_cd = dRue_dx.transpose() * lambda_ue_cd;
    VectorXd vjp_bl_cl = ComputeArcLengthVJP(state, lambda_bl_cl, lambda_ue_cl);
    VectorXd vjp_bl_cd = ComputeArcLengthVJP(state, lambda_bl_cd, lambda_ue_cd);
    VectorXd vjp_wake_cl = ComputeWakeGeomVJP(state, lambda_ue_cl);
    VectorXd vjp_wake_cd = ComputeWakeGeomVJP(state, lambda_ue_cd);
    VectorXd vjp_wgap_cl = ComputeWgapVJP(state, lambda_bl_cl);
    VectorXd vjp_wgap_cd = ComputeWgapVJP(state, lambda_bl_cd);
    VectorXd vjp_xichain_cl = ComputeXiVJP(state, lambda_bl_cl);
    VectorXd vjp_xichain_cd = ComputeXiVJP(state, lambda_bl_cd);

    // Split arc-length VJP into xi-only and sigma-only for diagnostics
    VectorXd vjp_xi_cl = ComputeArcLengthVJP(state, lambda_bl_cl, VectorXd::Zero(Nsys));
    VectorXd vjp_xi_cd = ComputeArcLengthVJP(state, lambda_bl_cd, VectorXd::Zero(Nsys));
    VectorXd vjp_sigma_cl = ComputeArcLengthVJP(state, VectorXd::Zero(3 * Nsys), lambda_ue_cl);
    VectorXd vjp_sigma_cd = ComputeArcLengthVJP(state, VectorXd::Zero(3 * Nsys), lambda_ue_cd);

    // Total C++ VJP
    VectorXd vjp_total_cl = vjp_bl_cl + vjp_ue_cl + vjp_wake_cl + vjp_wgap_cl + vjp_xichain_cl;
    VectorXd vjp_total_cd = vjp_bl_cd + vjp_ue_cd + vjp_wake_cd + vjp_wgap_cd + vjp_xichain_cd;

    // Also get the ViscousSensitivities result
    auto sens = ViscousSensitivities(state);

    // ==================================================================
    // Check 1: Total VJP — definitive chain-rule validation
    //   C++ vjp_total vs lambda_cpp^T @ dR/dx_CS at every node
    // ==================================================================
    std::cout << "\n=== Check 1: Total VJP (lambda^T @ dR/dx) ===" << std::endl;
    int nfail_vjp_cl = 0, nfail_vjp_cd = 0;
    double worst_vjp_cl = 0, worst_vjp_cd = 0;
    int worst_vjp_cl_k = -1, worst_vjp_cd_k = -1;
    int worst_vjp_cl_d = -1, worst_vjp_cd_d = -1;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            int idx = 2 * k + dim;

            // Reference: lambda_cpp^T @ dR/dx_CS
            double vjp_ref_cl = 0, vjp_ref_cd = 0;
            for (int j = 0; j < 4 * Nsys; ++j) {
                double dR = dR_dx_buf[idx3(j, dim, k, 2, N)];
                vjp_ref_cl += lambda_cl(j) * dR;
                vjp_ref_cd += lambda_cd(j) * dR;
            }

            double err_cl = std::abs(vjp_total_cl(idx) - vjp_ref_cl);
            double mag_cl = std::max(std::abs(vjp_ref_cl), 1e-10);
            double rel_cl = err_cl / mag_cl;
            if (rel_cl > TOL_VJP) {
                nfail_vjp_cl++;
                if (rel_cl > worst_vjp_cl) {
                    worst_vjp_cl = rel_cl;
                    worst_vjp_cl_k = k;
                    worst_vjp_cl_d = dim;
                }
            }

            double err_cd = std::abs(vjp_total_cd(idx) - vjp_ref_cd);
            double mag_cd = std::max(std::abs(vjp_ref_cd), 1e-10);
            double rel_cd = err_cd / mag_cd;
            if (rel_cd > TOL_VJP) {
                nfail_vjp_cd++;
                if (rel_cd > worst_vjp_cd) {
                    worst_vjp_cd = rel_cd;
                    worst_vjp_cd_k = k;
                    worst_vjp_cd_d = dim;
                }
            }
        }
    }

    std::cout << "  cl: " << nfail_vjp_cl << "/" << 2 * N
              << " failures, worst=" << worst_vjp_cl
              << " at k=" << worst_vjp_cl_k << " d=" << worst_vjp_cl_d << std::endl;
    std::cout << "  cd: " << nfail_vjp_cd << "/" << 2 * N
              << " failures, worst=" << worst_vjp_cd
              << " at k=" << worst_vjp_cd_k << " d=" << worst_vjp_cd_d << std::endl;

    // Diagnostic: for failing nodes, show equation-level and component breakdown
    if (nfail_vjp_cl > 0 || nfail_vjp_cd > 0) {
        std::cout << "\n  --- Equation-level diagnostic ---" << std::endl;
        std::cout << "  (BL-ref = lambda_bl^T @ dR_BL_CS, ue-ref = lambda_ue^T @ dR_ue_CS)"
                  << std::endl;
        std::cout << "  (xi = xi-chain, sig = sigma_m-chain, ue = AIC-chain, "
                     "wake = wake-geom, wgap = wake-gap)" << std::endl;

        for (int k = 0; k < N; ++k) {
            for (int dim = 0; dim < 2; ++dim) {
                int idx = 2 * k + dim;

                // Compute reference quantities
                double vjp_ref_cl = 0, vjp_ref_cd = 0;
                double vjp_bl_ref_cl = 0, vjp_ue_ref_cl = 0;
                double vjp_bl_ref_cd = 0, vjp_ue_ref_cd = 0;
                for (int j = 0; j < 3 * Nsys; ++j) {
                    double dRbl = dR_BL_dx_buf[idx3(j, dim, k, 2, N)];
                    vjp_bl_ref_cl += lambda_bl_cl(j) * dRbl;
                    vjp_bl_ref_cd += lambda_bl_cd(j) * dRbl;
                }
                for (int j = 0; j < Nsys; ++j) {
                    double dRue = dR_ue_dx_buf[idx3(j, dim, k, 2, N)];
                    vjp_ue_ref_cl += lambda_ue_cl(j) * dRue;
                    vjp_ue_ref_cd += lambda_ue_cd(j) * dRue;
                }
                vjp_ref_cl = vjp_bl_ref_cl + vjp_ue_ref_cl;
                vjp_ref_cd = vjp_bl_ref_cd + vjp_ue_ref_cd;

                double rel_cl = std::abs(vjp_total_cl(idx) - vjp_ref_cl) /
                                std::max(std::abs(vjp_ref_cl), 1e-10);
                double rel_cd = std::abs(vjp_total_cd(idx) - vjp_ref_cd) /
                                std::max(std::abs(vjp_ref_cd), 1e-10);

                if (rel_cl > TOL_VJP || rel_cd > TOL_VJP) {
                    // C++ BL-equation chains: xi + wgap + xichain
                    double cpp_bl_cl = vjp_xi_cl(idx) + vjp_wgap_cl(idx) + vjp_xichain_cl(idx);
                    double cpp_bl_cd = vjp_xi_cd(idx) + vjp_wgap_cd(idx) + vjp_xichain_cd(idx);
                    // C++ ue-equation chains: sigma + ue + wake
                    double cpp_ue_cl = vjp_sigma_cl(idx) + vjp_ue_cl(idx) + vjp_wake_cl(idx);
                    double cpp_ue_cd = vjp_sigma_cd(idx) + vjp_ue_cd(idx) + vjp_wake_cd(idx);

                    double bl_err_cl = std::abs(cpp_bl_cl - vjp_bl_ref_cl) /
                                       std::max(std::abs(vjp_bl_ref_cl), 1e-10);
                    double ue_err_cl = std::abs(cpp_ue_cl - vjp_ue_ref_cl) /
                                       std::max(std::abs(vjp_ue_ref_cl), 1e-10);

                    std::cout << "  k=" << k << " d=" << dim
                              << " total_rel=" << rel_cl
                              << " | BL: cpp=" << cpp_bl_cl
                              << " ref=" << vjp_bl_ref_cl
                              << " rel=" << bl_err_cl
                              << " (xi=" << vjp_xi_cl(idx)
                              << " wgap=" << vjp_wgap_cl(idx)
                              << " xichain=" << vjp_xichain_cl(idx) << ")"
                              << " | ue: cpp=" << cpp_ue_cl
                              << " ref=" << vjp_ue_ref_cl
                              << " rel=" << ue_err_cl
                              << " (sig=" << vjp_sigma_cl(idx)
                              << " aic=" << vjp_ue_cl(idx)
                              << " wake=" << vjp_wake_cl(idx) << ")"
                              << std::endl;
                }
            }
        }
    }

    EXPECT_EQ(nfail_vjp_cl, 0) << "Total VJP cl failures";
    EXPECT_EQ(nfail_vjp_cd, 0) << "Total VJP cd failures";

    // ==================================================================
    // Check 2: Total sensitivity (end-to-end)
    //   C++ dcl/dx vs CS-derived dcl/dx = partial_CS - lambda_cpp^T @ dR/dx_CS
    // ==================================================================
    std::cout << "\n=== Check 2: Total sensitivity ===" << std::endl;
    int nfail_total_cl = 0, nfail_total_cd = 0;
    double worst_total_cl = 0, worst_total_cd = 0;
    int worst_total_cl_k = -1, worst_total_cd_k = -1;

    for (int k = 0; k < N; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            double vjp_ref_cl = 0, vjp_ref_cd = 0;
            for (int j = 0; j < 4 * Nsys; ++j) {
                double dR = dR_dx_buf[idx3(j, dim, k, 2, N)];
                vjp_ref_cl += lambda_cl(j) * dR;
                vjp_ref_cd += lambda_cd(j) * dR;
            }
            // Use CS partials with C++ lambda for a consistent reference
            double total_ref_cl = dcl_partial_ref(dim, k) - vjp_ref_cl;
            double total_ref_cd = dcd_partial_ref(dim, k) - vjp_ref_cd;

            double total_cpp_cl = sens.dcl_dx(dim, k);
            double total_cpp_cd = sens.dcd_dx(dim, k);

            double rel_cl = std::abs(total_cpp_cl - total_ref_cl) /
                            std::max(std::abs(total_ref_cl), 1e-10);
            if (rel_cl > TOL_TOTAL) {
                nfail_total_cl++;
                if (rel_cl > worst_total_cl) {
                    worst_total_cl = rel_cl;
                    worst_total_cl_k = k;
                }
            }

            double rel_cd = std::abs(total_cpp_cd - total_ref_cd) /
                            std::max(std::abs(total_ref_cd), 1e-6);
            if (rel_cd > TOL_TOTAL) {
                nfail_total_cd++;
                if (rel_cd > worst_total_cd) {
                    worst_total_cd = rel_cd;
                    worst_total_cd_k = k;
                }
            }
        }
    }

    std::cout << "  cl: " << nfail_total_cl << "/" << 2 * N
              << " failures, worst=" << worst_total_cl
              << " at k=" << worst_total_cl_k << std::endl;
    std::cout << "  cd: " << nfail_total_cd << "/" << 2 * N
              << " failures, worst=" << worst_total_cd
              << " at k=" << worst_total_cd_k << std::endl;

    if (nfail_total_cl > 0 || nfail_total_cd > 0) {
        std::cout << "  Detailed total failures:" << std::endl;
        for (int k = 0; k < N; ++k) {
            for (int dim = 0; dim < 2; ++dim) {
                int idx = 2 * k + dim;
                double vjp_ref_cl = 0, vjp_ref_cd = 0;
                for (int j = 0; j < 4 * Nsys; ++j) {
                    double dR = dR_dx_buf[idx3(j, dim, k, 2, N)];
                    vjp_ref_cl += lambda_cl(j) * dR;
                    vjp_ref_cd += lambda_cd(j) * dR;
                }
                double total_ref_cl = dcl_partial_ref(dim, k) - vjp_ref_cl;
                double total_ref_cd = dcd_partial_ref(dim, k) - vjp_ref_cd;
                double total_cpp_cl = sens.dcl_dx(dim, k);
                double total_cpp_cd = sens.dcd_dx(dim, k);
                double rel_cl = std::abs(total_cpp_cl - total_ref_cl) /
                                std::max(std::abs(total_ref_cl), 1e-10);
                double rel_cd = std::abs(total_cpp_cd - total_ref_cd) /
                                std::max(std::abs(total_ref_cd), 1e-6);
                if (rel_cl > TOL_TOTAL || rel_cd > TOL_TOTAL) {
                    std::cout << "    k=" << k << " d=" << dim
                              << " cl: cpp=" << total_cpp_cl
                              << " ref=" << total_ref_cl
                              << " rel=" << rel_cl
                              << " | cd: cpp=" << total_cpp_cd
                              << " ref=" << total_ref_cd
                              << " rel=" << rel_cd
                              << std::endl;
                    std::cout << "      components: bl=" << vjp_bl_cl(idx)
                              << " ue=" << vjp_ue_cl(idx)
                              << " wake=" << vjp_wake_cl(idx)
                              << " wgap=" << vjp_wgap_cl(idx)
                              << " partial=" << state.post.dcl_dx_partial(dim, k)
                              << std::endl;
                }
            }
        }
    }

    EXPECT_EQ(nfail_total_cl, 0) << "total dcl failures";
    EXPECT_EQ(nfail_total_cd, 0) << "total dcd failures";

    // ==================================================================
    // Summary
    // ==================================================================
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "  VJP total:    cl=" << nfail_vjp_cl << "/" << 2 * N
              << "  cd=" << nfail_vjp_cd << "/" << 2 * N << std::endl;
    std::cout << "  Sensitivity:  cl=" << nfail_total_cl << "/" << 2 * N
              << "  cd=" << nfail_total_cd << "/" << 2 * N << std::endl;
}
