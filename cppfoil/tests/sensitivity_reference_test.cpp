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
#include <iostream>
#include <string>
#include <vector>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/data/state.h"

using namespace cppfoil;
using namespace Eigen;

// Tolerance: complex-step reference is machine-precision accurate.
// Any error beyond this indicates a bug in the C++ adjoint.
constexpr double TOL_INVISCID = 1e-10;   // inviscid: machine precision
constexpr double TOL_VISC_SCALAR = 1e-6; // viscous scalar sensitivities
constexpr double TOL_VISC_SHAPE = 1e-6;  // viscous shape sensitivities

// ---------------------------------------------------------------------------
// HDF5 helpers (same pattern as mfoil_reference_test.cpp)
// ---------------------------------------------------------------------------
template <typename T>
T readScalar(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    T value;
    dataset.read(&value, dataset.getDataType());
    return value;
}

VectorXd readVector(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(buffer.data(), dims[0]);
}

MatrixXd readMatrix(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    std::vector<double> buffer(dims[0] * dims[1]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd result(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i) {
        for (hsize_t j = 0; j < dims[1]; ++j) {
            result(i, j) = buffer[i * dims[1] + j];
        }
    }
    return result;
}

std::vector<int> readIntVector(H5::H5File& file, const std::string& path) {
    H5::DataSet dataset = file.openDataSet(path);
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    std::vector<int64_t> buffer(dims[0]);
    dataset.read(buffer.data(), H5::PredType::NATIVE_INT64);
    std::vector<int> result(dims[0]);
    for (size_t i = 0; i < dims[0]; ++i) {
        result[i] = static_cast<int>(buffer[i]);
    }
    return result;
}

// ---------------------------------------------------------------------------
// Test fixture: loads HDF5 once for all tests
// ---------------------------------------------------------------------------
class SensitivityReferenceTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool file_loaded_;

    static void SetUpTestSuite() {
        try {
            file_ = new H5::H5File(SENSITIVITY_REFERENCE_FILE, H5F_ACC_RDONLY);
            file_loaded_ = true;
            std::cout << "Loaded sensitivity reference file: "
                      << SENSITIVITY_REFERENCE_FILE << std::endl;
        } catch (H5::Exception& e) {
            std::cerr << "Failed to open HDF5 file: " << e.getCDetailMsg()
                      << std::endl;
            file_loaded_ = false;
        }
    }

    static void TearDownTestSuite() {
        if (file_) {
            file_->close();
            delete file_;
            file_ = nullptr;
        }
    }

    void SetUp() override {
        if (!file_loaded_) {
            GTEST_SKIP() << "Sensitivity reference file not loaded";
        }
    }

    H5::H5File& getFile() { return *file_; }
};

H5::H5File* SensitivityReferenceTest::file_ = nullptr;
bool SensitivityReferenceTest::file_loaded_ = false;

// ===========================================================================
// Inviscid Shape Sensitivity
//   NACA 2412, 80 panels, alpha=5
//   Replaces: InviscidSensitivity_cl (16 inviscid solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, InviscidShapeSensitivity) {
    // Read reference data
    auto test_nodes = readIntVector(getFile(), "inviscid_shape/test_nodes");
    MatrixXd dcl_dx_ref = readMatrix(getFile(), "inviscid_shape/dcl_dx");
    MatrixXd dcm_dx_ref = readMatrix(getFile(), "inviscid_shape/dcm_dx");
    MatrixXd x_base = readMatrix(getFile(), "inviscid_shape/x_base");
    double cl_ref = readScalar<double>(getFile(), "inviscid_shape/cl_base");

    // Set panel coordinates from reference — these are the exact coordinates
    // the complex-step derivatives were computed on
    int N_ref = x_base.cols();
    Cppfoil foil("2412", N_ref - 1);
    foil.SetOper(5.0);
    auto& state = foil.GetState();
    state.foil.x = x_base;
    state.foil.N = N_ref;
    RecomputeArcLength(state);
    // Clear cached inviscid solution so it resolves on new geometry
    state.isol.AIC.resize(0, 0);
    state.isol.gamref.resize(0, 0);
    state.isol.gam.resize(0);
    SolveInviscid(state);
    CalcForce(state);
    auto sens = InviscidSensitivities(state);

    // Verify base solution consistency
    double cl_err = std::abs(state.post.cl - cl_ref) /
                    std::max(std::abs(cl_ref), 1e-10);
    EXPECT_LT(cl_err, 1e-4) << "Base cl mismatch: cppfoil=" << state.post.cl
                             << " ref=" << cl_ref;

    // Compare adjoint vs complex-step reference at each test node
    int nfail = 0;
    for (int node : test_nodes) {
        if (node >= state.foil.N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            double dcl_adj = sens.dcl_dx(dim, node);
            double dcl_ref_val = dcl_dx_ref(dim, node);

            double err = std::abs(dcl_adj - dcl_ref_val);
            double ref_mag = std::max(std::abs(dcl_ref_val), 1e-10);
            double rel = err / ref_mag;

            if (rel > TOL_INVISCID) {
                nfail++;
                ADD_FAILURE() << "InviscidShape dcl node=" << node << " dim=" << dim
                              << " adj=" << dcl_adj << " ref=" << dcl_ref_val
                              << " rel=" << rel;
            }

            // Also check cm
            double dcm_adj = sens.dcm_dx(dim, node);
            double dcm_ref_val = dcm_dx_ref(dim, node);
            double err_cm = std::abs(dcm_adj - dcm_ref_val);
            double ref_mag_cm = std::max(std::abs(dcm_ref_val), 1e-10);
            double rel_cm = err_cm / ref_mag_cm;

            if (rel_cm > TOL_INVISCID) {
                nfail++;
                ADD_FAILURE() << "InviscidShape dcm node=" << node << " dim=" << dim
                              << " adj=" << dcm_adj << " ref=" << dcm_ref_val
                              << " rel=" << rel_cm;
            }
        }
    }
    std::cout << "InviscidShapeSensitivity: " << nfail << "/" << test_nodes.size() * 4
              << " failures" << std::endl;
}

// ===========================================================================
// Inviscid Alpha Sensitivity
//   NACA 2412, 80 panels, alpha=5
//   Replaces: InviscidAlphaSensitivity (4 inviscid solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, InviscidAlphaSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "inviscid_alpha/dcl_dalpha");
    double dcm_ref = readScalar<double>(getFile(), "inviscid_alpha/dcm_dalpha");

    Cppfoil foil("2412", 80);
    foil.SetOper(5.0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();

    double rel_cl = std::abs(sens.dcl_dalpha - dcl_ref) /
                    std::max(std::abs(dcl_ref), 1e-10);
    double rel_cm = std::abs(sens.dcm_dalpha - dcm_ref) /
                    std::max(std::abs(dcm_ref), 1e-10);

    std::cout << "InviscidAlpha: dcl adj=" << sens.dcl_dalpha << " ref=" << dcl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "InviscidAlpha: dcm adj=" << sens.dcm_dalpha << " ref=" << dcm_ref
              << " rel=" << rel_cm << std::endl;

    EXPECT_LT(rel_cl, TOL_INVISCID) << "dcl/dalpha inviscid";
    EXPECT_LT(rel_cm, TOL_INVISCID) << "dcm/dalpha inviscid";
}

// ===========================================================================
// Inviscid Mach Sensitivity
//   NACA 0012, 100 panels, alpha=3, Ma=0.3
//   Replaces: InviscidMachSensitivity (4 inviscid solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, InviscidMachSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "inviscid_mach/dcl_dMa");
    double dcm_ref = readScalar<double>(getFile(), "inviscid_mach/dcm_dMa");

    double Ma0 = 0.3;
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, Ma0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();

    double rel_cl = std::abs(sens.dcl_dMa - dcl_ref) /
                    std::max(std::abs(dcl_ref), 1e-10);
    double rel_cm = std::abs(sens.dcm_dMa - dcm_ref) /
                    std::max(std::abs(dcm_ref), 1e-10);

    std::cout << "InviscidMach: dcl adj=" << sens.dcl_dMa << " ref=" << dcl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "InviscidMach: dcm adj=" << sens.dcm_dMa << " ref=" << dcm_ref
              << " rel=" << rel_cm << std::endl;

    EXPECT_LT(rel_cl, TOL_INVISCID) << "dcl/dMa inviscid";
    EXPECT_LT(rel_cm, TOL_INVISCID) << "dcm/dMa inviscid";
}

// ===========================================================================
// Viscous Shape Sensitivity
//   NACA 0012, 100 panels, alpha=3, Re=1e6
//   Replaces: ViscousSensitivity (24 viscous solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousShapeSensitivity) {
    auto test_nodes = readIntVector(getFile(), "viscous_shape/test_nodes");
    MatrixXd dcl_dx_ref = readMatrix(getFile(), "viscous_shape/dcl_dx");
    MatrixXd dcd_dx_ref = readMatrix(getFile(), "viscous_shape/dcd_dx");
    MatrixXd x_base = readMatrix(getFile(), "viscous_shape/x_base");
    double cl_ref = readScalar<double>(getFile(), "viscous_shape/cl_base");
    double cd_ref = readScalar<double>(getFile(), "viscous_shape/cd_base");

    // Set panel coordinates from reference — these are the exact coordinates
    // the complex-step derivatives were computed on
    int N_ref = x_base.cols();
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
    auto sens = ViscousSensitivities(state);

    // --- Diagnostic: compare vjp_ue against CS reference ---
    {
        std::string vjp_file =
            std::string(SENSITIVITY_REFERENCE_FILE).substr(
                0, std::string(SENSITIVITY_REFERENCE_FILE).rfind('/')) +
            "/vjp_chain_reference.h5";
        try {
            H5::H5File vf(vjp_file, H5F_ACC_RDONLY);
            auto vjp_nodes = readIntVector(vf, "vjp_chain/test_nodes");

            // Read dRue_dx [Nsys, 2, ntest]
            int Nsys = state.glob.Nsys;
            int ntest = vjp_nodes.size();
            H5::DataSet ds = vf.openDataSet("vjp_chain/dRue_dx_fixed_wake");
            std::vector<double> buf(Nsys * 2 * ntest);
            ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);

            // Compute lambda_ue for cl and cd
            CalcForceLin(state);
            CalcCdLin(state);
            VectorXd dcl_dV = VectorXd::Zero(4 * Nsys);
            VectorXd dcd_dV = state.post.cd_V;
            for (int i = 0; i < state.foil.N; ++i) {
                dcl_dV(4 * i + 3) = state.post.cl_ue(i);
            }
            auto& lu = *state.glob.lu_factorization;
            VectorXd lam_cl = lu.transpose().solve(dcl_dV);
            VectorXd lam_cd = lu.transpose().solve(dcd_dV);
            VectorXd lam_ue_cl = lam_cl.tail(Nsys);
            VectorXd lam_ue_cd = lam_cd.tail(Nsys);

            // Also compute C++ vjp_ue directly
            MatrixXd dRue_dx_cpp = ComputeDRueDxFoil(state);

            std::cout << "VJP_ue comparison (C++ vs CS*lambda):" << std::endl;
            for (int ki = 0; ki < ntest; ++ki) {
                int k = vjp_nodes[ki];
                for (int dim = 0; dim < 2; ++dim) {
                    // CS reference: contract dRue_CS with lambda_ue
                    double vjp_ue_cs_cl = 0, vjp_ue_cs_cd = 0;
                    for (int j = 0; j < Nsys; ++j) {
                        int idx = j * (2 * ntest) + dim * ntest + ki;
                        vjp_ue_cs_cl += lam_ue_cl(j) * buf[idx];
                        vjp_ue_cs_cd += lam_ue_cd(j) * buf[idx];
                    }
                    // C++ vjp_ue
                    double vjp_ue_cpp_cl =
                        dRue_dx_cpp.col(2 * k + dim).dot(lam_ue_cl);
                    double vjp_ue_cpp_cd =
                        dRue_dx_cpp.col(2 * k + dim).dot(lam_ue_cd);

                    std::cout << "  k=" << k << " d=" << dim
                              << " cl: cpp=" << vjp_ue_cpp_cl
                              << " cs=" << vjp_ue_cs_cl
                              << " diff=" << vjp_ue_cpp_cl - vjp_ue_cs_cl
                              << " | cd: cpp=" << vjp_ue_cpp_cd
                              << " cs=" << vjp_ue_cs_cd
                              << " diff=" << vjp_ue_cpp_cd - vjp_ue_cs_cd
                              << std::endl;
                }
            }
            vf.close();
        } catch (...) {
            std::cout << "vjp_chain_reference.h5 not found, skipping VJP_ue check"
                      << std::endl;
        }
    }

    // Verify base solution consistency
    double cl_err = std::abs(state.post.cl - cl_ref) /
                    std::max(std::abs(cl_ref), 1e-10);
    double cd_err = std::abs(state.post.cd - cd_ref) /
                    std::max(std::abs(cd_ref), 1e-10);
    EXPECT_LT(cl_err, 1e-3) << "Base cl mismatch: cppfoil=" << state.post.cl
                             << " ref=" << cl_ref;
    EXPECT_LT(cd_err, 1e-2) << "Base cd mismatch: cppfoil=" << state.post.cd
                             << " ref=" << cd_ref;

    int nfail_cl = 0, nfail_cd = 0;
    for (int node : test_nodes) {
        if (node >= state.foil.N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            // cl shape sensitivity
            double dcl_adj = sens.dcl_dx(dim, node);
            double dcl_ref_val = dcl_dx_ref(dim, node);
            double rel_cl = std::abs(dcl_adj - dcl_ref_val) /
                            std::max(std::abs(dcl_ref_val), 1e-6);

            if (rel_cl > TOL_VISC_SHAPE) {
                nfail_cl++;
                ADD_FAILURE() << "ViscousShape dcl node=" << node << " dim=" << dim
                              << " adj=" << dcl_adj << " ref=" << dcl_ref_val
                              << " rel=" << rel_cl;
            } else {
                std::cout << "  ViscShape dcl OK   node=" << node << " dim=" << dim
                          << " adj=" << dcl_adj << " ref=" << dcl_ref_val
                          << " rel=" << rel_cl << std::endl;
            }

            // cd shape sensitivity
            double dcd_adj = sens.dcd_dx(dim, node);
            double dcd_ref_val = dcd_dx_ref(dim, node);
            double rel_cd = std::abs(dcd_adj - dcd_ref_val) /
                            std::max(std::abs(dcd_ref_val), 1e-3);

            if (rel_cd > TOL_VISC_SHAPE) {
                nfail_cd++;
                ADD_FAILURE() << "ViscousShape dcd node=" << node << " dim=" << dim
                              << " adj=" << dcd_adj << " ref=" << dcd_ref_val
                              << " rel=" << rel_cd;
            } else {
                std::cout << "  ViscShape dcd OK   node=" << node << " dim=" << dim
                          << " adj=" << dcd_adj << " ref=" << dcd_ref_val
                          << " rel=" << rel_cd << std::endl;
            }
        }
    }
    std::cout << "ViscousShapeSensitivity cl: " << nfail_cl << "/" << test_nodes.size() * 2
              << " failures" << std::endl;
    std::cout << "ViscousShapeSensitivity cd: " << nfail_cd << "/" << test_nodes.size() * 2
              << " failures" << std::endl;
}

// ===========================================================================
// Viscous Alpha Sensitivity
//   NACA 0012, 100 panels, alpha=3, Re=1e6
//   Replaces: ViscousAlphaSensitivity (4 viscous solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousAlphaSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_alpha/dcl_dalpha");
    double dcd_ref = readScalar<double>(getFile(), "viscous_alpha/dcd_dalpha");
    double dcm_ref = readScalar<double>(getFile(), "viscous_alpha/dcm_dalpha");

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto sens = foil.ComputeSensitivities();

    double rel_cl = std::abs(sens.dcl_dalpha - dcl_ref) /
                    std::max(std::abs(dcl_ref), 1e-6);
    double rel_cd = std::abs(sens.dcd_dalpha - dcd_ref) /
                    std::max(std::abs(dcd_ref), 1e-6);
    double rel_cm = std::abs(sens.dcm_dalpha - dcm_ref) /
                    std::max(std::abs(dcm_ref), 1e-6);

    std::cout << "ViscousAlpha: dcl adj=" << sens.dcl_dalpha << " ref=" << dcl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "ViscousAlpha: dcd adj=" << sens.dcd_dalpha << " ref=" << dcd_ref
              << " rel=" << rel_cd << std::endl;
    std::cout << "ViscousAlpha: dcm adj=" << sens.dcm_dalpha << " ref=" << dcm_ref
              << " rel=" << rel_cm << std::endl;

    EXPECT_LT(rel_cl, TOL_VISC_SCALAR) << "dcl/dalpha viscous";
    EXPECT_LT(rel_cd, TOL_VISC_SCALAR) << "dcd/dalpha viscous";
    EXPECT_LT(rel_cm, TOL_VISC_SCALAR) << "dcm/dalpha viscous";
}

// ===========================================================================
// Viscous Mach Sensitivity
//   NACA 0012, 100 panels, alpha=3, Re=1e6, Ma=0.3
//   Replaces: ViscousMachSensitivity (4 viscous solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousMachSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_mach/dcl_dMa");
    double dcd_ref = readScalar<double>(getFile(), "viscous_mach/dcd_dMa");
    double dcm_ref = readScalar<double>(getFile(), "viscous_mach/dcm_dMa");

    double Ma0 = 0.3;
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, Ma0);
    foil.Solve();
    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto sens = foil.ComputeSensitivities();

    double rel_cl = std::abs(sens.dcl_dMa - dcl_ref) /
                    std::max(std::abs(dcl_ref), 1e-6);
    double rel_cd = std::abs(sens.dcd_dMa - dcd_ref) /
                    std::max(std::abs(dcd_ref), 1e-6);
    double rel_cm = std::abs(sens.dcm_dMa - dcm_ref) /
                    std::max(std::abs(dcm_ref), 1e-6);

    std::cout << "ViscousMach: dcl adj=" << sens.dcl_dMa << " ref=" << dcl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "ViscousMach: dcd adj=" << sens.dcd_dMa << " ref=" << dcd_ref
              << " rel=" << rel_cd << std::endl;
    std::cout << "ViscousMach: dcm adj=" << sens.dcm_dMa << " ref=" << dcm_ref
              << " rel=" << rel_cm << std::endl;

    EXPECT_LT(rel_cl, TOL_VISC_SCALAR) << "dcl/dMa viscous";
    EXPECT_LT(rel_cd, TOL_VISC_SCALAR) << "dcd/dMa viscous";
    EXPECT_LT(rel_cm, TOL_VISC_SCALAR) << "dcm/dMa viscous";
}

// ===========================================================================
// Viscous Reynolds Sensitivity
//   NACA 0012, 100 panels, alpha=3, Re=1e6
//   Replaces: ViscousReSensitivity (4 viscous solves -> 1 + HDF5 read)
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousReSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_re/dcl_dRe");
    double dcd_ref = readScalar<double>(getFile(), "viscous_re/dcd_dRe");
    double dcm_ref = readScalar<double>(getFile(), "viscous_re/dcm_dRe");

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto sens = foil.ComputeSensitivities();

    double rel_cl = std::abs(sens.dcl_dRe - dcl_ref) /
                    std::max(std::abs(dcl_ref), 1e-6);
    double rel_cd = std::abs(sens.dcd_dRe - dcd_ref) /
                    std::max(std::abs(dcd_ref), 1e-6);
    double rel_cm = std::abs(sens.dcm_dRe - dcm_ref) /
                    std::max(std::abs(dcm_ref), 1e-6);

    std::cout << "ViscousRe: dcl adj=" << sens.dcl_dRe << " ref=" << dcl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "ViscousRe: dcd adj=" << sens.dcd_dRe << " ref=" << dcd_ref
              << " rel=" << rel_cd << std::endl;
    std::cout << "ViscousRe: dcm adj=" << sens.dcm_dRe << " ref=" << dcm_ref
              << " rel=" << rel_cm << std::endl;

    EXPECT_LT(rel_cl, TOL_VISC_SCALAR) << "dcl/dRe viscous";
    EXPECT_LT(rel_cd, TOL_VISC_SCALAR) << "dcd/dRe viscous";
    EXPECT_LT(rel_cm, TOL_VISC_SCALAR) << "dcm/dRe viscous";
}

// ===========================================================================
// Wake Geometry Intermediate Validation
//   Validates mu_x/mu_t (Steps 1/1b) against CS reference dRue/dwake data.
//   Same case: NACA 0012, 100 panels, alpha=3, Re=1e6
// ===========================================================================
TEST_F(SensitivityReferenceTest, WakeGeomIntermediates) {
    // Load CS wake geometry reference
    H5::H5File wg_file(
        std::string(SENSITIVITY_REFERENCE_FILE).substr(
            0, std::string(SENSITIVITY_REFERENCE_FILE).rfind('/')) +
            "/wake_geom_reference.h5",
        H5F_ACC_RDONLY);

    // Read dRue/dwake_x [Nsys, 2, Nw] and dRue/dwake_t [Nsys, 2, Nw]
    // NumPy stores row-major: index [j, d, i] at offset j*(2*Nw) + d*Nw + i
    auto read3D = [&](const std::string& path, int d0, int d1, int d2)
        -> std::vector<double> {
        H5::DataSet ds = wg_file.openDataSet(path);
        std::vector<double> buf(d0 * d1 * d2);
        ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
        return buf;
    };

    // Solve same case
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto& state = foil.GetState();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    std::cout << "WakeGeomIntermediates: N=" << N << " Nw=" << Nw
              << " Nsys=" << Nsys << std::endl;

    // Compare full state against mfoil reference
    {
        MatrixXd U_ref = readMatrix(wg_file, "wake_geom/U_full");
        VectorXd gam_ref = readVector(wg_file, "wake_geom/gam");
        VectorXd s_ref = readVector(wg_file, "wake_geom/foil_s");

        // U state: [4 x Nsys]
        double U_err = 0.0;
        int worst_r = -1, worst_i = -1;
        for (int i = 0; i < Nsys; ++i) {
            for (int r = 0; r < 4; ++r) {
                double ref = U_ref(r, i);
                double cpp = state.glob.U(r, i);
                double rel = std::abs(ref - cpp) / std::max(std::abs(ref), 1e-15);
                if (rel > U_err) { U_err = rel; worst_r = r; worst_i = i; }
            }
        }
        std::cout << "  max |U| rel_err = " << U_err
                  << " at row=" << worst_r << " col=" << worst_i
                  << " ref=" << U_ref(worst_r, worst_i)
                  << " cpp=" << state.glob.U(worst_r, worst_i) << std::endl;
        // Per-row max relative error (skip row 2 near-zero values)
        const char* row_names[] = {"th", "ds", "sa", "ue"};
        for (int r = 0; r < 4; ++r) {
            double max_rel = 0.0;
            int max_i = -1;
            for (int i = 0; i < Nsys; ++i) {
                double ref = U_ref(r, i);
                double cpp = state.glob.U(r, i);
                double denom = std::max(std::abs(ref), (r == 2) ? 1e-6 : 1e-15);
                double rel = std::abs(ref - cpp) / denom;
                if (rel > max_rel) { max_rel = rel; max_i = i; }
            }
            std::cout << "  U row " << row_names[r]
                      << ": max_rel=" << max_rel << " at i=" << max_i
                      << " ref=" << U_ref(r, max_i)
                      << " cpp=" << state.glob.U(r, max_i) << std::endl;
        }

        // gamma
        double gam_err = (state.isol.gam - gam_ref).norm() / gam_ref.norm();
        std::cout << "  gam rel_err = " << gam_err << std::endl;

        // foil.s
        double s_err = (state.foil.s - s_ref).norm() / s_ref.norm();
        std::cout << "  foil.s rel_err = " << s_err << std::endl;

        // cl, cd
        std::cout << "  cl=" << state.post.cl
                  << " ref=3.162623310113019e-01"
                  << " err=" << std::abs(state.post.cl - 3.162623310113019e-01) << std::endl;
        std::cout << "  cd=" << state.post.cd
                  << " ref=6.329823768574360e-03"
                  << " err=" << std::abs(state.post.cd - 6.329823768574360e-03) << std::endl;
    }

    // Compute cd adjoint → lambda_ue
    CalcForceLin(state);
    CalcCdLin(state);
    VectorXd dcd_dV = state.post.cd_V;
    auto& lu = *state.glob.lu_factorization;
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);
    VectorXd lambda_ue = lambda_cd.tail(Nsys);

    // Compute C++ mu_x/mu_t
    MatrixXd mu_x_cpp, mu_t_cpp;
    VectorXd vjp_wake = ComputeWakeGeomVJP(state, lambda_ue,
                                            mu_x_cpp, mu_t_cpp);

    // Read CS dRue/dwake_x and dRue/dwake_t
    auto dwx_buf = read3D("wake_geom/dRue_dwake_x", Nsys, 2, Nw);
    auto dwt_buf = read3D("wake_geom/dRue_dwake_t", Nsys, 2, Nw);

    // Compute reference mu_x and mu_t by contracting with lambda_ue:
    //   mu_x_ref(dim, i) = sum_j lambda_ue(j) * dRue_dwx(j, dim, i)
    // Note: Steps 1/1b compute the NEGATIVE of this because R_ue terms are
    // subtracted: R_ue = ue - ueinv - sgnue*ue_sigma*phi
    // and the mu_x/mu_t capture ∂R_ue/∂wake directly (sign included).
    MatrixXd mu_x_ref = MatrixXd::Zero(2, Nw);
    MatrixXd mu_t_ref = MatrixXd::Zero(2, Nw);

    for (int i = 0; i < Nw; ++i) {
        for (int d = 0; d < 2; ++d) {
            double sum_x = 0.0, sum_t = 0.0;
            for (int j = 0; j < Nsys; ++j) {
                int idx = j * (2 * Nw) + d * Nw + i;
                sum_x += lambda_ue(j) * dwx_buf[idx];
                sum_t += lambda_ue(j) * dwt_buf[idx];
            }
            mu_x_ref(d, i) = sum_x;
            mu_t_ref(d, i) = sum_t;
        }
    }

    // Compare mu_x, split into foil-row and wake-row contributions
    // mu_x_foil(d,i) = sum_{j<N} lambda_ue(j) * dRue(j, d, i)
    // mu_x_wake(d,i) = sum_{j>=N} lambda_ue(j) * dRue(j, d, i)
    MatrixXd mu_x_ref_foil = MatrixXd::Zero(2, Nw);
    MatrixXd mu_x_ref_wake = MatrixXd::Zero(2, Nw);
    for (int i = 0; i < Nw; ++i) {
        for (int d = 0; d < 2; ++d) {
            double sf = 0.0, sw = 0.0;
            for (int j = 0; j < N; ++j) {
                sf += lambda_ue(j) * dwx_buf[j * (2 * Nw) + d * Nw + i];
            }
            for (int j = N; j < Nsys; ++j) {
                sw += lambda_ue(j) * dwx_buf[j * (2 * Nw) + d * Nw + i];
            }
            mu_x_ref_foil(d, i) = sf;
            mu_x_ref_wake(d, i) = sw;
        }
    }

    std::cout << "  mu_x comparison (foil-row | wake-row | total):" << std::endl;
    int n_fail_x = 0;
    for (int i = 0; i < Nw; ++i) {
        double ref_norm = mu_x_ref.col(i).norm();
        double cpp_norm = mu_x_cpp.col(i).norm();
        double err_norm = (mu_x_ref.col(i) - mu_x_cpp.col(i)).norm();
        double rel = err_norm / std::max(ref_norm, 1e-15);
        double foil_ref = mu_x_ref_foil.col(i).norm();
        double wake_ref = mu_x_ref_wake.col(i).norm();
        if (rel > 1e-4 || i < 5 || i == Nw - 1) {
            std::cout << "    i=" << i
                      << " foil_ref=" << foil_ref
                      << " wake_ref=" << wake_ref
                      << " |ref|=" << ref_norm
                      << " |cpp|=" << cpp_norm
                      << " rel=" << rel
                      << (rel > 1e-4 ? " FAIL" : " ok") << std::endl;
        }
        if (rel > 1e-4 && ref_norm > 1e-12) n_fail_x++;
    }

    // Compare mu_t
    std::cout << "  mu_t comparison (C++ vs CS reference):" << std::endl;
    int n_fail_t = 0;
    for (int i = 0; i < Nw; ++i) {
        double ref_norm = mu_t_ref.col(i).norm();
        double cpp_norm = mu_t_cpp.col(i).norm();
        double err_norm = (mu_t_ref.col(i) - mu_t_cpp.col(i)).norm();
        double rel = err_norm / std::max(ref_norm, 1e-15);
        if (rel > 1e-4 || i < 3 || i == Nw - 1) {
            std::cout << "    i=" << i << " |ref|=" << ref_norm
                      << " |cpp|=" << cpp_norm << " rel=" << rel
                      << (rel > 1e-4 ? " FAIL" : " ok") << std::endl;
        }
        if (rel > 1e-4 && ref_norm > 1e-12) n_fail_t++;
    }

    EXPECT_EQ(n_fail_x, 0) << "mu_x mismatch at " << n_fail_x << " wake nodes";
    EXPECT_EQ(n_fail_t, 0) << "mu_t mismatch at " << n_fail_t << " wake nodes";
}
