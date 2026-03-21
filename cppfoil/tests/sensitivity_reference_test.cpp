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
