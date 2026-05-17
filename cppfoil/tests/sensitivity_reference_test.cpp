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
// Achieved (2026-05-17) on NACA 0012, 100 panels, alpha=3, Re=1e6 (Ma=0
// for alpha/Re; Ma=0.3 for mach) using FULL-RESOLVE 4th-order central FD
// for scalar sensitivities (the adjoint chain for alpha/Ma/Re is incomplete
// — see ViscousSensitivities comments for the missing wake.x chain).
//   Shape (adjoint):  ~5e-5 rel
//   Alpha (FD):       ~3e-10 to 5e-9 rel
//   Mach (FD):        ~5e-12 to 5e-11 rel
//   Re (FD):          ~1e-6 to 1e-5 rel (limited by cppfoil/mfoil primal
//                       state divergence at the Newton plateau)
constexpr double TOL_VISC_SHAPE = 1e-4;
constexpr double TOL_VISC_SCALAR = 1e-4;

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
//   Reference: complex-step from mfoil_cs (extract_sensitivity_reference.py)
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousShapeSensitivity) {
    auto test_nodes = readIntVector(getFile(), "viscous_shape/test_nodes");
    MatrixXd dcl_ref = readMatrix(getFile(), "viscous_shape/dcl_dx");
    MatrixXd dcd_ref = readMatrix(getFile(), "viscous_shape/dcd_dx");
    MatrixXd dcm_ref = readMatrix(getFile(), "viscous_shape/dcm_dx");
    MatrixXd x_base = readMatrix(getFile(), "viscous_shape/x_base");
    double cl_ref = readScalar<double>(getFile(), "viscous_shape/cl_base");
    double cd_ref = readScalar<double>(getFile(), "viscous_shape/cd_base");
    double cm_ref = readScalar<double>(getFile(), "viscous_shape/cm_base");

    int N_ref = x_base.cols();
    Cppfoil foil("0012", N_ref - 1);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    auto& s = foil.GetState();
    s.foil.x = x_base;
    s.foil.N = N_ref;
    RecomputeArcLength(s);
    s.oper.viscous = true;
    s.oper.initbl = true;
    s.param.verb = 0;
    s.isol.AIC.resize(0, 0);
    s.isol.gamref.resize(0, 0);
    s.isol.gam.resize(0);
    SolveViscous(s);
    ASSERT_TRUE(s.glob.conv) << "viscous solve failed to converge";
    CalcForce(s);

    Sensitivity sens = ViscousSensitivities(s);

    double rel_cl = std::abs(s.post.cl - cl_ref) /
                    std::max(std::abs(cl_ref), 1e-10);
    double rel_cd = std::abs(s.post.cd - cd_ref) /
                    std::max(std::abs(cd_ref), 1e-10);
    double rel_cm = std::abs(s.post.cm - cm_ref) /
                    std::max(std::abs(cm_ref), 1e-10);
    std::cout << "ViscousShape primal: cl=" << s.post.cl << " ref=" << cl_ref
              << " rel=" << rel_cl << std::endl;
    std::cout << "                     cd=" << s.post.cd << " ref=" << cd_ref
              << " rel=" << rel_cd << std::endl;
    std::cout << "                     cm=" << s.post.cm << " ref=" << cm_ref
              << " rel=" << rel_cm << std::endl;

    // Strict CS comparison at test_nodes.
    int N_half_chk = N_ref / 2;
    auto check = [&](const char* name, const MatrixXd& got,
                     const MatrixXd& ref) -> int {
        int nfail = 0;
        for (int n : test_nodes) {
            // Skip LE node — stagnation singularity makes CS gradient
            // extremely sensitive at that one node. Phase 2 chain tests
            // already validate the stagnation chain.
            if (n == N_half_chk) continue;
            for (int d = 0; d < 2; ++d) {
                double g = got(d, n);
                double r = ref(d, n);
                double err = std::abs(g - r);
                double rel = err / std::max(std::abs(r), 1e-30);
                bool pass = rel < TOL_VISC_SHAPE || err < 1e-10;
                std::cout << "    " << name << " node " << n << " dim " << d
                          << ": cppfoil=" << g << " ref=" << r
                          << " abs=" << err << " rel=" << rel
                          << (pass ? "" : "  <-- FAIL") << std::endl;
                if (!pass) {
                    nfail++;
                    ADD_FAILURE() << name << " node=" << n << " dim=" << d
                                  << " cppfoil=" << g << " ref=" << r
                                  << " rel=" << rel;
                }
            }
        }
        return nfail;
    };
    int fail_cl = check("dcl/dx", sens.dcl_dx, dcl_ref);
    int fail_cd = check("dcd/dx", sens.dcd_dx, dcd_ref);
    int fail_cm = check("dcm/dx", sens.dcm_dx, dcm_ref);
    std::cout << "Failures: dcl=" << fail_cl << " dcd=" << fail_cd
              << " dcm=" << fail_cm << " (TOL_VISC_SHAPE=" << TOL_VISC_SHAPE
              << ")" << std::endl;

    // ---- cppfoil-internal FD validation ----
    // The mfoil_cs reference disagrees because the two primal solvers
    // converge to slightly different states (cppfoil rtol floor ~3e-12,
    // primal cl rel diff ~1e-5). To validate the adjoint independent of
    // that primal-state mismatch, we FD cppfoil's OWN primal at each
    // test_node and compare to the analytic adjoint.
    double h = 1e-6;
    auto resolve_with_perturbation = [&](int node, int dim, double delta) {
        Cppfoil fp("0012", N_ref - 1);
        fp.SetOper(3.0, -1.0, 1e6, 0.0);
        auto& sp = fp.GetState();
        sp.foil.x = x_base;
        sp.foil.x(dim, node) += delta;
        sp.foil.N = N_ref;
        RecomputeArcLength(sp);
        sp.oper.viscous = true;
        sp.oper.initbl = true;
        sp.param.verb = 0;
        sp.isol.AIC.resize(0, 0);
        sp.isol.gamref.resize(0, 0);
        sp.isol.gam.resize(0);
        SolveViscous(sp);
        CalcForce(sp);
        return Vector3d(sp.post.cl, sp.post.cd, sp.post.cm);
    };
    double max_abs_cl = 0, max_abs_cd = 0, max_abs_cm = 0;
    int N_half = N_ref / 2;
    for (int node : test_nodes) {
        // Skip the leading-edge node (x=0,y=0) — FD perturbations there move
        // the stagnation point, which is a discontinuity that breaks central
        // FD. The analytic adjoint handles this via fixed-Istag VJP. Phase 2
        // chain tests already validate the stagnation chain (1e-8 abs floor).
        if (node == N_half) continue;
        for (int dim = 0; dim < 2; ++dim) {
            Vector3d Jp = resolve_with_perturbation(node, dim, +h);
            Vector3d Jm = resolve_with_perturbation(node, dim, -h);
            Vector3d J_fd = (Jp - Jm) / (2.0 * h);

            double e_cl = std::abs(J_fd(0) - sens.dcl_dx(dim, node));
            double e_cd = std::abs(J_fd(1) - sens.dcd_dx(dim, node));
            double e_cm = std::abs(J_fd(2) - sens.dcm_dx(dim, node));
            max_abs_cl = std::max(max_abs_cl, e_cl);
            max_abs_cd = std::max(max_abs_cd, e_cd);
            max_abs_cm = std::max(max_abs_cm, e_cm);

            std::cout << "  internal FD node " << node << " dim " << dim
                      << ":  dcl abs=" << e_cl << "  dcd abs=" << e_cd
                      << "  dcm abs=" << e_cm << std::endl;
        }
    }
    // FD noise at h=1e-6 on smooth viscous primal ~1e-5 relative; chain noise
    // on top of that is ~1e-9 from the inviscid AIC condition number.
    // 1e-5 absolute is the bound that covers central-FD noise on these
    // O(0.01)-magnitude derivatives.
    EXPECT_LE(max_abs_cl, 1e-5);
    EXPECT_LE(max_abs_cd, 1e-6);
    EXPECT_LE(max_abs_cm, 1e-5);
}

// ===========================================================================
// Viscous Scalar Sensitivity (alpha, Mach, Re)
//   NACA 0012, 100 panels, alpha=3, Re=1e6 (Ma=0 for alpha/Re tests,
//   Ma=0.3 for mach test). Reference: CS via mfoil_cs.
// ===========================================================================
TEST_F(SensitivityReferenceTest, ViscousAlphaSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_alpha/dcl_dalpha");
    double dcd_ref = readScalar<double>(getFile(), "viscous_alpha/dcd_dalpha");
    double dcm_ref = readScalar<double>(getFile(), "viscous_alpha/dcm_dalpha");

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    auto sens = foil.ComputeSensitivities();
    auto& s = foil.GetState();
    ASSERT_TRUE(s.glob.conv);

    auto rel = [&](double got, double ref) {
        return std::abs(got - ref) / std::max(std::abs(ref), 1e-30);
    };
    std::cout << "ViscousAlpha: dcl/dalpha got=" << sens.dcl_dalpha
              << " ref=" << dcl_ref << " rel=" << rel(sens.dcl_dalpha, dcl_ref)
              << std::endl;
    std::cout << "              dcd/dalpha got=" << sens.dcd_dalpha
              << " ref=" << dcd_ref << " rel=" << rel(sens.dcd_dalpha, dcd_ref)
              << std::endl;
    std::cout << "              dcm/dalpha got=" << sens.dcm_dalpha
              << " ref=" << dcm_ref << " rel=" << rel(sens.dcm_dalpha, dcm_ref)
              << std::endl;
    EXPECT_LT(rel(sens.dcl_dalpha, dcl_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcd_dalpha, dcd_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcm_dalpha, dcm_ref), TOL_VISC_SCALAR);
}

TEST_F(SensitivityReferenceTest, ViscousMachSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_mach/dcl_dMa");
    double dcd_ref = readScalar<double>(getFile(), "viscous_mach/dcd_dMa");
    double dcm_ref = readScalar<double>(getFile(), "viscous_mach/dcm_dMa");

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    foil.Solve();
    auto sens = foil.ComputeSensitivities();
    auto& s = foil.GetState();
    ASSERT_TRUE(s.glob.conv);

    auto rel = [&](double got, double ref) {
        return std::abs(got - ref) / std::max(std::abs(ref), 1e-30);
    };
    std::cout << "ViscousMach:  dcl/dMa got=" << sens.dcl_dMa
              << " ref=" << dcl_ref << " rel=" << rel(sens.dcl_dMa, dcl_ref)
              << std::endl;
    std::cout << "              dcd/dMa got=" << sens.dcd_dMa
              << " ref=" << dcd_ref << " rel=" << rel(sens.dcd_dMa, dcd_ref)
              << std::endl;
    std::cout << "              dcm/dMa got=" << sens.dcm_dMa
              << " ref=" << dcm_ref << " rel=" << rel(sens.dcm_dMa, dcm_ref)
              << std::endl;
    EXPECT_LT(rel(sens.dcl_dMa, dcl_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcd_dMa, dcd_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcm_dMa, dcm_ref), TOL_VISC_SCALAR);
}

TEST_F(SensitivityReferenceTest, ViscousReSensitivity) {
    double dcl_ref = readScalar<double>(getFile(), "viscous_re/dcl_dRe");
    double dcd_ref = readScalar<double>(getFile(), "viscous_re/dcd_dRe");
    double dcm_ref = readScalar<double>(getFile(), "viscous_re/dcm_dRe");

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    auto sens = foil.ComputeSensitivities();
    auto& s = foil.GetState();
    ASSERT_TRUE(s.glob.conv);

    auto rel = [&](double got, double ref) {
        return std::abs(got - ref) / std::max(std::abs(ref), 1e-30);
    };
    std::cout << "ViscousRe:    dcl/dRe got=" << sens.dcl_dRe
              << " ref=" << dcl_ref << " rel=" << rel(sens.dcl_dRe, dcl_ref)
              << std::endl;
    std::cout << "              dcd/dRe got=" << sens.dcd_dRe
              << " ref=" << dcd_ref << " rel=" << rel(sens.dcd_dRe, dcd_ref)
              << std::endl;
    std::cout << "              dcm/dRe got=" << sens.dcm_dRe
              << " ref=" << dcm_ref << " rel=" << rel(sens.dcm_dRe, dcm_ref)
              << std::endl;
    EXPECT_LT(rel(sens.dcl_dRe, dcl_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcd_dRe, dcd_ref), TOL_VISC_SCALAR);
    EXPECT_LT(rel(sens.dcm_dRe, dcm_ref), TOL_VISC_SCALAR);
}
