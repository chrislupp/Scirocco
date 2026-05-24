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
*/

#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "scirocco/core/inviscid.h"
#include "scirocco/core/viscous_geom_deriv.h"
#include "scirocco/core/wake.h"
#include "scirocco/data/state.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

// Phase 1 tolerance: complex-step gives machine-precision per-operation, but
// the inviscid-solve + wake-build chain accumulates ~N*eps absolute roundoff
// (~1e-13 to 1e-12 for an N=100-panel problem). We accept either an absolute
// or relative threshold, whichever is tighter — i.e.
//   |got - ref|  <=  REL * |ref|  +  ABS
//
// Pure algorithmic bugs would produce errors orders of magnitude larger than
// this floor; numerical-precision noise saturates near these values.
constexpr double TOL_REL = 1e-10;
// Absolute floor reflects the measured noise of the chains:
//   - closed-form chains  : ~5e-12  (cond(AIC) * eps_machine for the panel
//   solve)
//   - iterative wake chain: ~1.5e-11 (above floor grows ~linearly with wake
//                                     station index over Nw=21 steps)
// 2e-11 covers both with a small margin and is still 50x below the 1e-9
// Phase 7 end-to-end target.
constexpr double TOL_ABS = 2e-11;

static bool withinTol(double got, double ref) {
    return std::abs(got - ref) <= TOL_REL * std::abs(ref) + TOL_ABS;
}

// -------------------- HDF5 helpers (same as sensitivity test) ---------------
namespace {
template <typename T>
T readScalar(H5::H5File& file, const std::string& path) {
    H5::DataSet ds = file.openDataSet(path);
    T value;
    ds.read(&value, ds.getDataType());
    return value;
}

VectorXd readVector(H5::H5File& file, const std::string& path) {
    H5::DataSet ds = file.openDataSet(path);
    H5::DataSpace sp = ds.getSpace();
    hsize_t dims[1];
    sp.getSimpleExtentDims(dims);
    std::vector<double> buf(dims[0]);
    ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(buf.data(), dims[0]);
}

MatrixXd readMatrix(H5::H5File& file, const std::string& path) {
    H5::DataSet ds = file.openDataSet(path);
    H5::DataSpace sp = ds.getSpace();
    hsize_t dims[2];
    sp.getSimpleExtentDims(dims);
    std::vector<double> buf(dims[0] * dims[1]);
    ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd m(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j) m(i, j) = buf[i * dims[1] + j];
    return m;
}

std::vector<int> readIntVector(H5::H5File& file, const std::string& path) {
    H5::DataSet ds = file.openDataSet(path);
    H5::DataSpace sp = ds.getSpace();
    hsize_t dims[1];
    sp.getSimpleExtentDims(dims);
    std::vector<int64_t> buf(dims[0]);
    ds.read(buf.data(), H5::PredType::NATIVE_INT64);
    std::vector<int> out(dims[0]);
    for (size_t i = 0; i < dims[0]; ++i) out[i] = static_cast<int>(buf[i]);
    return out;
}
}  // namespace

// -------------------- Fixture: rebuild state matching mfoil_cs --------------
class GeomChainTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool loaded_;

    static void SetUpTestSuite() {
        try {
            file_ = new H5::H5File(GEOM_CHAIN_REFERENCE_FILE, H5F_ACC_RDONLY);
            loaded_ = true;
        } catch (H5::Exception& e) {
            std::cerr << "Failed to open " << GEOM_CHAIN_REFERENCE_FILE << ": "
                      << e.getCDetailMsg() << std::endl;
            loaded_ = false;
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
        if (!loaded_) GTEST_SKIP() << "geom_chain_reference.h5 not loaded";
        loadState();
    }
    H5::H5File& f() { return *file_; }

    // Persistent across all tests in the fixture
    Scirocco foil_{"0012", 100};
    MatrixXd dfoil_s_dx_;
    int N_{0};
    int Nw_{0};

    void loadState() {
        MatrixXd foil_x = readMatrix(f(), "geom_chain/foil_x");
        VectorXd foil_s = readVector(f(), "geom_chain/foil_s");
        double alpha = readScalar<double>(f(), "geom_chain/alpha");
        dfoil_s_dx_ = readMatrix(f(), "geom_chain/dfoil_s_dx");

        N_ = static_cast<int>(foil_x.cols());
        foil_.SetOper(alpha);
        auto& s = foil_.GetState();
        s.foil.x = foil_x;
        s.foil.N = N_;
        s.foil.s =
            foil_s;  // Match mfoil's foil.s exactly (no RecomputeArcLength)
        // Clear cached inviscid quantities so SolveInviscid rebuilds
        s.isol.AIC.resize(0, 0);
        s.isol.gamref.resize(0, 0);
        s.isol.gam.resize(0);
        SolveInviscid(s);  // builds AIC, gamref, gam (does NOT touch foil.s)
        StagPointFind(s);
        BuildWake(s);
        ExtendXiToWake(s);
        IdentifySurfaces(s);
        SetWakeGap(s);
        Nw_ = s.wake.N;
        // Populate glob.U(3, :) with inviscid ue so DSstagDFoilX (which
        // reads ue-based weights from glob.U) can run. At the inviscid
        // state ue = sgnue*gam, which matches the gamma-based stagpoint.
        int Nsys = N_ + Nw_;
        s.glob.U = MatrixXd::Zero(4, Nsys);
        VectorXd ueinv = GetUeInv(s);
        s.glob.U.row(3).head(Nsys) = ueinv.transpose();
    }

    void verifyBase() {
        auto& s = foil_.GetState();
        VectorXd foil_s_ref = readVector(f(), "geom_chain/foil_s");
        double sstag_ref = readScalar<double>(f(), "geom_chain/sstag");
        VectorXd xi_ref = readVector(f(), "geom_chain/xi");
        VectorXd wgap_ref = readVector(f(), "geom_chain/wgap");

        EXPECT_LT((s.foil.s - foil_s_ref).cwiseAbs().maxCoeff(), 1e-13)
            << "foil.s mismatch";
        EXPECT_LT(std::abs(s.isol.sstag - sstag_ref), 1e-13)
            << "sstag mismatch (scirocco=" << s.isol.sstag
            << ", mfoil=" << sstag_ref << ")";
        EXPECT_LT((s.isol.xi - xi_ref).cwiseAbs().maxCoeff(), 1e-12)
            << "xi mismatch";
        EXPECT_LT((s.vsol.wgap - wgap_ref).cwiseAbs().maxCoeff(), 1e-12)
            << "wgap mismatch";
    }

    // Compare element-wise with hybrid abs+rel tolerance, return #failures.
    template <typename A, typename B>
    int compareHybrid(const A& got, const B& ref, const std::string& label,
                      int max_print = 10) {
        int nfail = 0;
        int printed = 0;
        for (Eigen::Index i = 0; i < ref.size(); ++i) {
            double r = ref(i);
            double v = got(i);
            if (!withinTol(v, r)) {
                ++nfail;
                if (printed < max_print) {
                    double err = std::abs(v - r);
                    double rel = err / std::max(std::abs(r), 1e-300);
                    ADD_FAILURE()
                        << label << " idx=" << i << " got=" << v << " ref=" << r
                        << " abs=" << err << " rel=" << rel;
                    ++printed;
                }
            }
        }
        return nfail;
    }

    // 2D version for matrices.
    int compareHybrid2D(const MatrixXd& got, const MatrixXd& ref,
                        const std::string& label, int max_print = 10) {
        int nfail = 0;
        int printed = 0;
        for (Eigen::Index i = 0; i < ref.rows(); ++i) {
            for (Eigen::Index j = 0; j < ref.cols(); ++j) {
                double r = ref(i, j);
                double v = got(i, j);
                if (!withinTol(v, r)) {
                    ++nfail;
                    if (printed < max_print) {
                        double err = std::abs(v - r);
                        double rel = err / std::max(std::abs(r), 1e-300);
                        ADD_FAILURE()
                            << label << "[" << i << "," << j << "] got=" << v
                            << " ref=" << r << " abs=" << err << " rel=" << rel;
                        ++printed;
                    }
                }
            }
        }
        return nfail;
    }
};

H5::H5File* GeomChainTest::file_ = nullptr;
bool GeomChainTest::loaded_ = false;

// ============================================================================
// Sanity: base scirocco state matches the HDF5 reference
// ============================================================================
TEST_F(GeomChainTest, BaseStateMatchesReference) { verifyBase(); }

// ============================================================================
// DFoilSDFoilXCumulative — pure analytical, validate vs CS reference
// ============================================================================
TEST_F(GeomChainTest, DFoilSDFoilXCumulative) {
    auto& s = foil_.GetState();
    MatrixXd a = DFoilSDFoilXCumulative(s.foil);
    MatrixXd r = readMatrix(f(), "geom_chain/dfoil_s_dx");
    int nfail = compareHybrid2D(a, r, "dfoil_s");
    std::cout << "DFoilSDFoilXCumulative: " << nfail << "/" << r.size()
              << " failures" << std::endl;
}

// ============================================================================
// DSstagDFoilX
// ============================================================================
TEST_F(GeomChainTest, DSstagDFoilX) {
    auto& s = foil_.GetState();
    VectorXd a = scirocco::DSstagDFoilX(s, dfoil_s_dx_);
    VectorXd r = readVector(f(), "geom_chain/dsstag_dx");
    int nfail = compareHybrid(a, r, "dsstag_dx");
    std::cout << "DSstagDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}

// ============================================================================
// DWakeSDFoilX
// ============================================================================
TEST_F(GeomChainTest, DWakeSDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DWakeSDFoilX(s, dfoil_s_dx_);
    MatrixXd r = readMatrix(f(), "geom_chain/dwake_s_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    int nfail = compareHybrid2D(a, r, "dwake_s");
    std::cout << "DWakeSDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}

// ============================================================================
// DXiDFoilX (foil + wake combined)
// ============================================================================
TEST_F(GeomChainTest, DXiDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DXiDFoilX(s, dfoil_s_dx_);
    MatrixXd r = readMatrix(f(), "geom_chain/dxi_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    int nfail = compareHybrid2D(a, r, "dxi");
    std::cout << "DXiDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}

// ============================================================================
// DWgapDFoilX
// ============================================================================
TEST_F(GeomChainTest, DWgapDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DWgapDFoilX(s, dfoil_s_dx_);
    MatrixXd r = readMatrix(f(), "geom_chain/dwgap_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    int nfail = compareHybrid2D(a, r, "dwgap");
    std::cout << "DWgapDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}

// ============================================================================
// DWakeXDFoilX — currently STUBBED (returns zero). Re-enable when implemented.
// ============================================================================
TEST_F(GeomChainTest, DWakeXDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DWakeXDFoilX(s);
    MatrixXd r = readMatrix(f(), "geom_chain/dwake_x_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    int nfail = compareHybrid2D(a, r, "dwake_x");
    std::cout << "DWakeXDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}
