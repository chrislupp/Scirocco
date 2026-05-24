/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/

#include "scirocco/core/viscous_residual_deriv.h"

#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/viscous_geom_deriv.h"
#include "scirocco/core/wake.h"
#include "scirocco/data/state.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

// Same hybrid tolerance reasoning as the Phase 1 chain tests:
// Tolerance reasoning (measured 2026-05-17):
//   DRblDFoilX: max abs 3.13e-8 at the stagnation rows, but those entries
//     have |ref| ~ 8e3 so rel ~ 4e-12 — TOL_REL * |ref| absorbs them. Worst
//     rel anywhere with |ref|>1e-6 is 3.08e-8. Non-stag rows max abs 4e-10.
//   DRueDFoilX: max abs 1.22e-9 at foil rows (inviscid AIC noise floor),
//     9.65e-10 at wake rows, 1.4e-11 at last wake row.
// Hybrid bound |err| <= TOL_REL*|ref| + TOL_ABS with these values gives
// 0 failures on both tests. Cannot tighten TOL_REL below ~5e-8 without
// DRblDFoilX failures at the worst-rel entry.
constexpr double TOL_REL = 1e-7;
constexpr double TOL_ABS = 5e-9;

static bool withinTol(double got, double ref) {
    return std::abs(got - ref) <= TOL_REL * std::abs(ref) + TOL_ABS;
}

// -------------------- HDF5 helpers ------------------------------------------
namespace {
template <typename T>
T readScalar(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    T v;
    d.read(&v, d.getDataType());
    return v;
}
VectorXd readVector(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[1];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    return Map<VectorXd>(b.data(), dims[0]);
}
MatrixXd readMatrix(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[2];
    d.getSpace().getSimpleExtentDims(dims);
    std::vector<double> b(dims[0] * dims[1]);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd m(dims[0], dims[1]);
    for (hsize_t i = 0; i < dims[0]; ++i)
        for (hsize_t j = 0; j < dims[1]; ++j) m(i, j) = b[i * dims[1] + j];
    return m;
}

// Read a [D0, D1, D2] tensor stored as flat 1D, return [D0, D1*D2] matrix
// flattened in C order so column index = j*D2 + k. Useful when comparing the
// reference's [3*Nsys, 2, N] tensors to our analytical [3*Nsys, 2N] matrices
// with column convention 2*k + dim.
MatrixXd read3DAsMatrix(H5::H5File& f, const std::string& p) {
    H5::DataSet d = f.openDataSet(p);
    hsize_t dims[3];
    d.getSpace().getSimpleExtentDims(dims);
    size_t total = dims[0] * dims[1] * dims[2];
    std::vector<double> b(total);
    d.read(b.data(), H5::PredType::NATIVE_DOUBLE);
    // Source layout (C order): index (i, j, k) → i*dims[1]*dims[2] + j*dims[2]
    // + k Target layout: row = i, col = 2*k + j (i.e. swap j and k in the
    // column).
    MatrixXd m(dims[0], dims[1] * dims[2]);
    for (hsize_t i = 0; i < dims[0]; ++i) {
        for (hsize_t j = 0; j < dims[1]; ++j) {
            for (hsize_t k = 0; k < dims[2]; ++k) {
                m(i, 2 * k + j) = b[i * dims[1] * dims[2] + j * dims[2] + k];
            }
        }
    }
    return m;
}
}  // namespace

// ---------------------------------------------------------------------------
// Fixture: load comprehensive_reference.h5 and rebuild scirocco viscous state
// to the same point mfoil was at when the CS reference was extracted.
// ---------------------------------------------------------------------------
class ResidualDerivTest : public ::testing::Test {
   protected:
    static H5::H5File* file_;
    static bool loaded_;

    static void SetUpTestSuite() {
        try {
            file_ =
                new H5::H5File(COMPREHENSIVE_REFERENCE_FILE, H5F_ACC_RDONLY);
            loaded_ = true;
        } catch (H5::Exception& e) {
            std::cerr << "Failed to open " << COMPREHENSIVE_REFERENCE_FILE
                      << ": " << e.getCDetailMsg() << std::endl;
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
        if (!loaded_) GTEST_SKIP() << "comprehensive_reference.h5 not loaded";
        loadState();
    }
    H5::H5File& f() { return *file_; }

    Scirocco foil_{"0012", 100};
    MatrixXd dfoil_s_dx_;
    int N_{0};
    int Nw_{0};
    int Nsys_{0};

    void loadState() {
        MatrixXd foil_x = readMatrix(f(), "comprehensive/foil_x");
        MatrixXd U_ref = readMatrix(f(), "comprehensive/U");
        VectorXd wgap_ref = readVector(f(), "comprehensive/wgap");
        N_ = static_cast<int>(foil_x.cols());

        // IMPORTANT: extract_comprehensive_reference uses CUMULATIVE-panel-norm
        // foil.s during its CS perturbations (via recompute_geom), so the
        // captured d(R_BL)/d(foil.x) is along that path. The HDF5 base foil.s
        // is mfoil's SPLINE value, which is ~4e-5 different and leads to
        // 1e-8 chain mismatches at the TE if used in scirocco. Recompute foil.s
        // here as cumulative panel norms to match the CS path exactly.
        foil_.SetOper(3.0, -1.0, 1e6, 0.0);
        auto& s = foil_.GetState();
        s.foil.x = foil_x;
        s.foil.N = N_;
        RecomputeArcLength(s);  // cumulative panel norms (matches CS chain)
        s.oper.viscous = true;
        s.oper.initbl = true;
        s.param.verb = 0;
        s.isol.AIC.resize(0, 0);
        s.isol.gamref.resize(0, 0);
        s.isol.gam.resize(0);
        SolveViscous(s);  // converges to its own state; we'll override below
        ASSERT_TRUE(s.glob.conv) << "viscous solve did not converge";

        // Override U, wgap to match the HDF5 reference. UE-BASED xi:
        // after overriding U=U_ref, also run StagPointMove to refresh
        // sstag/sstag_ue/xi from U_ref's ue values at the current cumulative
        // foil.s — this mirrors what extract_comprehensive_reference.py does
        // in solve_base() (mfoil.stagpoint_move after recompute_geom).
        ASSERT_EQ(s.glob.U.cols(), U_ref.cols());
        s.glob.U = U_ref;
        ASSERT_EQ(static_cast<int>(wgap_ref.size()),
                  static_cast<int>(s.vsol.wgap.size()));
        s.vsol.wgap = wgap_ref;
        StagPointMove(s);
        s.glob.realloc = true;
        BuildGlobSys(s);

        Nw_ = s.wake.N;
        Nsys_ = N_ + Nw_;

        // d(foil.s)/d(foil.x) — closed-form cumulative version matches the
        // CS perturbation chain in extract_comprehensive_reference.py.
        dfoil_s_dx_ = DFoilSDFoilXCumulative(s.foil);
    }

    int compareHybrid2D(const MatrixXd& got, const MatrixXd& ref,
                        const std::string& label, int max_print = 5000) {
        int nfail = 0, printed = 0;
        for (Eigen::Index i = 0; i < ref.rows(); ++i) {
            for (Eigen::Index j = 0; j < ref.cols(); ++j) {
                if (!withinTol(got(i, j), ref(i, j))) {
                    ++nfail;
                    if (printed < max_print) {
                        double err = std::abs(got(i, j) - ref(i, j));
                        double rel =
                            err / std::max(std::abs(ref(i, j)), 1e-300);
                        ADD_FAILURE()
                            << label << "[" << i << "," << j
                            << "] got=" << got(i, j) << " ref=" << ref(i, j)
                            << " abs=" << err << " rel=" << rel;
                        ++printed;
                    }
                }
            }
        }
        return nfail;
    }
};

H5::H5File* ResidualDerivTest::file_ = nullptr;
bool ResidualDerivTest::loaded_ = false;

// ============================================================================
// DRblDFoilX — d(R_BL)/d(foil.x) at fixed U
// ============================================================================
TEST_F(ResidualDerivTest, DRblDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DRblDFoilX(s, dfoil_s_dx_);
    MatrixXd r = read3DAsMatrix(f(), "comprehensive/dR_BL_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    // Diagnostics: surface the actual max abs at the stagnation row vs.
    // elsewhere so we can size tolerances against real chain noise.
    double overall_max = 0, stag_row_max = 0, other_max = 0;
    int Istag = s.isol.Istag[0];
    int Ig_stag = 3 * Istag;
    Eigen::Index best_i = 0, best_j = 0;
    for (Eigen::Index i = 0; i < r.rows(); ++i) {
        for (Eigen::Index j = 0; j < r.cols(); ++j) {
            double e = std::abs(a(i, j) - r(i, j));
            if (e > overall_max) {
                overall_max = e;
                best_i = i;
                best_j = j;
            }
            if (i >= Ig_stag && i < Ig_stag + 6) {
                stag_row_max = std::max(stag_row_max, e);
            } else {
                other_max = std::max(other_max, e);
            }
        }
    }
    std::cout << "DRblDFoilX: max abs = " << overall_max << " (stag rows "
              << Ig_stag << ".." << Ig_stag + 5 << " max=" << stag_row_max
              << ", other rows max=" << other_max << ")" << std::endl;
    double worst_ref = r(best_i, best_j);
    std::cout << "  Worst entry: [" << best_i << "," << best_j
              << "] ref=" << worst_ref
              << " rel=" << overall_max / std::max(std::abs(worst_ref), 1e-300)
              << std::endl;
    // Worst rel over all entries with |ref| > 1e-6
    double worst_rel = 0;
    Eigen::Index wr_i = 0, wr_j = 0;
    for (Eigen::Index i = 0; i < r.rows(); ++i) {
        for (Eigen::Index j = 0; j < r.cols(); ++j) {
            double ref = std::abs(r(i, j));
            if (ref < 1e-6) continue;
            double rel = std::abs(a(i, j) - r(i, j)) / ref;
            if (rel > worst_rel) {
                worst_rel = rel;
                wr_i = i;
                wr_j = j;
            }
        }
    }
    std::cout << "  Worst rel (|ref|>1e-6): " << worst_rel << " at [" << wr_i
              << "," << wr_j << "] ref=" << r(wr_i, wr_j)
              << " abs=" << std::abs(a(wr_i, wr_j) - r(wr_i, wr_j))
              << std::endl;
    int nfail = compareHybrid2D(a, r, "dR_BL");
    std::cout << "DRblDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}

// ============================================================================
// DRueDFoilX — full analytical (A1, A2, B1, B2a-c, C). Reports per-region
// max abs error and validates against CS reference using the same hybrid
// tolerance as DRblDFoilX.
// ============================================================================
TEST_F(ResidualDerivTest, DRueDFoilX) {
    auto& s = foil_.GetState();
    MatrixXd a = scirocco::DRueDFoilX(s, dfoil_s_dx_);
    MatrixXd r = read3DAsMatrix(f(), "comprehensive/dR_ue_dx");
    ASSERT_EQ(a.rows(), r.rows());
    ASSERT_EQ(a.cols(), r.cols());
    // Diagnostic: report per-row max abs error so we see what's missing.
    double overall_max = 0, foil_max = 0, wake_max = 0, wake_inner_max = 0;
    int n_foil_above = 0, n_wake_above = 0;
    int last_row = r.rows() - 1;
    double last_row_max =
        (a.row(last_row) - r.row(last_row)).cwiseAbs().maxCoeff();
    for (int i = 0; i < r.rows(); ++i) {
        double row_max = (a.row(i) - r.row(i)).cwiseAbs().maxCoeff();
        overall_max = std::max(overall_max, row_max);
        if (i < N_) {
            foil_max = std::max(foil_max, row_max);
            if (row_max > 1e-6) n_foil_above++;
        } else {
            wake_max = std::max(wake_max, row_max);
            if (i < last_row)
                wake_inner_max = std::max(wake_inner_max, row_max);
            if (row_max > 1e-6) n_wake_above++;
        }
    }
    std::cout << "  Last wake row [" << last_row
              << "] (predictor chain): " << last_row_max << std::endl;
    std::cout << "  Wake rows excluding last: " << wake_inner_max << std::endl;
    std::cout << "DRueDFoilX: overall max abs error = " << overall_max
              << std::endl;
    std::cout << "  Foil rows: max=" << foil_max
              << ", fails(>1e-6)=" << n_foil_above << "/" << N_ << std::endl;
    std::cout << "  Wake rows: max=" << wake_max
              << ", fails(>1e-6)=" << n_wake_above << "/" << Nw_ << std::endl;
    // Locate worst entry
    Eigen::Index best_i = 0, best_j = 0;
    double maxabs = 0;
    for (int i = N_; i < r.rows(); ++i) {
        for (int j = 0; j < r.cols(); ++j) {
            double e = std::abs(a(i, j) - r(i, j));
            if (e > maxabs) {
                maxabs = e;
                best_i = i;
                best_j = j;
            }
        }
    }
    std::cout << "  Worst wake entry: [" << best_i << "," << best_j
              << "] got=" << a(best_i, best_j) << " ref=" << r(best_i, best_j)
              << " abs=" << maxabs << std::endl;
    int nfail = compareHybrid2D(a, r, "dR_ue");
    std::cout << "DRueDFoilX: " << nfail << "/" << r.size() << " failures"
              << std::endl;
}
