/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").

    Replaced original FD convergence-rate ping tests with mfoil_cs
    complex-step references (extract_sensitivity_ping_cs_reference.py).
*/
#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "cppfoil/config.h"
#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/panel_influence_deriv.h"
#include "cppfoil/cppfoil.h"

using namespace cppfoil;
using namespace Eigen;

static constexpr double TOL_REL = 1e-9;
static constexpr double TOL_ABS = 1e-11;

static bool withinTol(double got, double ref) {
    return std::abs(got - ref) <= TOL_REL * std::abs(ref) + TOL_ABS;
}

namespace {
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
}  // namespace

class SensitivityPingTest : public ::testing::Test {
   protected:
    static std::unique_ptr<H5::H5File> file_;
    static void SetUpTestSuite() {
        file_ = std::make_unique<H5::H5File>(SENSITIVITY_PING_REFERENCE_FILE,
                                             H5F_ACC_RDONLY);
    }
    static void TearDownTestSuite() {
        if (file_) file_->close();
    }
};
std::unique_ptr<H5::H5File> SensitivityPingTest::file_;

static const Matrix2d kXj = (Matrix2d() << 0.3, 0.6, 0.05, -0.02).finished();
static const Vector2d kXi(0.45, 0.15);

TEST_F(SensitivityPingTest, PanelLinVortexStreamDeriv_Xj) {
    auto deriv = PanelLinVortexStreamDeriv(kXj, kXi);
    VectorXd cs_a = readVector(*file_, "panels/LinVortexStream/a_Xj");
    VectorXd cs_b = readVector(*file_, "panels/LinVortexStream/b_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        if (!withinTol(deriv.a_Xj(k), cs_a(k))) {
            ++nfail;
            ADD_FAILURE() << "a_Xj[" << k << "] analytic=" << deriv.a_Xj(k)
                          << " cs=" << cs_a(k);
        }
        if (!withinTol(deriv.b_Xj(k), cs_b(k))) {
            ++nfail;
            ADD_FAILURE() << "b_Xj[" << k << "] analytic=" << deriv.b_Xj(k)
                          << " cs=" << cs_b(k);
        }
    }
    std::cout << "PanelLinVortexStream Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelLinVortexStreamDeriv_xi) {
    auto deriv = PanelLinVortexStreamDeriv(kXj, kXi);
    VectorXd cs_a = readVector(*file_, "panels/LinVortexStream/a_xi");
    VectorXd cs_b = readVector(*file_, "panels/LinVortexStream/b_xi");
    int nfail = 0;
    for (int k = 0; k < 2; ++k) {
        if (!withinTol(deriv.a_xi(k), cs_a(k))) {
            ++nfail;
            ADD_FAILURE() << "a_xi[" << k << "] analytic=" << deriv.a_xi(k)
                          << " cs=" << cs_a(k);
        }
        if (!withinTol(deriv.b_xi(k), cs_b(k))) {
            ++nfail;
            ADD_FAILURE() << "b_xi[" << k << "] analytic=" << deriv.b_xi(k)
                          << " cs=" << cs_b(k);
        }
    }
    std::cout << "PanelLinVortexStream xi: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelConstSourceStreamDeriv_Xj) {
    auto deriv = PanelConstSourceStreamDeriv(kXj, kXi);
    VectorXd cs = readVector(*file_, "panels/ConstSourceStream/val_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        if (!withinTol(deriv.val_Xj(k), cs(k))) {
            ++nfail;
            ADD_FAILURE() << "val_Xj[" << k << "] analytic=" << deriv.val_Xj(k)
                          << " cs=" << cs(k);
        }
    }
    std::cout << "PanelConstSourceStream Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelLinVortexVelocityDeriv_Xj) {
    auto deriv = PanelLinVortexVelocityDeriv(kXj, kXi);
    MatrixXd cs_a = readMatrix(*file_, "panels/LinVortexVelocity/a_Xj");
    MatrixXd cs_b = readMatrix(*file_, "panels/LinVortexVelocity/b_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        for (int d = 0; d < 2; ++d) {
            if (!withinTol(deriv.a_Xj(d, k), cs_a(d, k))) {
                ++nfail;
                ADD_FAILURE() << "a_Xj[" << d << "," << k
                              << "] analytic=" << deriv.a_Xj(d, k)
                              << " cs=" << cs_a(d, k);
            }
            if (!withinTol(deriv.b_Xj(d, k), cs_b(d, k))) {
                ++nfail;
                ADD_FAILURE() << "b_Xj[" << d << "," << k
                              << "] analytic=" << deriv.b_Xj(d, k)
                              << " cs=" << cs_b(d, k);
            }
        }
    }
    std::cout << "PanelLinVortexVelocity Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelConstSourceVelocityDeriv_Xj) {
    auto deriv = PanelConstSourceVelocityDeriv(kXj, kXi);
    MatrixXd cs = readMatrix(*file_, "panels/ConstSourceVelocity/val_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        for (int d = 0; d < 2; ++d) {
            if (!withinTol(deriv.val_Xj(d, k), cs(d, k))) {
                ++nfail;
                ADD_FAILURE() << "val_Xj[" << d << "," << k
                              << "] analytic=" << deriv.val_Xj(d, k)
                              << " cs=" << cs(d, k);
            }
        }
    }
    std::cout << "PanelConstSourceVelocity Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelLinSourceStreamDeriv_Xj) {
    auto deriv = PanelLinSourceStreamDeriv(kXj, kXi);
    VectorXd cs_a = readVector(*file_, "panels/LinSourceStream/a_Xj");
    VectorXd cs_b = readVector(*file_, "panels/LinSourceStream/b_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        if (!withinTol(deriv.a_Xj(k), cs_a(k))) {
            ++nfail;
            ADD_FAILURE() << "a_Xj[" << k << "] analytic=" << deriv.a_Xj(k)
                          << " cs=" << cs_a(k);
        }
        if (!withinTol(deriv.b_Xj(k), cs_b(k))) {
            ++nfail;
            ADD_FAILURE() << "b_Xj[" << k << "] analytic=" << deriv.b_Xj(k)
                          << " cs=" << cs_b(k);
        }
    }
    std::cout << "PanelLinSourceStream Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, PanelLinSourceVelocityDeriv_Xj) {
    auto deriv = PanelLinSourceVelocityDeriv(kXj, kXi);
    MatrixXd cs_a = readMatrix(*file_, "panels/LinSourceVelocity/a_Xj");
    MatrixXd cs_b = readMatrix(*file_, "panels/LinSourceVelocity/b_Xj");
    int nfail = 0;
    for (int k = 0; k < 4; ++k) {
        for (int d = 0; d < 2; ++d) {
            if (!withinTol(deriv.a_Xj(d, k), cs_a(d, k))) {
                ++nfail;
                ADD_FAILURE() << "a_Xj[" << d << "," << k
                              << "] analytic=" << deriv.a_Xj(d, k)
                              << " cs=" << cs_a(d, k);
            }
            if (!withinTol(deriv.b_Xj(d, k), cs_b(d, k))) {
                ++nfail;
                ADD_FAILURE() << "b_Xj[" << d << "," << k
                              << "] analytic=" << deriv.b_Xj(d, k)
                              << " cs=" << cs_b(d, k);
            }
        }
    }
    std::cout << "PanelLinSourceVelocity Xj: " << nfail << " failures\n";
}

TEST_F(SensitivityPingTest, AllInfluenceDeriv_xi) {
    int nfail = 0;
    // ConstSourceStream xi
    {
        auto d = PanelConstSourceStreamDeriv(kXj, kXi);
        VectorXd cs = readVector(*file_, "panels/ConstSourceStream/val_xi");
        for (int k = 0; k < 2; ++k) {
            if (!withinTol(d.val_xi(k), cs(k))) {
                ++nfail;
                ADD_FAILURE()
                    << "CSS val_xi[" << k << "] analytic=" << d.val_xi(k)
                    << " cs=" << cs(k);
            }
        }
    }
    // LinVortexVelocity xi
    {
        auto d = PanelLinVortexVelocityDeriv(kXj, kXi);
        MatrixXd cs_a = readMatrix(*file_, "panels/LinVortexVelocity/a_xi");
        MatrixXd cs_b = readMatrix(*file_, "panels/LinVortexVelocity/b_xi");
        for (int k = 0; k < 2; ++k) {
            for (int dd = 0; dd < 2; ++dd) {
                if (!withinTol(d.a_xi(dd, k), cs_a(dd, k))) {
                    ++nfail;
                    ADD_FAILURE() << "LVV a_xi[" << dd << "," << k << "]";
                }
                if (!withinTol(d.b_xi(dd, k), cs_b(dd, k))) {
                    ++nfail;
                    ADD_FAILURE() << "LVV b_xi[" << dd << "," << k << "]";
                }
            }
        }
    }
    // ConstSourceVelocity xi
    {
        auto d = PanelConstSourceVelocityDeriv(kXj, kXi);
        MatrixXd cs = readMatrix(*file_, "panels/ConstSourceVelocity/val_xi");
        for (int k = 0; k < 2; ++k) {
            for (int dd = 0; dd < 2; ++dd) {
                if (!withinTol(d.val_xi(dd, k), cs(dd, k))) {
                    ++nfail;
                    ADD_FAILURE() << "CSV val_xi[" << dd << "," << k << "]";
                }
            }
        }
    }
    // LinSourceStream xi
    {
        auto d = PanelLinSourceStreamDeriv(kXj, kXi);
        VectorXd cs_a = readVector(*file_, "panels/LinSourceStream/a_xi");
        VectorXd cs_b = readVector(*file_, "panels/LinSourceStream/b_xi");
        for (int k = 0; k < 2; ++k) {
            if (!withinTol(d.a_xi(k), cs_a(k))) {
                ++nfail;
                ADD_FAILURE() << "LSS a_xi[" << k << "]";
            }
            if (!withinTol(d.b_xi(k), cs_b(k))) {
                ++nfail;
                ADD_FAILURE() << "LSS b_xi[" << k << "]";
            }
        }
    }
    // LinSourceVelocity xi
    {
        auto d = PanelLinSourceVelocityDeriv(kXj, kXi);
        MatrixXd cs_a = readMatrix(*file_, "panels/LinSourceVelocity/a_xi");
        MatrixXd cs_b = readMatrix(*file_, "panels/LinSourceVelocity/b_xi");
        for (int k = 0; k < 2; ++k) {
            for (int dd = 0; dd < 2; ++dd) {
                if (!withinTol(d.a_xi(dd, k), cs_a(dd, k))) {
                    ++nfail;
                    ADD_FAILURE() << "LSV a_xi[" << dd << "," << k << "]";
                }
                if (!withinTol(d.b_xi(dd, k), cs_b(dd, k))) {
                    ++nfail;
                    ADD_FAILURE() << "LSV b_xi[" << dd << "," << k << "]";
                }
            }
        }
    }
    std::cout << "AllInfluenceDeriv_xi: " << nfail << " failures\n";
}

// ---------------------------------------------------------------------------
// Force partials: dcl/dx at fixed gamma, NACA 2412 α=5
// ---------------------------------------------------------------------------
TEST_F(SensitivityPingTest, ForcePartials_cl) {
    Cppfoil foil("2412", 80);
    foil.SetOper(5.0);
    foil.Solve();
    auto& state = foil.GetState();
    CalcForce(state);
    CalcForceLin(state);

    int N = state.foil.N;
    MatrixXd cs = readMatrix(*file_, "force_partial/dcl_dx_partial");
    ASSERT_EQ(cs.rows(), 2);
    ASSERT_EQ(cs.cols(), N);

    int nfail = 0;
    double max_abs = 0;
    for (int d = 0; d < 2; ++d) {
        for (int i = 0; i < N; ++i) {
            double got = state.post.dcl_dx_partial(d, i);
            double ref = cs(d, i);
            double err = std::abs(got - ref);
            if (err > max_abs) max_abs = err;
            if (!withinTol(got, ref)) {
                if (nfail < 5)
                    ADD_FAILURE() << "dcl_dx[" << d << "," << i
                                  << "] analytic=" << got << " cs=" << ref;
                ++nfail;
            }
        }
    }
    std::cout << "ForcePartials_cl: " << nfail
              << " failures, max|.|=" << max_abs << "\n";
}

// ---------------------------------------------------------------------------
// Inviscid Ma=0 and Re-inviscid sanity (no FD)
// ---------------------------------------------------------------------------
TEST_F(SensitivityPingTest, InviscidMachSensitivityZero) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, 0.0);
    foil.Solve();
    auto sens = foil.ComputeSensitivities();
    EXPECT_NEAR(sens.dcl_dMa, 0.0, 1e-12);
    EXPECT_NEAR(sens.dcm_dMa, 0.0, 1e-12);
    EXPECT_NEAR(sens.dcd_dMa, 0.0, 1e-12);
}

TEST_F(SensitivityPingTest, InviscidReSensitivityZero) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, 0.0);
    foil.Solve();
    auto sens = foil.ComputeSensitivities();
    EXPECT_NEAR(sens.dcl_dRe, 0.0, 1e-12);
    EXPECT_NEAR(sens.dcd_dRe, 0.0, 1e-12);
    EXPECT_NEAR(sens.dcm_dRe, 0.0, 1e-12);
}
