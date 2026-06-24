/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/

#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/wake.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
using namespace Eigen;

/**
 * @brief Phase 3 ping: validate viscous output partials cl_V, cm_V, cd_V
 * populated by CalcForceLin against finite differences inside scirocco's own
 * CalcForce.
 *
 * Self-contained — no HDF5/mfoil reference needed. The Squire-Young cd_U
 * derivation in CalcForceLin is the only new analytic chain at this phase.
 * cl_V/cm_V are aliases of cl_ue/cm_ue already validated by other tests but
 * we sanity-check the structural mapping too.
 */
class ViscousOutputPingTest : public ::testing::Test {
   protected:
    Scirocco foil_{"0012", 100};

    void SetUp() override {
        foil_.SetOper(3.0, -1.0, 1e6, 0.0);
        auto& s = foil_.GetState();
        s.oper.viscous = true;
        s.oper.initbl = true;
        s.param.verb = 0;
        SolveViscous(s);
        ASSERT_TRUE(s.glob.conv) << "viscous solve did not converge";
        CalcForce(s);
        CalcForceLin(s);
    }
};

TEST_F(ViscousOutputPingTest, ClVMatchesClUe) {
    auto& s = foil_.GetState();
    int N = s.foil.N;
    int Nsys = static_cast<int>(s.glob.U.cols());
    ASSERT_EQ(s.post.cl_V.size(), 4 * Nsys);
    /* cl_V[4*i + 3] == cl_ue[i] for foil nodes; other entries zero. */
    for (int i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 3), s.post.cl_ue(i));
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 0), 0.0);
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 1), 0.0);
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 2), 0.0);
    }
    /* Wake entries: all zero. */
    for (int i = N; i < Nsys; ++i) {
        for (int k = 0; k < 4; ++k) {
            EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + k), 0.0);
        }
    }
}

TEST_F(ViscousOutputPingTest, CmVMatchesCmUe) {
    auto& s = foil_.GetState();
    int N = s.foil.N;
    int Nsys = static_cast<int>(s.glob.U.cols());
    ASSERT_EQ(s.post.cm_V.size(), 4 * Nsys);
    for (int i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(s.post.cm_V(4 * i + 3), s.post.cm_ue(i));
    }
    for (int i = N; i < Nsys; ++i) {
        for (int k = 0; k < 4; ++k) {
            EXPECT_DOUBLE_EQ(s.post.cm_V(4 * i + k), 0.0);
        }
    }
}

/**
 * @brief CS reference: ∂cd/∂Uw at iw computed via mfoil_cs complex-step on
 * calc_force at perturbed glob.U[:, iw]. Reference file:
 * extract_viscous_output_ping_cs.py.
 */
TEST_F(ViscousOutputPingTest, CdUMatchesCS) {
    auto& s = foil_.GetState();
    int Nsys = static_cast<int>(s.glob.U.cols());

    ASSERT_FALSE(s.vsol.Is.empty());
    ASSERT_FALSE(s.vsol.Is[2].empty());
    int iw = s.vsol.Is[2].back();

    /* Load CS reference */
    H5::H5File file(VISCOUS_OUTPUT_PING_REFERENCE_FILE, H5F_ACC_RDONLY);
    int iw_ref = 0;
    file.openDataSet("iw").read(&iw_ref, H5::PredType::NATIVE_INT);
    ASSERT_EQ(iw, iw_ref);
    std::vector<double> cd_U_cs(4);
    file.openDataSet("cd_U").read(cd_U_cs.data(), H5::PredType::NATIVE_DOUBLE);

    /* Compare analytic to mfoil-CS at machine precision. */
    for (int k = 0; k < 4; ++k) {
        double got = s.post.cd_U(k);
        double ref = cd_U_cs[k];
        double err = std::abs(got - ref);
        double tol = 1e-9 * std::abs(ref) + 1e-12;
        EXPECT_LE(err, tol) << "cd_U[" << k << "]: got=" << got << " cs=" << ref
                            << " abs=" << err;
    }

    /* cd_V must be sparse at the end-of-wake column iw only. */
    ASSERT_EQ(s.post.cd_V.size(), 4 * Nsys);
    for (int i = 0; i < Nsys; ++i) {
        for (int k = 0; k < 4; ++k) {
            double v = s.post.cd_V(4 * i + k);
            if (i == iw) {
                EXPECT_DOUBLE_EQ(v, s.post.cd_U(k));
            } else {
                EXPECT_DOUBLE_EQ(v, 0.0);
            }
        }
    }
}

TEST_F(ViscousOutputPingTest, DcdDxPartialIsZero) {
    auto& s = foil_.GetState();
    int N = s.foil.N;
    ASSERT_EQ(s.post.dcd_dx_partial.rows(), 2);
    ASSERT_EQ(s.post.dcd_dx_partial.cols(), N);
    EXPECT_EQ(s.post.dcd_dx_partial.cwiseAbs().maxCoeff(), 0.0)
        << "Squire-Young cd depends only on Uw, no direct foil.x dependence "
           "→ dcd/dx_partial must be exactly zero.";
}

/**
 * @brief Phase 5: ∂R_V/∂alpha FD validation. Holds U and foil.x fixed,
 * perturbs alpha, recomputes the global system, takes (R(α+h) - R(α-h))/(2h)
 * and compares to analytic DRVDAlpha. Central FD at h=1e-6 → noise ~1e-12.
 *
 * CS reference: ∂R_V/∂alpha at fixed glob.U via mfoil_cs complex-step.
 * Reference matches scirocco's ue-based convention (StagPointMove with
 * fixed glob.U[3]).
 */
TEST_F(ViscousOutputPingTest, DRVDAlphaMatchesCS) {
    auto& s = foil_.GetState();
    int Nsys = s.glob.Nsys;

    VectorXd dR_analytic = DRVDAlpha(s);

    H5::H5File file(VISCOUS_OUTPUT_PING_REFERENCE_FILE, H5F_ACC_RDONLY);
    int Nsys_ref = 0;
    file.openDataSet("Nsys").read(&Nsys_ref, H5::PredType::NATIVE_INT);
    ASSERT_EQ(Nsys, Nsys_ref);
    std::vector<double> dRV_dalpha_cs(4 * Nsys);
    file.openDataSet("dRV_dalpha")
        .read(dRV_dalpha_cs.data(), H5::PredType::NATIVE_DOUBLE);
    VectorXd dR_cs = Eigen::Map<VectorXd>(dRV_dalpha_cs.data(), 4 * Nsys);

    double max_abs = (dR_analytic - dR_cs).cwiseAbs().maxCoeff();
    int worst_i = 0;
    (dR_analytic - dR_cs).cwiseAbs().maxCoeff(&worst_i);
    std::cout << "DRVDAlpha: max abs vs CS = " << max_abs << " at row "
              << worst_i << " (analytic=" << dR_analytic(worst_i)
              << " cs=" << dR_cs(worst_i) << ")" << std::endl;
    double max_rel = 0;
    for (int i = 0; i < 4 * Nsys; ++i) {
        double mag = std::max(std::abs(dR_analytic(i)), std::abs(dR_cs(i)));
        if (mag < 1e-6) continue;
        double rel = std::abs(dR_analytic(i) - dR_cs(i)) / mag;
        if (rel > max_rel) max_rel = rel;
    }
    std::cout << "DRVDAlpha: max rel (|val|>1e-6) = " << max_rel << std::endl;
    EXPECT_LE(max_abs, 1e-9);
    EXPECT_LE(max_rel, 1e-9);
}

/**
 * @brief Phase 4: forward/reverse adjoint consistency.
 *   For any b, g ∈ ℝ^{4·Nsys}:
 *     bᵀ · (R_V^{-1} g) == (R_V^{-T} b)ᵀ · g
 *   Pure linear-algebra identity. Tests the transpose solve through
 *   SolveViscousAdjoint vs. SparseLU's forward solve.
 */
TEST_F(ViscousOutputPingTest, ForwardReverseConsistency) {
    auto& s = foil_.GetState();
    StoreAdjointFactorization(s);
    ASSERT_TRUE(s.glob.lu_factorization) << "R_V LU factorization unavailable";

    int Nsys = s.glob.Nsys;
    int NN = 4 * Nsys;
    VectorXd b = VectorXd::Random(NN);
    VectorXd g = VectorXd::Random(NN);

    VectorXd y = s.glob.lu_factorization->solve(g);  /* R_V^{-1} g */
    VectorXd lambda = SolveViscousAdjoint(s, b);     /* R_V^{-T} b */

    double lhs = b.dot(y);
    double rhs = lambda.dot(g);
    double abs_err = std::abs(lhs - rhs);
    double rel_err = abs_err / std::max(std::abs(lhs), 1.0);

    /* Linear-algebra identity should hold near machine precision. */
    EXPECT_LE(rel_err, 1e-12)
        << "Forward/reverse adjoint inconsistency: lhs=" << lhs
        << " rhs=" << rhs << " abs=" << abs_err << " rel=" << rel_err;
}
