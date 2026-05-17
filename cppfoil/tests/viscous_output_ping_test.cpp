/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cmath>

#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/util/math_utils.h"

using namespace cppfoil;
using namespace Eigen;

// ============================================================================
// Phase 3 ping: validate viscous output partials cl_V, cm_V, cd_V populated
// by CalcForceLin against finite differences inside cppfoil's own CalcForce.
//
// Self-contained — no HDF5/mfoil reference needed. The Squire-Young cd_U
// derivation in CalcForceLin is the only new analytic chain at this phase.
// cl_V/cm_V are aliases of cl_ue/cm_ue already validated by other tests but
// we sanity-check the structural mapping too.
// ============================================================================
class ViscousOutputPingTest : public ::testing::Test {
   protected:
    Cppfoil foil_{"0012", 100};

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
    // cl_V[4*i + 3] == cl_ue[i] for foil nodes; other entries zero.
    for (int i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 3), s.post.cl_ue(i));
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 0), 0.0);
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 1), 0.0);
        EXPECT_DOUBLE_EQ(s.post.cl_V(4 * i + 2), 0.0);
    }
    // Wake entries: all zero.
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

TEST_F(ViscousOutputPingTest, CdUMatchesFD) {
    auto& s = foil_.GetState();
    int Nsys = static_cast<int>(s.glob.U.cols());

    ASSERT_FALSE(s.vsol.Is.empty());
    ASSERT_FALSE(s.vsol.Is[2].empty());
    int iw = s.vsol.Is[2].back();

    // Reference cd from current state
    double cd0 = s.post.cd;

    // Forward-difference each Uw[k] and compare to cd_U[k].
    // Squire-Young is smooth and well-conditioned; h=1e-6 gives ~1e-7 FD noise.
    const double h = 1e-6;
    Vector4d cd_U_fd;
    for (int k = 0; k < 4; ++k) {
        FoilState s_perturbed = s;
        s_perturbed.glob.U(k, iw) += h;
        CalcForce(s_perturbed);
        cd_U_fd(k) = (s_perturbed.post.cd - cd0) / h;
    }

    // Hybrid tolerance: rel 1e-5, abs 1e-9 (FD noise floor).
    for (int k = 0; k < 4; ++k) {
        double got = s.post.cd_U(k);
        double ref = cd_U_fd(k);
        double err = std::abs(got - ref);
        double tol = 1e-5 * std::abs(ref) + 1e-9;
        EXPECT_LE(err, tol)
            << "cd_U[" << k << "]: got=" << got << " fd=" << ref
            << " abs=" << err;
    }

    // cd_V must be sparse at the end-of-wake column iw only.
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

// ============================================================================
// Phase 5: ∂R_V/∂alpha FD validation. Holds U and foil.x fixed, perturbs
// alpha, recomputes the global system, takes (R(α+h) - R(α-h))/(2h) and
// compares to analytic DRVDAlpha. Central FD at h=1e-6 → noise ~1e-12.
// ============================================================================
TEST_F(ViscousOutputPingTest, DRVDAlphaMatchesFD) {
    auto& s = foil_.GetState();
    int Nsys = s.glob.Nsys;

    // Analytical
    VectorXd dR_analytic = DRVDAlpha(s);

    // R_V at fixed (U, foil.x, wake.x, ue_sigma); only alpha varies through
    // gamma → ueinv. This matches DRVDAlpha's "partial holding U and the
    // geometry-derived matrices fixed" definition, and is what the Newton
    // solve sees in its constrained-cl block.
    int N = s.foil.N;
    VectorXd ds_state(Nsys);
    VectorXd ue_state(Nsys);
    for (int j = 0; j < Nsys; ++j) {
        ds_state(j) = s.glob.U(1, j);
        ue_state(j) = s.glob.U(3, j);
    }
    VectorXd mass = ds_state.cwiseProduct(ue_state);
    VectorXd source_term = s.vsol.ue_m * mass;  // sgue * ue_sigma * phi
    // R_BL doesn't depend on alpha → that block is zero. Only assemble R_ue.
    auto eval_Rue = [&](double alpha) -> VectorXd {
        double ca = Cosd(alpha);
        double sa = Sind(alpha);
        VectorXd gam(N);
        for (int i = 0; i < N; ++i) {
            gam(i) = s.isol.gamref(i, 0) * ca + s.isol.gamref(i, 1) * sa;
        }
        VectorXd uewi(Nsys - N);
        for (int j = 0; j < Nsys - N; ++j) {
            uewi(j) = s.isol.uewiref(j, 0) * ca + s.isol.uewiref(j, 1) * sa;
        }
        VectorXd R_ue(Nsys);
        for (int i = 0; i < Nsys; ++i) {
            double ueinv_i;
            if (i < N) {
                ueinv_i = s.isol.sgnue(i) * gam(i);
            } else if (i - N == 0) {
                ueinv_i = s.isol.sgnue(N - 1) * gam(N - 1);
            } else {
                ueinv_i = uewi(i - N);
            }
            R_ue(i) = ue_state(i) - ueinv_i - source_term(i);
        }
        return R_ue;
    };

    double alpha0 = s.oper.alpha;
    double h = 1e-6;
    VectorXd Rup = eval_Rue(alpha0 + h);
    VectorXd Rum = eval_Rue(alpha0 - h);
    VectorXd dRue_fd = (Rup - Rum) / (2.0 * h);

    // Analytic R_ue block lives in rows 3*Nsys..4*Nsys-1.
    double max_err = 0.0;
    int worst_i = 0;
    for (int i = 0; i < Nsys; ++i) {
        double e = std::abs(dR_analytic(3 * Nsys + i) - dRue_fd(i));
        if (e > max_err) { max_err = e; worst_i = i; }
    }
    std::cout << "DRVDAlpha: max abs vs FD = " << max_err
              << " at R_ue row " << worst_i << std::endl;
    // FD noise floor ~1e-9 at central-difference h=1e-6 on smooth R_ue.
    EXPECT_LE(max_err, 1e-9);

    // R_BL block (rows 0..3*Nsys-1) must be exactly zero analytically.
    EXPECT_EQ(dR_analytic.head(3 * Nsys).cwiseAbs().maxCoeff(), 0.0);
}

// ============================================================================
// Phase 4: forward/reverse adjoint consistency.
//   For any b, g ∈ ℝ^{4·Nsys}:
//     bᵀ · (R_V^{-1} g) == (R_V^{-T} b)ᵀ · g
//   Pure linear-algebra identity. Tests the transpose solve through
//   SolveViscousAdjoint vs. SparseLU's forward solve.
// ============================================================================
TEST_F(ViscousOutputPingTest, ForwardReverseConsistency) {
    auto& s = foil_.GetState();
    StoreAdjointFactorization(s);
    ASSERT_TRUE(s.glob.lu_factorization) << "R_V LU factorization unavailable";

    int Nsys = s.glob.Nsys;
    int NN = 4 * Nsys;
    VectorXd b = VectorXd::Random(NN);
    VectorXd g = VectorXd::Random(NN);

    VectorXd y = s.glob.lu_factorization->solve(g);           // R_V^{-1} g
    VectorXd lambda = SolveViscousAdjoint(s, b);              // R_V^{-T} b

    double lhs = b.dot(y);
    double rhs = lambda.dot(g);
    double abs_err = std::abs(lhs - rhs);
    double rel_err = abs_err / std::max(std::abs(lhs), 1.0);

    // Linear-algebra identity should hold near machine precision.
    EXPECT_LE(rel_err, 1e-12)
        << "Forward/reverse adjoint inconsistency: lhs=" << lhs
        << " rhs=" << rhs << " abs=" << abs_err << " rel=" << rel_err;
}
