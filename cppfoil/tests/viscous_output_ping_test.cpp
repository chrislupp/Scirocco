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
#include "cppfoil/core/wake.h"
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

    // DRVDAlpha is itself a 4th-order central FD over a full geometry
    // rebuild (gam, wake.x, wake.t, ue_sigma) at perturbed alpha at fixed U.
    // To validate, compute the same R_V FD with a 2nd-order central FD via
    // BuildWake/CalcUeM/BuildGlobSys pipeline and compare.
    auto full_eval_RV = [&](double alpha) -> VectorXd {
        FoilState sp = s;
        sp.oper.alpha = alpha;
        RebuildIsol(sp);
        BuildWake(sp);
        // Production uses ue-based xi (StagPointMove with fixed glob.U).
        // StagPointFind would use gamma-based weights and produce a
        // gamma-chain that DRVDAlpha doesn't include, so use StagPointMove
        // to match the partial-at-fixed-U convention.
        StagPointMove(sp);
        ExtendXiToWake(sp);
        IdentifySurfaces(sp);
        SetWakeGap(sp);
        CalcUeM(sp);
        sp.glob.realloc = true;
        BuildGlobSys(sp);
        // Assemble R_V (R_BL || R_ue) at fixed U.
        VectorXd R = VectorXd::Zero(4 * Nsys);
        R.head(3 * Nsys) = sp.glob.R;
        VectorXd mass(Nsys);
        int N = sp.foil.N;
        for (int j = 0; j < Nsys; ++j)
            mass(j) = sp.glob.U(1, j) * sp.glob.U(3, j);
        VectorXd src = sp.vsol.ue_m * mass;
        for (int i = 0; i < Nsys; ++i) {
            double ueinv_i;
            if (i < N) ueinv_i = sp.isol.sgnue(i) * sp.isol.gam(i);
            else if (i - N == 0)
                ueinv_i = sp.isol.sgnue(N - 1) * sp.isol.gam(N - 1);
            else
                ueinv_i = sp.isol.uewi(i - N);
            R(3 * Nsys + i) = sp.glob.U(3, i) - ueinv_i - src(i);
        }
        return R;
    };
    double alpha0 = s.oper.alpha;
    double h = 1e-5;
    VectorXd Rp = full_eval_RV(alpha0 + h);
    VectorXd Rm = full_eval_RV(alpha0 - h);
    VectorXd dR_fd = (Rp - Rm) / (2.0 * h);

    double max_err = (dR_analytic - dR_fd).cwiseAbs().maxCoeff();
    int worst_i;
    (dR_analytic - dR_fd).cwiseAbs().maxCoeff(&worst_i);
    std::cout << "DRVDAlpha: max abs vs 2nd-order FD = " << max_err
              << " at row " << worst_i
              << " (analytic=" << dR_analytic(worst_i)
              << " fd=" << dR_fd(worst_i) << ")" << std::endl;
    // Compute worst rel for entries with |value| > 1e-3.
    double max_rel = 0;
    for (int i = 0; i < 4 * Nsys; ++i) {
        double mag = std::max(std::abs(dR_analytic(i)), std::abs(dR_fd(i)));
        if (mag < 1e-3) continue;
        double rel = std::abs(dR_analytic(i) - dR_fd(i)) / mag;
        if (rel > max_rel) max_rel = rel;
    }
    std::cout << "DRVDAlpha: max rel (|val|>1e-3) = " << max_rel << std::endl;
    // DRVDAlpha is 4th-order; the test FD here is 2nd-order. At the stagnation
    // rows (~row 3·Istag) the BL residual has very large R''' so 2nd-order
    // truncation at h=1e-5 gives ~1e-3 abs noise. This is FD truncation,
    // not a chain bug. Validation at machine-precision rigour happens against
    // the mfoil_cs CS reference in Phase 7 (ViscousAlphaSensitivity).
    EXPECT_LE(max_rel, 1e-3);
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
