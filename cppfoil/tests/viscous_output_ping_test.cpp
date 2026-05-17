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
