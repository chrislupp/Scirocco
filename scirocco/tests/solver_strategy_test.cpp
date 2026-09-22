#include <gtest/gtest.h>

#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/scirocco.h"

using namespace scirocco;

TEST(SolverStrategy, ScalingPreservesNewtonEquation) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 1e6, 0.2);
    foil.GetState().param.use_ptc = false;
    ASSERT_TRUE(foil.Solve());
    auto original = foil.GetState();
    original.glob.U.row(0) *= 1.03;
    BuildGlobSys(original);
    CalcForce(original);
    auto scaled = original;
    original.param.scale_coupled_system = false;
    scaled.param.scale_coupled_system = true;
    SolveGlob(original);
    SolveGlob(scaled);
    const VectorXd delta = Eigen::Map<const VectorXd>(scaled.glob.dU.data(),
                                                      scaled.glob.dU.size());
    const VectorXd defect =
        scaled.glob.R_V * delta + scaled.glob.coupled_residual;
    EXPECT_LT(defect.norm(), 1e-11 * scaled.glob.coupled_residual.norm());
    EXPECT_LT((original.glob.dU - scaled.glob.dU).norm(),
              1e-9 * original.glob.dU.norm());
    EXPECT_GT((scaled.glob.variable_scale.array() != 1.).count(), 0);
    EXPECT_GT((scaled.glob.residual_scale.array() != 1.).count(), 0);
}

TEST(SolverStrategy, TrustStepRejectsUphillWithoutChangingState) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    auto state = foil.GetState();
    state.param.use_trust_region = true;
    state.param.trust_max_trials = 1;
    state.glob.trust_radius = 1e6;
    state.glob.U.row(0) *= 1.001;
    BuildGlobSys(state);
    CalcForce(state);
    SolveGlob(state);
    state.glob.dU *= -1.0;
    const auto before = state;
    EXPECT_FALSE(TakeTrustRegionStep(state));
    EXPECT_EQ(state.glob.U, before.glob.U);
    EXPECT_EQ(state.isol.xi, before.isol.xi);
    EXPECT_EQ(state.vsol.turb, before.vsol.turb);
    EXPECT_EQ(state.post.cl, before.post.cl);
    EXPECT_EQ(state.oper.alpha, before.oper.alpha);
    EXPECT_EQ(state.glob.trust_rejected, before.glob.trust_rejected + 1);
}

TEST(SolverStrategy, TrustStepAcceptsLocalNewtonDirection) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    auto state = foil.GetState();
    state.param.use_trust_region = true;
    state.glob.trust_radius = 10.;
    state.glob.U.row(0) *= 1.001;
    BuildGlobSys(state);
    CalcForce(state);
    SolveGlob(state);
    EXPECT_TRUE(TakeTrustRegionStep(state));
    EXPECT_GT(state.glob.trust_last_ratio, 0.9);
    EXPECT_EQ(state.glob.trust_accepted, 1);
}

TEST(SolverStrategy, AcceptedTrialHasConsistentStagnationCoordinates) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    auto state = foil.GetState();
    state.param.use_trust_region = true;
    state.param.trust_update_topology = true;
    state.glob.trust_radius = 10.;
    state.glob.U.row(0) *= 1.01;
    BuildGlobSys(state);
    CalcForce(state);
    SolveGlob(state);
    ASSERT_TRUE(TakeTrustRegionStep(state));
    const int i = state.isol.Istag[0], j = state.isol.Istag[1];
    const double ui = state.glob.U(3, i), uj = state.glob.U(3, j);
    const double s = (uj * state.foil.s(i) + ui * state.foil.s(j)) / (ui + uj);
    EXPECT_NEAR(state.isol.sstag, s, 1e-14);
    for (int k = 0; k < state.foil.N; ++k)
        EXPECT_NEAR(state.isol.xi(k), std::abs(state.foil.s(k) - s), 1e-14);
    const auto accepted_r = state.glob.R;
    BuildGlobSys(state);
    EXPECT_LT((state.glob.R - accepted_r).norm(), 1e-14);
}
