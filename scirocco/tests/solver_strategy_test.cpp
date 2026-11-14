#include <gtest/gtest.h>

#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/scirocco.h"

using namespace scirocco;

namespace {

/// A converged state with the momentum thickness perturbed, so the Newton
/// step is nontrivial.
FoilState PerturbedState(double cl_target) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, cl_target, 1e6, 0.2);
    foil.GetState().param.use_ptc = false;
    EXPECT_TRUE(foil.Solve());
    FoilState state = foil.GetState();
    state.glob.U.row(0) *= 1.03;
    BuildGlobSys(state);
    CalcForce(state);
    return state;
}

VectorXd Packed(const FoilState& state) {
    VectorXd v(state.glob.dU.size() + (state.oper.givencl ? 1 : 0));
    v.head(state.glob.dU.size()) =
        Eigen::Map<const VectorXd>(state.glob.dU.data(), state.glob.dU.size());
    if (state.oper.givencl) v(v.size() - 1) = state.glob.dalpha;
    return v;
}

}  // namespace

TEST(SolverStrategy, ScalingPreservesNewtonEquation) {
    /* The assembled path (trust region) equilibrates its sparse LU. */
    auto original = PerturbedState(-1.0);
    original.param.use_trust_region = true;
    auto scaled = original;
    original.param.scale_coupled_system = false;
    scaled.param.scale_coupled_system = true;
    SolveGlob(original);
    SolveGlob(scaled);
    const VectorXd defect =
        scaled.glob.R_V * Packed(scaled) + scaled.glob.coupled_residual;
    EXPECT_LT(defect.norm(), 1e-11 * scaled.glob.coupled_residual.norm());
    EXPECT_LT((original.glob.dU - scaled.glob.dU).norm(),
              1e-9 * original.glob.dU.norm());
    EXPECT_GT((scaled.glob.variable_scale.array() != 1.).count(), 0);
    EXPECT_GT((scaled.glob.residual_scale.array() != 1.).count(), 0);
}

/* The Schur-complement step is the exact Newton step: it satisfies the
 * assembled Newton equation, at fixed alpha and with a Cl target. */
TEST(SolverStrategy, SchurStepSolvesAssembledNewtonEquation) {
    for (const double cl_target : {-1.0, 0.5}) {
        auto schur = PerturbedState(cl_target);
        auto assembled = schur;
        assembled.param.use_trust_region = true;
        SolveGlob(schur);
        SolveGlob(assembled);
        const VectorXd defect = assembled.glob.R_V * Packed(schur) +
                                assembled.glob.coupled_residual;
        EXPECT_LT(defect.norm(), 1e-11 * assembled.glob.coupled_residual.norm())
            << "cl_target " << cl_target;
        EXPECT_LT((Packed(schur) - Packed(assembled)).norm(),
                  1e-9 * Packed(assembled).norm())
            << "cl_target " << cl_target;
    }
}

/* Pseudo-time terms sit on each equation's own variable. With a tiny tau
 * the step reduces to dU = -tau R, node by node: boundary-layer equation k
 * of node i moves variable k of node i, and the edge-velocity equation moves
 * ue. A shifted diagonal would move other nodes' variables instead. */
TEST(SolverStrategy, PseudoTimeTermsActOnOwnVariables) {
    auto state = PerturbedState(-1.0);
    state.param.use_ptc = true;
    state.glob.tau = 1e-9;
    SolveGlob(state);
    const int n = state.glob.Nsys;
    const VectorXd& r = state.glob.coupled_residual;
    MatrixXd expected(4, n);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < 3; ++k)
            expected(k, i) = -state.glob.tau * r(3 * i + k);
        expected(3, i) = -r(3 * n + i) / (1.0 + 1.0 / state.glob.tau);
    }
    EXPECT_LT((state.glob.dU - expected).norm(), 1e-3 * expected.norm());
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

TEST(SolverStrategy, ContinuationPreservesTargetAndFrozenWake) {
    Scirocco foil("2412", 99);
    foil.SetOper(4.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    const auto target = foil.GetState();
    auto result = target;
    result.glob.conv = false;
    ASSERT_TRUE(SolveContinuation(target, result));
    EXPECT_TRUE(result.glob.conv);
    EXPECT_EQ(result.oper.alpha, target.oper.alpha);
    EXPECT_EQ(result.oper.Re, target.oper.Re);
    EXPECT_EQ(result.oper.Ma, target.oper.Ma);
    EXPECT_EQ(result.wake.x, target.wake.x);
    EXPECT_NEAR(result.post.cl, target.post.cl, 1e-10);
    EXPECT_NEAR(result.post.cd, target.post.cd, 1e-11);
    EXPECT_NEAR(result.post.cm, target.post.cm, 1e-10);
    EXPECT_GT(result.glob.continuation_steps, 0);
    EXPECT_LE(result.glob.continuation_iterations,
              target.param.continuation_max_iterations);
}

TEST(SolverStrategy, FailedContinuationRestoresRequestedPoint) {
    Scirocco foil("2412", 99);
    foil.SetOper(4.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    auto target = foil.GetState();
    target.param.continuation_max_iterations = 1;
    auto result = target;
    result.glob.conv = false;
    EXPECT_FALSE(SolveContinuation(target, result));
    EXPECT_FALSE(result.glob.conv);
    EXPECT_EQ(result.oper.alpha, target.oper.alpha);
    EXPECT_EQ(result.glob.U, target.glob.U);
    EXPECT_EQ(result.isol.xi, target.isol.xi);
    EXPECT_EQ(result.vsol.turb, target.vsol.turb);
    EXPECT_EQ(result.wake.x, target.wake.x);
    EXPECT_EQ(result.post.cl, target.post.cl);
    EXPECT_EQ(result.glob.continuation_iterations, 1);
}

TEST(SolverStrategy, ContinuationReachesExactLiftTarget) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, 0.6, 1e6, 0.2);
    ASSERT_TRUE(foil.Solve());
    const auto target = foil.GetState();
    auto result = target;
    result.glob.conv = false;
    ASSERT_TRUE(SolveContinuation(target, result));
    EXPECT_TRUE(result.oper.givencl);
    EXPECT_EQ(result.oper.cltgt, 0.6);
    EXPECT_NEAR(result.post.cl, 0.6, target.param.rtol);
    EXPECT_EQ(result.wake.x, target.wake.x);
    EXPECT_NEAR(result.oper.alpha, target.oper.alpha, 1e-9);
}
