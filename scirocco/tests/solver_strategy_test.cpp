#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include "scirocco/core/bl_march.h"
#include "scirocco/core/coupled_linear.h"
#include "scirocco/core/coupled_solver.h"
#include "scirocco/core/forces.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/paneling.h"
#include "scirocco/core/transition.h"
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

/// A small coupled system: each station's boundary-layer rows depend on its
/// own variables and the station upstream, every ue on every mass.
struct SmallSystem {
    SparseMatrixd bl;
    MatrixXd ue_m;
    VectorXd ue, ds, d;
};

/// @param empty_station Station left without boundary-layer equations
SmallSystem MakeSmallSystem(int n, int empty_station = -1) {
    SmallSystem sys;
    std::vector<Eigen::Triplet<double>> t;
    for (int i = 0; i < n; ++i)
        for (int r = 0; r < 3 && i != empty_station; ++r)
            for (int k = 0; k < 4; ++k) {
                const double own = (r == k ? 4.0 : 0.0) + 0.1 * (r + 2 * k + i);
                t.emplace_back(3 * i + r, 4 * i + k, own);
                if (i > 0)
                    t.emplace_back(3 * i + r, 4 * (i - 1) + k,
                                   -0.3 + 0.05 * (r - k));
            }
    sys.bl.resize(3 * n, 4 * n);
    sys.bl.setFromTriplets(t.begin(), t.end());
    sys.ue_m.resize(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            sys.ue_m(i, j) = 0.2 / (1.0 + std::abs(i - j));
    sys.ue = VectorXd::LinSpaced(n, 0.8, 1.2);
    sys.ds = VectorXd::LinSpaced(n, 0.002, 0.01);
    sys.d = VectorXd::Ones(n);
    return sys;
}

/// Solve and transposed solve against the assembled matrix.
void ExpectSolvesAssembled(CoupledFactorization& f, int n) {
    const MatrixXd j = MatrixXd(f.Assemble());
    const VectorXd rhs = VectorXd::LinSpaced(4 * n, -1.0, 2.0);
    EXPECT_LT((j * f.Solve(rhs) - rhs).norm(), 1e-12 * rhs.norm());
    EXPECT_LT((j.transpose() * f.SolveTranspose(rhs) - rhs).norm(),
              1e-12 * rhs.norm());
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

/* NACA 0012 at alpha 5, 199 panels: the first Newton step used to flip the
 * whole upper surface laminar, and the step limiter then held the step to a
 * few percent for some 40 iterations. The boundary-layer march keeps the
 * boundary layer and transition consistent with ue between Newton steps, so
 * Newton converges directly, to the same solution. */
TEST(SolverStrategy, BoundaryLayerMarchConvergesWithoutStalling) {
    Scirocco marched("0012", 199);
    marched.SetOper(5.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(marched.Solve());
    EXPECT_LE(marched.GetState().glob.coupled_iterations, 20);
    EXPECT_EQ(marched.GetState().glob.continuation_attempts, 0);

    Scirocco plain("0012", 199);
    plain.SetOper(5.0, -1.0, 1e6, 0.0);
    plain.GetState().param.bl_march = false;
    plain.GetState().param.step_limiter = StepLimiter::Absolute;
    ASSERT_TRUE(plain.Solve());
    EXPECT_GT(plain.GetState().glob.coupled_iterations, 30);
    EXPECT_NEAR(marched.GetCl(), plain.GetCl(), 1e-8);
    EXPECT_NEAR(marched.GetCd(), plain.GetCd(), 1e-9);
}

/* The march solves the global residual equations, and places transition
 * with UpdateTransition's test, so a converged state is a fixed point. */
TEST(SolverStrategy, BoundaryLayerMarchLeavesConvergedStateUnchanged) {
    Scirocco foil("2412", 99);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    FoilState state = foil.GetState();
    EXPECT_EQ(BoundaryLayerMarch(state), 0);
    EXPECT_LT((state.glob.U - foil.GetState().glob.U).norm(),
              1e-9 * state.glob.U.norm());
    EXPECT_EQ(state.vsol.turb, foil.GetState().vsol.turb);
}

/* XFOIL's rule: a change in theta of +300% is scaled to +150%, a change of
 * -90% to -50%. The mfoil rule leaves theta increases unlimited. */
TEST(SolverStrategy, RelativeLimiterBoundsNormalizedChanges) {
    auto state = PerturbedState(-1.0);
    const int k = state.glob.Nsys / 3;
    for (const double change : {3.0, -0.9}) {
        FoilState trial = state;
        trial.param.step_limiter = StepLimiter::Relative;
        trial.glob.dU.setZero();
        trial.glob.dU(0, k) = change * trial.glob.U(0, k);
        const double before = trial.glob.U(0, k);
        UpdateState(trial);
        const double expected = change > 0.0 ? 1.5 : -0.5;
        EXPECT_NEAR(trial.glob.U(0, k) / before - 1.0, expected, 1e-12);
    }
    FoilState absolute = state;
    absolute.param.step_limiter = StepLimiter::Absolute;
    absolute.glob.dU.setZero();
    absolute.glob.dU(0, k) = 3.0 * absolute.glob.U(0, k);
    const double before = absolute.glob.U(0, k);
    UpdateState(absolute);
    EXPECT_NEAR(absolute.glob.U(0, k) / before, 4.0, 1e-12);
}

/* NACA 0012 at alpha 0 has its stagnation point on a node; the residual
 * stalls near 1e-10, above rtol. Newton used to run out its 200 iterations
 * and fall back to continuation. The round-off floor test ends it. (cl is
 * -4.5e-8, not 0: the stagnation node breaks the symmetry.) */
TEST(SolverStrategy, RoundOffFloorEndsNewton) {
    Scirocco foil("0012", 130);
    foil.SetOper(0.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_LT(foil.GetState().glob.coupled_iterations, 30);
    EXPECT_EQ(foil.GetState().glob.continuation_attempts, 0);
    EXPECT_NEAR(foil.GetCl(), 0.0, 1e-7);
}

/* A second solve on the same object starts from the first solution: fewer
 * iterations, same answer as a cold start. */
TEST(SolverStrategy, WarmStartReusesPreviousSolution) {
    Scirocco warm("0012", 200);
    warm.GetState().param.warm_start = true;
    warm.SetOper(8.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(warm.Solve());
    warm.SetOper(8.5, -1.0, 1e6, 0.0);
    ASSERT_TRUE(warm.Solve());

    Scirocco cold("0012", 200);
    cold.SetOper(8.5, -1.0, 1e6, 0.0);
    cold.GetState().param.warm_start = false;
    ASSERT_TRUE(cold.Solve());
    EXPECT_LT(warm.GetState().glob.coupled_iterations,
              cold.GetState().glob.coupled_iterations);
    EXPECT_NEAR(warm.GetCl(), cold.GetCl(), 1e-8);
    EXPECT_NEAR(warm.GetCd(), cold.GetCd(), 1e-9);
}

/* With the default settings a solve does not depend on earlier solves on
 * the same object. NACA 0012 at Re 2e5 and 14.25 deg has two converged
 * boundary layers (Cl = 0.71 and 0.93), and a warm sweep may land on either
 * one: which depends on round-off, so it varies between builds. */
TEST(SolverStrategy, DefaultSolveIgnoresEarlierSolves) {
    EXPECT_FALSE(Param{}.warm_start);
    Scirocco fresh("0012", 160);
    fresh.SetOper(14.25, -1.0, 2e5, 0.0);
    ASSERT_TRUE(fresh.Solve());

    Scirocco reused("0012", 160);
    for (const double alpha : {13.5, 13.75, 14.0, 14.25}) {
        reused.SetOper(alpha, -1.0, 2e5, 0.0);
        ASSERT_TRUE(reused.Solve());
    }
    EXPECT_NEAR(reused.GetCl(), fresh.GetCl(), 1e-10);
    EXPECT_NEAR(reused.GetCd(), fresh.GetCd(), 1e-10);
}

/* A cold solve on a reused object follows the same Newton path as one on a
 * new object: nothing from the earlier solve (such as the stagnation-point
 * linearization isol.sstag_ue) reaches its first Jacobian. */
TEST(SolverStrategy, ColdSolveIgnoresEarlierSolves) {
    Scirocco fresh("0012", 160);
    fresh.GetState().param.warm_start = false;
    fresh.SetOper(-4.0, -1.0, 2e5, 0.0);
    ASSERT_TRUE(fresh.Solve());

    for (const double before : {-4.0, 4.0, 12.0}) {
        Scirocco reused("0012", 160);
        reused.GetState().param.warm_start = false;
        reused.SetOper(before, -1.0, 2e5, 0.0);
        ASSERT_TRUE(reused.Solve());
        reused.SetOper(-4.0, -1.0, 2e5, 0.0);
        ASSERT_TRUE(reused.Solve());
        EXPECT_EQ(reused.GetState().glob.coupled_iterations,
                  fresh.GetState().glob.coupled_iterations)
            << "after a solve at alpha " << before;
        EXPECT_EQ(reused.GetCl(), fresh.GetCl()) << "alpha " << before;
        EXPECT_EQ(reused.GetCd(), fresh.GetCd()) << "alpha " << before;
    }
}

/* The influence matrix and its factorization are reused while the panel
 * nodes are unchanged, and rebuilt as soon as any node moves. */
TEST(SolverStrategy, InviscidCacheFollowsGeometry) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    auto& state = foil.GetState();
    const auto lu = state.isol.AIC_lu;
    const MatrixXd bp = state.isol.Bp_foil;
    foil.SetOper(4.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_EQ(state.isol.AIC_lu, lu);
    EXPECT_EQ(state.isol.Bp_foil, bp);

    state.foil.x(1, 20) += 1e-4;
    RefreshPanelGeometry(state.foil);
    FoilState fresh = state;
    fresh.isol.AIC_x.resize(0, 0);
    fresh.isol.Bp_foil_x.resize(0, 0);
    BuildGamma(state, 4.0);
    BuildGamma(fresh, 4.0);
    EXPECT_NE(state.isol.AIC_lu, lu);
    EXPECT_EQ(state.isol.gamref, fresh.isol.gamref);
}

/* Newton iterations skip dR/dRe and dR/dMa; the converged state must still
 * carry them, exactly as a full assembly gives them. */
TEST(SolverStrategy, ConvergedStateCarriesScalarDerivatives) {
    Scirocco foil("2412", 99);
    foil.SetOper(3.0, -1.0, 1e6, 0.3);
    ASSERT_TRUE(foil.Solve());
    FoilState state = foil.GetState();
    const VectorXd re = state.glob.dR_dRe_times_Re;
    const VectorXd ma = state.glob.dR_dMa_times_Ma;
    ASSERT_GT(re.norm(), 0.0);
    ASSERT_GT(ma.norm(), 0.0);
    BuildGlobSys(state);
    EXPECT_EQ(re, state.glob.dR_dRe_times_Re);
    EXPECT_EQ(ma, state.glob.dR_dMa_times_Ma);
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

/* The Schur and direct paths solve the same system. A zero edge velocity
 * makes the mass transform singular, so it must fall back to the direct
 * factorization. */
TEST(CoupledFactorization, SchurAndDirectPathsSolveAssembledSystem) {
    const int n = 5;
    SmallSystem sys = MakeSmallSystem(n);
    CoupledFactorization schur;
    ASSERT_TRUE(schur.Compute(sys.bl, sys.ue_m, sys.ue, sys.ds, sys.d));
    EXPECT_TRUE(schur.UsesSchur());
    ExpectSolvesAssembled(schur, n);

    sys.ue(2) = 0.0;
    CoupledFactorization direct;
    ASSERT_TRUE(direct.Compute(sys.bl, sys.ue_m, sys.ue, sys.ds, sys.d));
    EXPECT_FALSE(direct.UsesSchur());
    ExpectSolvesAssembled(direct, n);
}

/* A station with no boundary-layer equations makes J singular: both
 * factorizations fail, and solving reports it instead of returning junk. */
TEST(CoupledFactorization, SingularSystemFailsLoudly) {
    const int n = 4;
    const SmallSystem sys = MakeSmallSystem(n, 1);
    CoupledFactorization f;
    EXPECT_FALSE(f.Compute(sys.bl, sys.ue_m, sys.ue, sys.ds, sys.d));
    const VectorXd rhs = VectorXd::Ones(4 * n);
    EXPECT_THROW((void)f.Solve(rhs), std::runtime_error);
    EXPECT_THROW((void)f.SolveTranspose(rhs), std::runtime_error);
}

TEST(CoupledFactorization, RejectsInconsistentSizes) {
    const int n = 4;
    const SmallSystem sys = MakeSmallSystem(n);
    CoupledFactorization f;
    EXPECT_THROW(f.Compute(sys.bl, sys.ue_m, sys.ue, sys.ds.head(n - 1), sys.d),
                 std::invalid_argument);
    ASSERT_TRUE(f.Compute(sys.bl, sys.ue_m, sys.ue, sys.ds, sys.d));
    const VectorXd short_rhs = VectorXd::Ones(4 * n - 1);
    EXPECT_THROW((void)f.Solve(short_rhs), std::invalid_argument);
    EXPECT_THROW((void)f.SolveTranspose(short_rhs), std::invalid_argument);
}

/* Without a cached LU, SolveAIC factors the influence matrix itself. */
TEST(SolverStrategy, SolveAICWithoutCachedFactorization) {
    Scirocco foil("2412", 99);
    foil.SetOper(2.0, -1.0, 0.0, 0.0);
    ASSERT_TRUE(foil.Solve());
    FoilState state = foil.GetState();
    const MatrixXd rhs = MatrixXd::Identity(state.isol.AIC.rows(), 3);
    const MatrixXd cached = SolveAIC(state, rhs);
    state.isol.AIC_lu.reset();
    EXPECT_LT((SolveAIC(state, rhs) - cached).norm(), 1e-12 * cached.norm());
}

/* Stations whose local solve fails are extrapolated from upstream or kept,
 * so a bad edge velocity never leaves a non-physical boundary layer for
 * the Newton step. */
TEST(SolverStrategy, BoundaryLayerMarchSurvivesFailedStations) {
    Scirocco foil("2412", 99);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    FoilState state = foil.GetState();
    for (int j = 0; j < state.glob.Nsys; j += 2) state.glob.U(3, j) *= 3.0;
    EXPECT_GT(BoundaryLayerMarch(state), 0);
    EXPECT_TRUE(state.glob.U.allFinite());
    EXPECT_GT(state.glob.U.row(0).minCoeff(), 0.0);
    EXPECT_GT(state.glob.U.row(1).minCoeff(), 0.0);
}

/* With a Cl target the step borders the factorization with the lift
 * constraint; a singular constraint must be reported, not stepped on. */
TEST(SolverStrategy, SingularLiftConstraintThrows) {
    auto state = PerturbedState(0.5);
    state.post.cl_alpha = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(SolveGlob(state), std::runtime_error);
}

/* A warm start that does not converge falls back to a cold start, which
 * must give exactly what a fresh object gives. */
TEST(SolverStrategy, FailedWarmStartMatchesColdStart) {
    Scirocco warm("0012", 99);
    warm.GetState().param.warm_start = true;
    warm.SetOper(2.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(warm.Solve());
    warm.SetOper(16.0, -1.0, 1e6, 0.0);
    const bool warm_ok = warm.Solve();

    Scirocco cold("0012", 99);
    cold.SetOper(16.0, -1.0, 1e6, 0.0);
    EXPECT_EQ(warm_ok, cold.Solve());
    if (warm_ok) EXPECT_NEAR(warm.GetCl(), cold.GetCl(), 1e-6);
}

/// The observer is called once per Newton iteration and leaves the solve
/// unchanged.
TEST(SolverStrategy, IterationObserver) {
    Scirocco plain("2412", 99);
    plain.SetOper(3.0, -1.0, 1e6, 0.2);
    ASSERT_TRUE(plain.Solve());

    Scirocco watched("2412", 99);
    watched.SetOper(3.0, -1.0, 1e6, 0.2);
    std::vector<int> iterations;
    std::vector<double> relax;
    watched.GetState().glob.iteration_observer = [&](const FoilState& s,
                                                     int iter) {
        iterations.push_back(iter);
        relax.push_back(s.glob.relax);
    };
    ASSERT_TRUE(watched.Solve());

    ASSERT_FALSE(iterations.empty());
    EXPECT_EQ(static_cast<int>(iterations.size()),
              watched.GetState().glob.coupled_iterations);
    for (std::size_t k = 0; k < iterations.size(); ++k) {
        EXPECT_EQ(iterations[k], static_cast<int>(k));
        EXPECT_GT(relax[k], 0.0);
        EXPECT_LE(relax[k], 1.0);
    }
    EXPECT_EQ(relax.front(), 1.0);
    EXPECT_EQ(watched.GetCl(), plain.GetCl());
    EXPECT_EQ(watched.GetCd(), plain.GetCd());
}

namespace {

/// NACA 2412 at Re 1e6 with trips on the upper and lower surfaces.
Scirocco TrippedFoil(double alpha, double xtrip_upper, double xtrip_lower) {
    Scirocco foil("2412", 160);
    foil.SetOper(alpha, -1.0, 1e6, 0.0);
    foil.GetState().param.xtrip_upper = xtrip_upper;
    foil.GetState().param.xtrip_lower = xtrip_lower;
    return foil;
}

}  // namespace

/* A trip ahead of free transition moves transition to the trip and adds
 * drag; the lower surface, laminar to the trailing edge when free, turns
 * turbulent at its trip. */
TEST(ForcedTransition, TripSetsTransitionLocation) {
    Scirocco free_foil = TrippedFoil(4.0, 1.0, 1.0);
    ASSERT_TRUE(free_foil.Solve());
    const auto& f = free_foil.GetState();
    EXPECT_FALSE(f.vsol.tripped[0]);
    EXPECT_FALSE(f.vsol.tripped[1]);
    EXPECT_FALSE(f.vsol.turb[f.vsol.Is[0].back()]);
    EXPECT_GT(f.vsol.Xt(1, 1), 0.3);

    Scirocco foil = TrippedFoil(4.0, 0.2, 0.6);
    ASSERT_TRUE(foil.Solve());
    const auto& s = foil.GetState();
    EXPECT_TRUE(s.vsol.tripped[0]);
    EXPECT_TRUE(s.vsol.tripped[1]);
    EXPECT_NEAR(s.vsol.Xt(1, 1), 0.2, 1e-12);
    EXPECT_NEAR(s.vsol.Xt(0, 1), 0.6, 1e-12);
    EXPECT_GT(foil.GetCd(), 1.2 * free_foil.GetCd());

    /* Laminar up to the trip interval, turbulent from it on */
    for (int si = 0; si < 2; ++si) {
        const Trip trip = FindTrip(s, si);
        const auto& Is = s.vsol.Is[si];
        ASSERT_GT(trip.station, 0);
        for (int i = 0; i < static_cast<int>(Is.size()); ++i)
            EXPECT_EQ(s.vsol.turb[Is[i]], i >= trip.station)
                << "surface " << si << " station " << i;
    }
}

/* XFOIL's rule: whichever comes first. A trip behind free transition
 * changes nothing. */
TEST(ForcedTransition, TripBehindFreeTransitionIsInactive) {
    Scirocco free_foil = TrippedFoil(4.0, 1.0, 1.0);
    ASSERT_TRUE(free_foil.Solve());
    Scirocco foil = TrippedFoil(4.0, 0.7, 1.0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_FALSE(foil.GetState().vsol.tripped[1]);
    EXPECT_NEAR(foil.GetCl(), free_foil.GetCl(), 1e-9);
    EXPECT_NEAR(foil.GetCd(), free_foil.GetCd(), 1e-10);
    EXPECT_NEAR(foil.GetState().vsol.Xt(1, 1),
                free_foil.GetState().vsol.Xt(1, 1), 1e-8);
}

/* The solution is continuous where the trip takes over from free
 * transition: a trip just ahead of the free location barely changes it. */
TEST(ForcedTransition, ContinuousAtFreeTransitionLocation) {
    Scirocco free_foil = TrippedFoil(4.0, 1.0, 1.0);
    ASSERT_TRUE(free_foil.Solve());
    const double xt = free_foil.GetState().vsol.Xt(1, 1);

    Scirocco ahead = TrippedFoil(4.0, xt - 1e-6, 1.0);
    ASSERT_TRUE(ahead.Solve());
    EXPECT_TRUE(ahead.GetState().vsol.tripped[1]);
    EXPECT_NEAR(ahead.GetCd(), free_foil.GetCd(), 1e-7);
    EXPECT_NEAR(ahead.GetCl(), free_foil.GetCl(), 1e-5);

    Scirocco behind = TrippedFoil(4.0, xt + 1e-6, 1.0);
    ASSERT_TRUE(behind.Solve());
    EXPECT_FALSE(behind.GetState().vsol.tripped[1]);
    EXPECT_NEAR(behind.GetCd(), free_foil.GetCd(), 1e-10);
}

TEST(ForcedTransition, FindTripLocatesInterval) {
    Scirocco foil = TrippedFoil(4.0, 0.2, 0.6);
    ASSERT_TRUE(foil.Solve());
    FoilState s = foil.GetState();
    for (int si = 0; si < 2; ++si) {
        const double xtrip = si == 0 ? 0.6 : 0.2;
        const Trip trip = FindTrip(s, si);
        const auto& Is = s.vsol.Is[si];
        ASSERT_GT(trip.station, 0);
        EXPECT_FALSE(trip.at_node);
        const double x1 = s.foil.x(0, Is[trip.station - 1]);
        const double x2 = s.foil.x(0, Is[trip.station]);
        EXPECT_LT(x1, xtrip);
        EXPECT_LE(xtrip, x2);
        EXPECT_NEAR(x1 + trip.frac * (x2 - x1), xtrip, 1e-14);
    }

    /* No trip: at 1 or more, or behind the trailing-edge node */
    s.param.xtrip_upper = 1.0;
    EXPECT_LT(FindTrip(s, 1).station, 0);
    s.param.xtrip_upper = 5.0;
    EXPECT_LT(FindTrip(s, 1).station, 0);
    EXPECT_LT(FindTrip(s, 2).station, 0); /* The wake has none */

    /* At alpha = 4 the stagnation point is on the lower side behind
     * x/c = 0.001: past that trip for the lower surface (XFOIL leaves the
     * side free), while the upper surface, which runs around the leading
     * edge, finds it on its own side. */
    s.param.xtrip_lower = 0.001;
    s.param.xtrip_upper = 0.001;
    EXPECT_LT(FindTrip(s, 0).station, 0);
    const Trip upper = FindTrip(s, 1);
    ASSERT_GT(upper.station, 0);
    EXPECT_GT(s.foil.x(1, s.vsol.Is[1][upper.station]), 0.0);

    /* Ahead of the leading edge: at the leading-edge node */
    s.param.xtrip_upper = -0.1;
    const Trip nose = FindTrip(s, 1);
    ASSERT_GT(nose.station, 0);
    EXPECT_TRUE(nose.at_node);
    EXPECT_EQ(nose.frac, 1.0);
}

/* A new trip on a solved object (warm start) gives the cold-start result. */
TEST(ForcedTransition, WarmStartFollowsNewTrip) {
    Scirocco warm = TrippedFoil(4.0, 1.0, 1.0);
    ASSERT_TRUE(warm.Solve());
    warm.GetState().param.xtrip_upper = 0.2;
    warm.GetState().param.xtrip_lower = 0.6;
    ASSERT_TRUE(warm.Solve());

    Scirocco cold = TrippedFoil(4.0, 0.2, 0.6);
    ASSERT_TRUE(cold.Solve());
    EXPECT_NEAR(warm.GetCl(), cold.GetCl(), 1e-8);
    EXPECT_NEAR(warm.GetCd(), cold.GetCd(), 1e-9);

    /* And back to free transition */
    warm.GetState().param.xtrip_upper = 1.0;
    warm.GetState().param.xtrip_lower = 1.0;
    ASSERT_TRUE(warm.Solve());
    Scirocco free_foil = TrippedFoil(4.0, 1.0, 1.0);
    ASSERT_TRUE(free_foil.Solve());
    EXPECT_NEAR(warm.GetCd(), free_foil.GetCd(), 1e-9);
}

/* Without the boundary-layer march (mfoil's iteration) the trip gives the
 * same solution. */
TEST(ForcedTransition, SameSolutionWithoutMarch) {
    Scirocco march = TrippedFoil(2.0, 0.3, 0.3);
    ASSERT_TRUE(march.Solve());
    Scirocco plain = TrippedFoil(2.0, 0.3, 0.3);
    plain.GetState().param.bl_march = false;
    ASSERT_TRUE(plain.Solve());
    EXPECT_NEAR(march.GetCl(), plain.GetCl(), 1e-8);
    EXPECT_NEAR(march.GetCd(), plain.GetCd(), 1e-9);
}

TEST(ForcedTransition, CapturedModelRejectsTrip) {
    Scirocco foil = TrippedFoil(4.0, 0.2, 1.0);
    foil.GetState().param.transition = TransitionModel::Captured;
    EXPECT_THROW(foil.Solve(), std::invalid_argument);
}
