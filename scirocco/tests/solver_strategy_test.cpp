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
