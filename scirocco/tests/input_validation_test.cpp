#include <gtest/gtest.h>

#include <limits>

#include "scirocco/scirocco.h"

using namespace scirocco;

TEST(InputValidation, InvalidOperatingConditionsInvalidateSolution) {
    Scirocco foil("0012", 60);
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity()}) {
        for (int parameter = 0; parameter < 4; ++parameter) {
            foil.SetOper(3.0);
            ASSERT_TRUE(foil.Solve());
            foil.ComputeSensitivities();
            double values[] = {3.0, -1.0, 0.0, 0.0};
            values[parameter] = invalid;
            EXPECT_THROW(
                foil.SetOper(values[0], values[1], values[2], values[3]),
                std::invalid_argument);
            EXPECT_FALSE(foil.IsConverged());
            EXPECT_THROW((void)foil.GetDclDx(), std::runtime_error);
            EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
        }
    }
    EXPECT_THROW(foil.SetOper(0, -1, 0, -0.1), std::invalid_argument);
    EXPECT_THROW(foil.SetOper(0, -1, 0, 1.0), std::invalid_argument);
    // Preserve finite negative sentinels for inviscid/alpha mode.
    foil.SetOper(3, -1, -1, 0);
    EXPECT_TRUE(foil.Solve());
}

TEST(InputValidation, InvalidInternalStateCannotReportConvergence) {
    Scirocco foil("0012", 60);
    foil.SetOper(3);
    ASSERT_TRUE(foil.Solve());
    foil.GetState().oper.alpha = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(foil.Solve(), std::invalid_argument);
    EXPECT_FALSE(foil.IsConverged());
    foil.SetOper(3);
    // Force nonfinite force coefficients through the advanced state API.
    foil.GetState().geom.chord = 0;
    EXPECT_THROW(foil.Solve(), std::runtime_error);
    EXPECT_FALSE(foil.IsConverged());
}

TEST(InputValidation, RejectsUnsafePanelCountsAndCoordinates) {
    for (int n : {-1, 0, 1, 3, 4097, std::numeric_limits<int>::max()})
        EXPECT_THROW(Scirocco("0012", n), std::invalid_argument);
    EXPECT_THROW(Scirocco("ab12", 60), std::invalid_argument);
    EXPECT_THROW(Scirocco(MatrixXd::Zero(2, 3), 60), std::invalid_argument);
    EXPECT_THROW(Scirocco(MatrixXd(2, 0), 60), std::invalid_argument);
    EXPECT_THROW(Scirocco(MatrixXd::Zero(3, 6), 60), std::invalid_argument);
    EXPECT_THROW(Scirocco(MatrixXd::Zero(2, 6), 60), std::invalid_argument);
    Scirocco valid("0012", 60);
    MatrixXd coords = valid.GetFoilCoords();
    coords.col(2) = coords.col(1);
    EXPECT_THROW(Scirocco(coords, 60), std::invalid_argument);
    coords = valid.GetFoilCoords();
    coords(0, 2) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(Scirocco(coords, 60), std::invalid_argument);
}

TEST(SolverReuse, ViscousInviscidViscousMatchesFreshSolvers) {
    Scirocco reused("0012", 100);
    reused.SetOper(3, -1, 1e6);
    ASSERT_TRUE(reused.Solve());
    ASSERT_GT(reused.GetCf().norm(), 0);
    reused.SetOper(4);
    ASSERT_TRUE(reused.Solve());
    EXPECT_EQ(reused.GetCf().size(), 0);
    EXPECT_EQ(reused.GetTheta().size(), 0);
    EXPECT_EQ(reused.GetDeltaStar().size(), 0);
    EXPECT_EQ(reused.GetWakeCoords().cols(), 0);
    EXPECT_EQ(reused.GetState().oper.Re, 0);
    Scirocco fresh("0012", 100);
    fresh.SetOper(4);
    ASSERT_TRUE(fresh.Solve());
    EXPECT_NEAR(reused.GetCl(), fresh.GetCl(), 1e-12);
    EXPECT_NEAR(reused.GetCm(), fresh.GetCm(), 1e-12);
    EXPECT_EQ(reused.GetCd(), fresh.GetCd());
    reused.SetOper(2, -1, 1e6);
    fresh.SetOper(2, -1, 1e6);
    ASSERT_TRUE(reused.Solve());
    ASSERT_TRUE(fresh.Solve());
    EXPECT_NEAR(reused.GetCl(), fresh.GetCl(), 1e-12);
    EXPECT_NEAR(reused.GetCd(), fresh.GetCd(), 1e-12);
    EXPECT_LT((reused.GetCf() - fresh.GetCf()).norm(), 1e-12);
}
