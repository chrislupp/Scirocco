#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include "scirocco/core/wake.h"
#include "scirocco/scirocco.h"

using namespace scirocco;

namespace {
/* NACA 0012, 161 cosine-spaced points ordered lower TE -> LE -> upper TE,
 * with the last thickness coefficient set so the total TE gap is `gap`
 * (chord 1). gap = 0 closes the trailing edge exactly. */
MatrixXd Naca0012WithGap(double gap) {
    const int n = 81;
    const double a4 = -0.1036 + gap / 1.2;
    auto yt = [a4](double x) {
        return 0.6 * (0.2969 * std::sqrt(x) - 0.126 * x - 0.3516 * x * x +
                      0.2843 * x * x * x + a4 * x * x * x * x);
    };
    MatrixXd coords(2, 2 * n - 1);
    for (int i = 0; i < n; ++i) {
        const double x = 0.5 * (1.0 - std::cos(M_PI * i / (n - 1)));
        coords.col(n - 1 - i) << x, -yt(x);
        coords.col(n - 1 + i) << x, yt(x);
    }
    return coords;
}
}  // namespace

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

TEST(InputValidation, CamberRejectsMalformedDataWithoutChangingGeometry) {
    Scirocco foil("0012", 60);
    foil.SetOper(3);
    ASSERT_TRUE(foil.Solve());
    const MatrixXd before = foil.GetFoilCoords();
    const MatrixXd source = foil.GetState().geom.xpoint;
    std::vector<MatrixXd> invalid = {MatrixXd(2, 0), MatrixXd(0, 2),
                                     MatrixXd::Zero(3, 3), MatrixXd::Zero(2, 1),
                                     MatrixXd::Zero(1, 2)};
    MatrixXd camber(2, 3);
    camber << 0, 0.5, 1, 0, 0.02, 0;
    for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::quiet_NaN()}) {
        MatrixXd data = camber;
        data(0, 1) = bad;
        invalid.push_back(data);
    }
    camber(1, 1) = std::numeric_limits<double>::quiet_NaN();
    invalid.push_back(camber);
    for (const auto& data : invalid) {
        EXPECT_THROW(foil.GeomAddCamber(data), std::invalid_argument);
        EXPECT_EQ(foil.GetFoilCoords(), before);
        EXPECT_EQ(foil.GetState().geom.xpoint, source);
    }
}

TEST(InputValidation, CamberAcceptsBothLayoutsAndTwoPointInput) {
    Scirocco rows("0012", 60), columns("0012", 60);
    MatrixXd camber(2, 3);
    camber << 0, 0.5, 1, 0, 0.02, 0;
    rows.GeomAddCamber(camber);
    columns.GeomAddCamber(camber.transpose());
    EXPECT_EQ(rows.GetFoilCoords(), columns.GetFoilCoords());
    rows.SetOper(3);
    EXPECT_TRUE(rows.Solve());
    MatrixXd two(2, 2);
    two << 0, 1, 0, 0.01;
    EXPECT_NO_THROW(columns.GeomAddCamber(two));
    columns.SetOper(3);
    EXPECT_TRUE(columns.Solve());
}

/* Issue #18: the viscous wake model degenerates on a closed trailing edge.
 * Without the check these solves fail or, on some coordinates, report
 * convergence at a wrong Cl (0.03 instead of 0.21 at 199 panels, alpha 2). */
TEST(InputValidation, ViscousRejectsClosedTrailingEdge) {
    for (double gap : {0.0, 1e-9, 0.5 * kMinViscousTEGap}) {
        for (int npanel : {120, 199}) {
            Scirocco foil(Naca0012WithGap(gap), npanel);
            foil.SetOper(2.0, -1.0, 1e6, 0.0);
            EXPECT_THROW(foil.Solve(), std::invalid_argument)
                << "gap " << gap << ", " << npanel << " panels";
            EXPECT_FALSE(foil.IsConverged());
            EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
        }
    }
}

TEST(InputValidation, ClosedTrailingEdgeStillSolvesInviscid) {
    Scirocco foil(Naca0012WithGap(0.0), 160);
    foil.SetOper(2.0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_NEAR(foil.GetCl(), 0.2414, 1e-3);
}

TEST(InputValidation, ViscousAcceptsSmallTrailingEdgeGap) {
    /* 1e-5 chord: the resolution of a 5-decimal coordinate file. */
    Scirocco foil(Naca0012WithGap(1e-5), 160);
    foil.SetOper(2.0, -1.0, 1e6, 0.0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_NEAR(foil.GetCl(), 0.1982, 1e-3);
}
