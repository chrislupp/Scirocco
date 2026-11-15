/*
    Scirocco: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include <H5Cpp.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "scirocco/core/adjoint.h"
#include "scirocco/core/paneling.h"
#include "scirocco/core/spline.h"
#include "scirocco/core/viscous_geom_deriv.h"
#include "scirocco/core/viscous_residual_deriv.h"
#include "scirocco/scirocco.h"

using namespace scirocco;
namespace {
MatrixXd ReadMatrix(H5::Group& group, const char* name) {
    auto dataset = group.openDataSet(name);
    hsize_t dims[2];
    dataset.getSpace().getSimpleExtentDims(dims);
    std::vector<double> values(dims[0] * dims[1]);
    dataset.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    MatrixXd out(dims[0], dims[1]);
    for (size_t i = 0; i < dims[0]; ++i)
        for (size_t j = 0; j < dims[1]; ++j)
            out(i, j) = values[i * dims[1] + j];
    return out;
}
double Attribute(H5::Group& group, const char* name) {
    double out;
    group.openAttribute(name).read(H5::PredType::NATIVE_DOUBLE, &out);
    return out;
}
std::vector<MatrixXd> Directions(H5::Group& group) {
    auto data = group.openDataSet("directions");
    hsize_t dims[3];
    data.getSpace().getSimpleExtentDims(dims);
    std::vector<double> values(dims[0] * dims[1] * dims[2]);
    data.read(values.data(), H5::PredType::NATIVE_DOUBLE);
    std::vector<MatrixXd> out;
    for (size_t p = 0; p < dims[0]; ++p) {
        MatrixXd direction(2, dims[2]);
        for (size_t d = 0; d < 2; ++d)
            for (size_t i = 0; i < dims[2]; ++i)
                direction(d, i) = values[(p * 2 + d) * dims[2] + i];
        out.push_back(direction);
    }
    return out;
}
void Check(double got, double reference) {
    EXPECT_TRUE(std::isfinite(got));
    EXPECT_NEAR(got, reference, 2e-11 + 2e-8 * std::abs(reference));
}
class GradientContractTest : public ::testing::TestWithParam<const char*> {};

TEST_P(GradientContractTest, PythonSplineArcLengthChain) {
    H5::H5File file(GRADIENT_CONTRACT_REFERENCE_FILE, H5F_ACC_RDONLY);
    auto group = file.openGroup(GetParam());
    const MatrixXd x = ReadMatrix(group, "x");
    const MatrixXd jac = SplineArcLengthDeriv(x);
    const MatrixXd ref = ReadMatrix(group, "arc_length_derivatives");
    const auto directions = Directions(group);
    for (size_t p = 0; p < directions.size(); ++p) {
        SCOPED_TRACE(p);
        VectorXd direction =
            Eigen::Map<const VectorXd>(directions[p].data(), x.size());
        VectorXd got = jac * direction;
        EXPECT_LT((got - ref.row(p).transpose()).cwiseAbs().maxCoeff(), 3e-11);
    }
}

TEST_P(GradientContractTest, OutputsAndGradientsMatchPython) {
    H5::H5File file(GRADIENT_CONTRACT_REFERENCE_FILE, H5F_ACC_RDONLY);
    auto group = file.openGroup(GetParam());
    MatrixXd x = ReadMatrix(group, "x");
    Scirocco foil("0012", x.cols() - 1);
    auto& s = foil.GetState();
    s.foil.x = x;
    // This is the actual optimization input path, not the cumulative test
    // helper.
    RefreshPanelGeometry(s.foil);
    foil.SetOper(Attribute(group, "alpha"), -1, Attribute(group, "Re"),
                 Attribute(group, "Ma"));
    ASSERT_TRUE(foil.Solve());
    auto sens = foil.ComputeSensitivities();
    std::array<double, 3> ref;
    group.openDataSet("outputs").read(ref.data(), H5::PredType::NATIVE_DOUBLE);
    EXPECT_NEAR(foil.GetCl(), ref[0], 2e-11);
    EXPECT_NEAR(foil.GetCd(), ref[1], 2e-11);
    EXPECT_NEAR(foil.GetCm(), ref[2], 2e-11);
    const auto directions = Directions(group);
    MatrixXd shape = ReadMatrix(group, "shape_derivatives");
    for (size_t p = 0; p < directions.size(); ++p) {
        SCOPED_TRACE(p);
        Check((sens.dcl_dx.array() * directions[p].array()).sum(), shape(p, 0));
        Check((sens.dcd_dx.array() * directions[p].array()).sum(), shape(p, 1));
        Check((sens.dcm_dx.array() * directions[p].array()).sum(), shape(p, 2));
    }
    MatrixXd scalar = ReadMatrix(group, "scalar_derivatives");
    MatrixXd got(3, 3);
    got << sens.dcl_dalpha, sens.dcd_dalpha, sens.dcm_dalpha, sens.dcl_dMa,
        sens.dcd_dMa, sens.dcm_dMa, sens.dcl_dRe, sens.dcd_dRe, sens.dcm_dRe;
    for (int p = 0; p < 3; ++p)
        for (int j = 0; j < 3; ++j) {
            // Keep tiny Reynolds derivatives meaningful rather than hidden by
            // atol.
            double scale = p == 2 && s.oper.viscous ? s.oper.Re : 1;
            Check(got(p, j) * scale, scalar(p, j) * scale);
        }
}
INSTANTIATE_TEST_SUITE_P(Mfoil, GradientContractTest,
                         ::testing::Values("inviscid_te", "viscous_spline",
                                           "viscous_mach"));

TEST(SplineDerivativeTest, ScalingAndRigidMotionIdentities) {
    for (int n : {2, 3, 8}) {
        SCOPED_TRACE(n);
        MatrixXd points(2, n);
        for (int i = 0; i < n; ++i) {
            double t = static_cast<double>(i) / (n - 1);
            points.col(i) << t, .2 * std::sin(2 * t);
        }
        const MatrixXd jac = SplineArcLengthDeriv(points);
        const VectorXd s = Spline2DCreate(points).X.s;
        const VectorXd dilation =
            Eigen::Map<const VectorXd>(points.data(), 2 * n);
        EXPECT_LT((jac * dilation - s).norm(), 2e-12);
        MatrixXd direction = MatrixXd::Zero(2, n);
        for (int d = 0; d < 2; ++d) {
            direction.setZero();
            direction.row(d).setOnes();
            EXPECT_LT(
                (jac * Eigen::Map<const VectorXd>(direction.data(), 2 * n))
                    .norm(),
                2e-12);
        }
        direction.row(0) = -points.row(1);
        direction.row(1) = points.row(0);
        EXPECT_LT(
            (jac * Eigen::Map<const VectorXd>(direction.data(), 2 * n)).norm(),
            2e-12);
    }
}

TEST(GradientValidityTest, UnsolvedStateIsRejected) {
    Scirocco foil("0012", 60);
    EXPECT_FALSE(foil.GetState().glob.conv);
    EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
}
TEST(GradientValidityTest, ExceptionInvalidatesPreviousSolve) {
    Scirocco foil("0012", 60);
    foil.SetOper(3);
    ASSERT_TRUE(foil.Solve());
    foil.GetState().foil.N = 0;  // Forces SolveInviscid to reject the geometry.
    EXPECT_THROW(foil.Solve(), std::runtime_error);
    EXPECT_FALSE(foil.GetState().glob.conv);
    EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
}
TEST(GradientValidityTest, UnconvergedSolveIsRejected) {
    Scirocco foil("0012", 100);
    foil.GetState().param.niglob = 1;
    foil.SetOper(3, -1, 1e6, 0);
    ASSERT_FALSE(foil.Solve());
    EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
    EXPECT_THROW(ViscousSensitivities(foil.GetState()), std::runtime_error);
}
TEST(GradientValidityTest, OperatingChangeInvalidatesPreviousSolve) {
    Scirocco foil("0012", 60);
    foil.SetOper(3);
    ASSERT_TRUE(foil.Solve());
    foil.SetOper(4);
    EXPECT_THROW(foil.ComputeSensitivities(), std::runtime_error);
}
TEST(GradientValidityTest, CachedAccessorsRejectStaleGradients) {
    Scirocco foil("0012", 60);
    EXPECT_THROW((void)foil.GetDclDalpha(), std::runtime_error);
    foil.SetOper(3);
    ASSERT_TRUE(foil.Solve());
    const auto sens = foil.ComputeSensitivities();
    EXPECT_DOUBLE_EQ(foil.GetDclDalpha(), sens.dcl_dalpha);
    foil.SetOper(4);
    EXPECT_THROW((void)foil.GetDclDx(), std::runtime_error);
    ASSERT_TRUE(foil.Solve());
    EXPECT_THROW((void)foil.GetDclDalpha(), std::runtime_error);
    foil.ComputeSensitivities();
    foil.GetState().oper.alpha = 5;
    EXPECT_THROW((void)foil.GetDclDalpha(), std::runtime_error);
}
TEST(GradientValidityTest, InviscidDragIsIdenticallyZero) {
    Scirocco foil("2412", 80);
    foil.SetOper(5, -1, 0, .3);
    ASSERT_TRUE(foil.Solve());
    auto sens = foil.ComputeSensitivities();
    EXPECT_DOUBLE_EQ(foil.GetCd(), 0);
    EXPECT_DOUBLE_EQ(sens.dcd_dalpha, 0);
    EXPECT_DOUBLE_EQ(sens.dcd_dMa, 0);
    EXPECT_DOUBLE_EQ(sens.dcd_dx.norm(), 0);
}
TEST(GradientValidityTest, TrimLiftDerivativesRespectConstraint) {
    Scirocco foil("0012", 100);
    foil.SetOper(3, .5, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    EXPECT_NEAR(foil.GetCl(), .5, 1e-10);
    const auto sens = foil.ComputeSensitivities();
    EXPECT_DOUBLE_EQ(sens.dcl_dx.norm(), 0);
    EXPECT_DOUBLE_EQ(sens.dcl_dalpha, 0);
    EXPECT_DOUBLE_EQ(sens.dcl_dMa, 0);
    EXPECT_DOUBLE_EQ(sens.dcl_dRe, 0);
    EXPECT_DOUBLE_EQ(sens.dcl_dcltarget, 1);
}

/* A shared RueVJPContext must give exactly what a per-call build gives, for
 * several lambdas, including when the wake was traced at another alpha. */
void ExpectContextMatchesPerCallBuild(const FoilState& state) {
    const MatrixXd dfoil_s_dx = DFoilSDFoilX(state.foil);
    const RueVJPContext ctx = BuildRueVJPContext(state, dfoil_s_dx);
    const int Nsys = state.glob.Nsys;
    for (int k = 0; k < 3; ++k) {
        VectorXd lambda(Nsys);
        for (int i = 0; i < Nsys; ++i) lambda(i) = std::sin(1.0 + i * (k + 1));
        const VectorXd per_call = ComputeRueVJP(state, lambda, dfoil_s_dx);
        const VectorXd shared = ComputeRueVJP(state, lambda, ctx);
        ASSERT_EQ(per_call.size(), shared.size());
        for (int i = 0; i < per_call.size(); ++i)
            EXPECT_EQ(per_call(i), shared(i)) << "k=" << k << " i=" << i;
    }
}
TEST(RueVJPContextTest, MatchesPerCallBuild) {
    Scirocco foil("2412", 80);
    foil.SetOper(4, -1, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    ASSERT_FALSE(WakeTracedAtOtherAlpha(foil.GetState()));
    ExpectContextMatchesPerCallBuild(foil.GetState());
}
TEST(RueVJPContextTest, MatchesPerCallBuildWithFrozenWake) {
    Scirocco foil("2412", 80);
    foil.SetOper(3, .5, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    ASSERT_TRUE(WakeTracedAtOtherAlpha(foil.GetState()));
    ExpectContextMatchesPerCallBuild(foil.GetState());
}
/* The reverse wake sweep must reproduce the forward wake Jacobians, which
 * viscous_chain_test validates against mfoil complex step:
 *   foil_bar + dgam^T gam_bar == dwake_x^T xbar + dwake_t^T tbar. */
void ExpectWakeSweepVJPMatchesForward(const FoilState& state) {
    const int N = state.foil.N;
    const int Nw = state.wake.N;
    const int K = 3;
    const MatrixXd dfoil_s_dx = DFoilSDFoilX(state.foil);
    MatrixXd xbar(2 * Nw, K), tbar(2 * Nw, K);
    for (int r = 0; r < 2 * Nw; ++r)
        for (int k = 0; k < K; ++k) {
            xbar(r, k) = std::sin(1.3 * r + 0.7 * k + 0.1);
            tbar(r, k) = std::cos(0.9 * r - 1.1 * k);
        }
    const WakeSweepLinearization lin = LinearizeWakeSweep(state, dfoil_s_dx);
    MatrixXd foil_bar = MatrixXd::Zero(2 * N, K);
    MatrixXd gam_bar = MatrixXd::Zero(N + 1, K);
    WakeSweepVJP(lin, xbar, tbar, foil_bar, gam_bar);
    const MatrixXd reverse =
        foil_bar + DGammaDFoilXForWake(state).transpose() * gam_bar;
    const MatrixXd forward =
        DWakeXDFoilX(state, dfoil_s_dx).transpose() * xbar +
        DWakeTDFoilX(state, dfoil_s_dx).transpose() * tbar;
    const double scale = forward.cwiseAbs().maxCoeff();
    EXPECT_LT((reverse - forward).cwiseAbs().maxCoeff(), 1e-12 * scale);
}
TEST(WakeSweepVJPTest, MatchesForwardJacobians) {
    Scirocco foil("2412", 80);
    foil.SetOper(4, -1, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    ExpectWakeSweepVJPMatchesForward(foil.GetState());
}
TEST(WakeSweepVJPTest, MatchesForwardJacobiansWithFrozenWake) {
    Scirocco foil("2412", 80);
    foil.SetOper(3, .5, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    ASSERT_TRUE(
        LinearizeWakeSweep(foil.GetState(), DFoilSDFoilX(foil.GetState().foil))
            .traced_at_other_alpha);
    ExpectWakeSweepVJPMatchesForward(foil.GetState());
}
TEST(RueVJPContextTest, BatchColumnsMatchSingleVJPs) {
    Scirocco foil("2412", 80);
    foil.SetOper(4, -1, 1e6, 0);
    ASSERT_TRUE(foil.Solve());
    const FoilState& state = foil.GetState();
    const RueVJPContext ctx =
        BuildRueVJPContext(state, DFoilSDFoilX(state.foil));
    const int Nsys = state.glob.Nsys;
    MatrixXd lambdas(Nsys, 3);
    for (int i = 0; i < Nsys; ++i)
        for (int k = 0; k < 3; ++k) lambdas(i, k) = std::sin(0.3 + i * (k + 2));
    const MatrixXd batch = ComputeRueVJPBatch(state, lambdas, ctx);
    for (int k = 0; k < 3; ++k) {
        const VectorXd single =
            ComputeRueVJP(state, VectorXd(lambdas.col(k)), ctx);
        const double scale = single.cwiseAbs().maxCoeff();
        EXPECT_LT((batch.col(k) - single).cwiseAbs().maxCoeff(), 1e-13 * scale);
    }
}
TEST(RueVJPContextTest, RejectsContextFromAnotherState) {
    Scirocco a("2412", 80), b("2412", 100);
    a.SetOper(4, -1, 1e6, 0);
    b.SetOper(4, -1, 1e6, 0);
    ASSERT_TRUE(a.Solve());
    ASSERT_TRUE(b.Solve());
    const RueVJPContext ctx =
        BuildRueVJPContext(a.GetState(), DFoilSDFoilX(a.GetState().foil));
    const VectorXd lambda = VectorXd::Ones(b.GetState().glob.Nsys);
    EXPECT_THROW(ComputeRueVJP(b.GetState(), lambda, ctx),
                 std::invalid_argument);
}
}  // namespace
