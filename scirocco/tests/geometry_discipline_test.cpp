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

    GeometryDiscipline tests: NACA-0012 reproduction, LE invariant, and
    central-FD vs analytic-partial agreement.
*/

#include "scirocco/geometry_discipline.h"

#include <google/protobuf/struct.pb.h>
#include <gtest/gtest.h>
#include <variable.h>

#include <Eigen/Core>
#include <cmath>
#include <utility>
#include <vector>

#include "scirocco/core/geometry.h"
#include "scirocco/data/geom.h"

using namespace Eigen;

namespace {

google::protobuf::Struct BuildOptions(int n_camber_cp, int n_thickness_cp,
                                      int npanel, bool closed_te) {
    google::protobuf::Struct opts;
    (*opts.mutable_fields())["n_camber_cp"].set_number_value(n_camber_cp);
    (*opts.mutable_fields())["n_thickness_cp"].set_number_value(
        n_thickness_cp);
    (*opts.mutable_fields())["npanel"].set_number_value(npanel);
    (*opts.mutable_fields())["closed_te"].set_bool_value(closed_te);
    return opts;
}

void RunLifecycle(scirocco::GeometryDiscipline& disc,
                  const google::protobuf::Struct& opts,
                  philote::Variables& inputs, philote::Variables& outputs,
                  philote::Partials& partials) {
    disc.Initialize();
    disc.SetOptions(opts);
    disc.Setup();
    disc.SetupPartials();

    for (const auto& meta : disc.var_meta()) {
        std::vector<size_t> shape;
        for (auto dim : meta.shape()) shape.push_back(static_cast<size_t>(dim));

        if (meta.type() == philote::kInput)
            inputs[meta.name()] = philote::Variable(philote::kInput, shape);
        else if (meta.type() == philote::kOutput)
            outputs[meta.name()] = philote::Variable(philote::kOutput, shape);
    }

    for (const auto& meta : disc.partials_meta()) {
        std::vector<size_t> shape;
        for (auto dim : meta.shape()) shape.push_back(static_cast<size_t>(dim));

        partials[std::make_pair(meta.name(), meta.subname())] =
            philote::Variable(philote::kPartial, shape);
    }
}

// NACA finite-TE thickness law (same as scirocco::NacaPoints): half-thickness
// scaled so the digits "tt" (last two NACA digits) give tmax = 0.01·tt.
double NacaT(double x, double tmax) {
    const double t0 = 0.2969 * std::sqrt(x) - 0.126 * x - 0.3516 * x * x +
                      0.2843 * x * x * x - 0.1015 * x * x * x * x;
    return t0 * (tmax / 0.2);
}

}  // namespace

// ---------------------------------------------------------------------------
// Basic invariants
// ---------------------------------------------------------------------------

TEST(GeometryDisciplineTest, SetupShapesAreConsistent) {
    scirocco::GeometryDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 8, 200, true), inputs, outputs,
                 partials);

    EXPECT_EQ(disc.npanel(), 200);
    EXPECT_EQ(disc.npoint(), 201);
    EXPECT_EQ(inputs.at("camber_cp").Size(), 6u);
    EXPECT_EQ(inputs.at("thickness_cp").Size(), 8u);
    EXPECT_EQ(outputs.at("x").Size(), 2u * 201u);
}

TEST(GeometryDisciplineTest, OddNpanelRejected) {
    scirocco::GeometryDiscipline disc;
    disc.Initialize();
    EXPECT_THROW(disc.SetOptions(BuildOptions(6, 8, 199, true)),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// LE / TE invariants
// ---------------------------------------------------------------------------

TEST(GeometryDisciplineTest, LeadingEdgeAtOriginForZeroCamber) {
    scirocco::GeometryDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 8, 200, true), inputs, outputs,
                 partials);

    auto& cam = inputs.at("camber_cp");
    for (size_t i = 0; i < cam.Size(); ++i) cam(i) = 0.0;
    auto& thk = inputs.at("thickness_cp");
    for (size_t i = 0; i < thk.Size(); ++i) thk(i) = 0.05;

    disc.Compute(inputs, outputs);

    // LE index = N_per_side = npanel/2 = 100 (in both x-block and y-block)
    const int N = disc.npanel() / 2;
    const int npoint = disc.npoint();
    EXPECT_NEAR(outputs.at("x")(N), 0.0, 1e-14);                 // LE x
    EXPECT_NEAR(outputs.at("x")(npoint + N), 0.0, 1e-14);         // LE y
    // TE points (closed TE → both lower and upper TE at (1, 0))
    EXPECT_NEAR(outputs.at("x")(0), 1.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(npoint + 0), 0.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(npoint - 1), 1.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(2 * npoint - 1), 0.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Symmetric thickness reproduces NACA-0012 (open TE so spline can match
// the small NACA finite-TE gap exactly at xi = 1)
// ---------------------------------------------------------------------------

TEST(GeometryDisciplineTest, SymmetricMatchesNaca0012) {
    const int n_thk = 20;
    const int npanel = 200;
    scirocco::GeometryDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, n_thk, npanel, /*closed_te=*/false),
                 inputs, outputs, partials);

    // Zero camber
    auto& cam = inputs.at("camber_cp");
    for (size_t i = 0; i < cam.Size(); ++i) cam(i) = 0.0;

    // Thickness CPs = NACA-0012 half-thickness sampled at the free knot
    // positions. Free knot k corresponds to full-knot index 1+k.
    const auto& xi_knots = disc.xi_thickness_knots();
    auto& thk = inputs.at("thickness_cp");
    for (int k = 0; k < n_thk; ++k) {
        const double xi = xi_knots(1 + k);
        const double x = xi * xi;
        thk(k) = NacaT(x, 0.12);
    }

    disc.Compute(inputs, outputs);

    // Reference: scirocco NACA-0012 with same N=100 per side
    scirocco::Geom ref;
    scirocco::NacaPoints(ref, "0012");
    ASSERT_EQ(ref.npoint, disc.npoint());

    const int npoint = disc.npoint();
    double max_err = 0.0;
    for (int i = 0; i < npoint; ++i) {
        const double dx = outputs.at("x")(i) - ref.xpoint(0, i);
        const double dy = outputs.at("x")(npoint + i) - ref.xpoint(1, i);
        max_err = std::max(max_err, std::sqrt(dx * dx + dy * dy));
    }
    // Spline interpolates samples but not points in between; tolerance
    // reflects spline-vs-polynomial mismatch.
    EXPECT_LT(max_err, 5e-4);
}

// ---------------------------------------------------------------------------
// Analytic partials match central finite differences
// ---------------------------------------------------------------------------

TEST(GeometryDisciplineTest, PartialsMatchCentralFD) {
    scirocco::GeometryDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 8, 40, true), inputs, outputs, partials);

    const int n_cam = disc.n_camber_cp();
    const int n_thk = disc.n_thickness_cp();
    const int npoint = disc.npoint();

    // Non-trivial baseline: smooth low-amplitude camber, NACA-ish thickness
    const auto& xcc = disc.x_camber_knots();
    auto& cam = inputs.at("camber_cp");
    for (int j = 0; j < n_cam; ++j) {
        cam(j) = 0.04 * std::sin(M_PI * xcc(j));
    }
    const auto& xi_knots = disc.xi_thickness_knots();
    auto& thk = inputs.at("thickness_cp");
    for (int k = 0; k < n_thk; ++k) {
        const double xi = xi_knots(1 + k);
        thk(k) = NacaT(xi * xi, 0.12);
    }

    disc.ComputePartials(inputs, partials);
    const auto& d_c = partials[{"x", "camber_cp"}];
    const auto& d_t = partials[{"x", "thickness_cp"}];

    // Central-FD helper: perturb one CP, run Compute, return output vector
    auto eval_with = [&](const std::string& input_name, int idx,
                         double delta) {
        const double save = inputs.at(input_name)(idx);
        inputs.at(input_name)(idx) = save + delta;
        philote::Variables out_local = outputs;
        disc.Compute(inputs, out_local);
        VectorXd v(2 * npoint);
        for (int i = 0; i < 2 * npoint; ++i) v(i) = out_local.at("x")(i);
        inputs.at(input_name)(idx) = save;
        return v;
    };

    const double h = 1e-7;
    double max_err_c = 0.0;
    for (int j = 0; j < n_cam; ++j) {
        VectorXd up = eval_with("camber_cp", j, +h);
        VectorXd dn = eval_with("camber_cp", j, -h);
        VectorXd fd = (up - dn) / (2.0 * h);
        for (int row = 0; row < 2 * npoint; ++row) {
            const double an = d_c(row * n_cam + j);
            max_err_c = std::max(max_err_c, std::abs(an - fd(row)));
        }
    }
    double max_err_t = 0.0;
    for (int k = 0; k < n_thk; ++k) {
        VectorXd up = eval_with("thickness_cp", k, +h);
        VectorXd dn = eval_with("thickness_cp", k, -h);
        VectorXd fd = (up - dn) / (2.0 * h);
        for (int row = 0; row < 2 * npoint; ++row) {
            const double an = d_t(row * n_thk + k);
            max_err_t = std::max(max_err_t, std::abs(an - fd(row)));
        }
    }

    EXPECT_LT(max_err_c, 1e-6);
    EXPECT_LT(max_err_t, 1e-6);
}
