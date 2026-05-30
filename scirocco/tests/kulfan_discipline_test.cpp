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

    KulfanDiscipline tests: NACA-0012 CST-fit reproduction, LE/TE invariants,
    Bernstein-positivity guarantee, and central-FD vs analytic-partial
    agreement.
*/

#include "scirocco/kulfan_discipline.h"

#include <google/protobuf/struct.pb.h>
#include <gtest/gtest.h>
#include <variable.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <vector>

#include "scirocco/core/geometry.h"
#include "scirocco/data/geom.h"

using namespace Eigen;

namespace {

google::protobuf::Struct BuildOptions(int n_cst_upper, int n_cst_lower,
                                      int npanel, double te_upper,
                                      double te_lower) {
    google::protobuf::Struct opts;
    (*opts.mutable_fields())["n_cst_upper"].set_number_value(n_cst_upper);
    (*opts.mutable_fields())["n_cst_lower"].set_number_value(n_cst_lower);
    (*opts.mutable_fields())["npanel"].set_number_value(npanel);
    (*opts.mutable_fields())["te_upper"].set_number_value(te_upper);
    (*opts.mutable_fields())["te_lower"].set_number_value(te_lower);
    return opts;
}

void RunLifecycle(scirocco::KulfanDiscipline& disc,
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

double Binomial(int n, int i) {
    if (i < 0 || i > n) return 0.0;
    if (i > n - i) i = n - i;
    double result = 1.0;
    for (int k = 0; k < i; ++k) {
        result *= static_cast<double>(n - k);
        result /= static_cast<double>(k + 1);
    }
    return result;
}

double NacaT(double x, double tmax) {
    const double t0 = 0.2969 * std::sqrt(x) - 0.126 * x - 0.3516 * x * x +
                      0.2843 * x * x * x - 0.1015 * x * x * x * x;
    return t0 * (tmax / 0.2);
}

// Least-squares fit of CST upper-surface coefficients to a target
// half-thickness profile y_target(psi) with te_offset accounted for, so the
// fitted Bernstein polynomial S(psi) satisfies
//
//   C(psi) · S(psi) + psi · te_offset ≈ y_target(psi).
//
VectorXd FitCstSurface(int n_cst, double te_offset,
                       const std::vector<double>& psi_samples,
                       const std::vector<double>& y_samples) {
    const int order = n_cst - 1;
    const int M = static_cast<int>(psi_samples.size());

    MatrixXd A(M, n_cst);
    VectorXd b(M);
    VectorXd binom(n_cst);
    for (int i = 0; i < n_cst; ++i) binom(i) = Binomial(order, i);

    for (int j = 0; j < M; ++j) {
        const double psi = psi_samples[j];
        const double one_minus = 1.0 - psi;
        const double C = std::sqrt(psi) * one_minus;  // N1=0.5, N2=1.0
        for (int i = 0; i < n_cst; ++i) {
            A(j, i) = C * binom(i) * std::pow(psi, i) *
                      std::pow(one_minus, order - i);
        }
        b(j) = y_samples[j] - psi * te_offset;
    }
    // Normal equations (well-conditioned for small n_cst on a uniform sample)
    return (A.transpose() * A).ldlt().solve(A.transpose() * b);
}

}  // namespace

// ---------------------------------------------------------------------------
// Basic invariants
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, SetupShapesAreConsistent) {
    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 6, 200, 0.0, 0.0), inputs, outputs,
                 partials);

    EXPECT_EQ(disc.npanel(), 200);
    EXPECT_EQ(disc.npoint(), 201);
    EXPECT_EQ(inputs.at("cst_upper").Size(), 6u);
    EXPECT_EQ(inputs.at("cst_lower").Size(), 6u);
    EXPECT_EQ(outputs.at("x").Size(), 2u * 201u);
}

TEST(KulfanDisciplineTest, OddNpanelRejected) {
    scirocco::KulfanDiscipline disc;
    disc.Initialize();
    EXPECT_THROW(disc.SetOptions(BuildOptions(6, 6, 199, 0.0, 0.0)),
                 std::invalid_argument);
}

// ---------------------------------------------------------------------------
// LE / TE invariants
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, LeadingAndTrailingEdgeAtChordEndpoints) {
    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 6, 200, 0.0, 0.0), inputs, outputs,
                 partials);

    auto& up = inputs.at("cst_upper");
    auto& lo = inputs.at("cst_lower");
    for (size_t i = 0; i < up.Size(); ++i) up(i) = 0.15;
    for (size_t i = 0; i < lo.Size(); ++i) lo(i) = 0.15;

    disc.Compute(inputs, outputs);

    const int N = disc.npanel() / 2;
    const int npoint = disc.npoint();
    // LE at index N (shared)
    EXPECT_NEAR(outputs.at("x")(N), 0.0, 1e-14);
    EXPECT_NEAR(outputs.at("x")(npoint + N), 0.0, 1e-14);
    // TE at indices 0 (lower) and 2N (upper)
    EXPECT_NEAR(outputs.at("x")(0), 1.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(npoint + 0), 0.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(2 * N), 1.0, 1e-12);
    EXPECT_NEAR(outputs.at("x")(npoint + 2 * N), 0.0, 1e-12);
}

// ---------------------------------------------------------------------------
// NACA-0012 reproduction via least-squares CST fit
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, Naca0012ReproductionViaCstFit) {
    const int n_cst = 8;
    const int npanel = 200;

    // NACA-0012 has a small open TE; account for it via te_upper / te_lower.
    const double te_thk = NacaT(1.0, 0.12);  // half-thickness at psi = 1

    // Sample the NACA-0012 thickness on a dense grid for the fit
    const int n_fit = 200;
    std::vector<double> psi_samples(n_fit), y_upper_samples(n_fit),
        y_lower_samples(n_fit);
    for (int j = 0; j < n_fit; ++j) {
        const double psi = static_cast<double>(j) / (n_fit - 1);
        psi_samples[j] = psi;
        y_upper_samples[j] = NacaT(psi, 0.12);
        // y_lower(psi) = -C·S_l + psi·te_lower (with te_lower = -te_thk)
        // and we want y_lower = -NacaT (mirrored). Internally the discipline
        // computes -C·S_l, so we fit S_l against +NacaT(psi).
        y_lower_samples[j] = NacaT(psi, 0.12);
    }

    const VectorXd cst_u_fit =
        FitCstSurface(n_cst, te_thk, psi_samples, y_upper_samples);
    const VectorXd cst_l_fit =
        FitCstSurface(n_cst, te_thk, psi_samples, y_lower_samples);

    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(n_cst, n_cst, npanel, te_thk, -te_thk),
                 inputs, outputs, partials);

    auto& up = inputs.at("cst_upper");
    auto& lo = inputs.at("cst_lower");
    for (int i = 0; i < n_cst; ++i) up(i) = cst_u_fit(i);
    for (int i = 0; i < n_cst; ++i) lo(i) = cst_l_fit(i);

    disc.Compute(inputs, outputs);

    // Reference: NACA-0012 generated by the scirocco helper.
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
    // Order-7 (n_cst=8) CST captures NACA-0012 to well under 1e-3 with the
    // open-TE convention. Tighter tolerances need more coefficients.
    EXPECT_LT(max_err, 1e-3);
}

// ---------------------------------------------------------------------------
// Bernstein positivity: positive CPs → positive thickness everywhere, not
// just at panel nodes. Cross-check on a dense intermediate-ψ grid.
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, PositiveCpsGivePositiveThicknessEverywhere) {
    const int n_cst = 6;
    const int npanel = 40;

    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(n_cst, n_cst, npanel, 0.0, 0.0), inputs,
                 outputs, partials);

    // Wildly non-uniform but strictly positive coefficients (the kind a
    // cubic spline interpolant could overshoot through). Bernstein basis
    // should still produce a non-negative interpolant everywhere.
    std::vector<double> raw = {0.30, 0.05, 0.40, 0.02, 0.35, 0.04};
    auto& up = inputs.at("cst_upper");
    auto& lo = inputs.at("cst_lower");
    for (int i = 0; i < n_cst; ++i) {
        up(i) = raw[i];
        lo(i) = raw[i];
    }
    disc.Compute(inputs, outputs);

    // Panel nodes: lower y must be ≤ upper y (and at the LE / TE they're
    // equal). Walk over the lower-then-upper layout: for each chord station
    // q ∈ [0, N], lower is at index (N − q), upper is at index (N + q).
    const int N = disc.npanel() / 2;
    const int npoint = disc.npoint();
    for (int q = 1; q < N; ++q) {
        const double y_l = outputs.at("x")(npoint + (N - q));
        const double y_u = outputs.at("x")(npoint + (N + q));
        EXPECT_GT(y_u, y_l) << "Thickness non-positive at q=" << q;
    }

    // Dense intermediate-ψ cross-check by re-evaluating the closed-form CST
    // expression with the same coefficients, on 5000 uniform ψ samples.
    auto binom = [](int n, int i) {
        if (i < 0 || i > n) return 0.0;
        if (i > n - i) i = n - i;
        double r = 1.0;
        for (int k = 0; k < i; ++k) {
            r *= static_cast<double>(n - k);
            r /= static_cast<double>(k + 1);
        }
        return r;
    };
    const int order = n_cst - 1;
    std::vector<double> binoms(n_cst);
    for (int i = 0; i < n_cst; ++i) binoms[i] = binom(order, i);

    const int n_dense = 5000;
    double min_thk = 1e9;
    for (int j = 0; j < n_dense; ++j) {
        const double psi = static_cast<double>(j) / (n_dense - 1);
        const double one_minus = 1.0 - psi;
        const double C = std::sqrt(psi) * one_minus;
        double S = 0.0;
        for (int i = 0; i < n_cst; ++i) {
            S += raw[i] * binoms[i] * std::pow(psi, i) *
                 std::pow(one_minus, order - i);
        }
        const double thk = C * S;  // upper has +C·S, lower has -C·S
        min_thk = std::min(min_thk, thk);
    }
    EXPECT_GE(min_thk, 0.0);
}

// ---------------------------------------------------------------------------
// Analytical partials match central FD
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, PartialsMatchCentralFD) {
    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 6, 40, 0.001, -0.001), inputs, outputs,
                 partials);

    const int n_u = disc.n_cst_upper();
    const int n_l = disc.n_cst_lower();
    const int npoint = disc.npoint();

    // Non-trivial baseline
    std::vector<double> raw_u = {0.20, 0.18, 0.16, 0.20, 0.18, 0.16};
    std::vector<double> raw_l = {0.18, 0.16, 0.14, 0.18, 0.16, 0.14};
    auto& up = inputs.at("cst_upper");
    auto& lo = inputs.at("cst_lower");
    for (int i = 0; i < n_u; ++i) up(i) = raw_u[i];
    for (int i = 0; i < n_l; ++i) lo(i) = raw_l[i];

    disc.ComputePartials(inputs, partials);
    const auto& d_u = partials[{"x", "cst_upper"}];
    const auto& d_l = partials[{"x", "cst_lower"}];

    auto eval_with = [&](const std::string& name, int idx, double delta) {
        const double save = inputs.at(name)(idx);
        inputs.at(name)(idx) = save + delta;
        philote::Variables out_local = outputs;
        disc.Compute(inputs, out_local);
        VectorXd v(2 * npoint);
        for (int i = 0; i < 2 * npoint; ++i) v(i) = out_local.at("x")(i);
        inputs.at(name)(idx) = save;
        return v;
    };

    const double h = 1e-7;
    double max_err_u = 0.0;
    for (int i = 0; i < n_u; ++i) {
        VectorXd up_d = eval_with("cst_upper", i, +h);
        VectorXd dn_d = eval_with("cst_upper", i, -h);
        VectorXd fd = (up_d - dn_d) / (2.0 * h);
        for (int row = 0; row < 2 * npoint; ++row) {
            const double an = d_u(row * n_u + i);
            max_err_u = std::max(max_err_u, std::abs(an - fd(row)));
        }
    }
    double max_err_l = 0.0;
    for (int i = 0; i < n_l; ++i) {
        VectorXd up_d = eval_with("cst_lower", i, +h);
        VectorXd dn_d = eval_with("cst_lower", i, -h);
        VectorXd fd = (up_d - dn_d) / (2.0 * h);
        for (int row = 0; row < 2 * npoint; ++row) {
            const double an = d_l(row * n_l + i);
            max_err_l = std::max(max_err_l, std::abs(an - fd(row)));
        }
    }

    EXPECT_LT(max_err_u, 1e-8);
    EXPECT_LT(max_err_l, 1e-8);
}

// ---------------------------------------------------------------------------
// Smoothness output: D_i = A_{i+2} − 2 A_{i+1} + A_i, with constant
// linear partials. Verify both the values and the [1, −2, 1] Jacobian.
// ---------------------------------------------------------------------------

TEST(KulfanDisciplineTest, SmoothnessSecondDifferenceAndPartials) {
    scirocco::KulfanDiscipline disc;
    philote::Variables inputs, outputs;
    philote::Partials partials;
    RunLifecycle(disc, BuildOptions(6, 6, 40, 0.0, 0.0), inputs, outputs,
                 partials);

    const int n_u = disc.n_cst_upper();
    const int n_l = disc.n_cst_lower();
    std::vector<double> raw_u = {0.30, 0.05, 0.40, 0.02, 0.35, 0.04};
    std::vector<double> raw_l = {0.10, 0.20, 0.10, 0.20, 0.10, 0.20};
    for (int i = 0; i < n_u; ++i) inputs.at("cst_upper")(i) = raw_u[i];
    for (int i = 0; i < n_l; ++i) inputs.at("cst_lower")(i) = raw_l[i];

    disc.Compute(inputs, outputs);

    // Values match the closed-form second differences.
    for (int i = 0; i < n_u - 2; ++i) {
        const double expected =
            raw_u[i + 2] - 2.0 * raw_u[i + 1] + raw_u[i];
        EXPECT_NEAR(outputs.at("smoothness_upper")(i), expected, 1e-14);
    }
    for (int i = 0; i < n_l - 2; ++i) {
        const double expected =
            raw_l[i + 2] - 2.0 * raw_l[i + 1] + raw_l[i];
        EXPECT_NEAR(outputs.at("smoothness_lower")(i), expected, 1e-14);
    }

    // Constant Jacobian — verify the [1, −2, 1] band exactly.
    disc.ComputePartials(inputs, partials);
    const auto& d_su = partials[{"smoothness_upper", "cst_upper"}];
    for (int i = 0; i < n_u - 2; ++i) {
        for (int j = 0; j < n_u; ++j) {
            const double expected = (j == i)       ? 1.0
                                    : (j == i + 1) ? -2.0
                                    : (j == i + 2) ? 1.0
                                                   : 0.0;
            EXPECT_NEAR(d_su(i * n_u + j), expected, 1e-14)
                << "row=" << i << " col=" << j;
        }
    }
}
