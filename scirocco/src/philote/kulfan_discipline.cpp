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

    Kulfan / CST Geometry Discipline Implementation
*/

#include "scirocco/kulfan_discipline.h"

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

namespace scirocco {

namespace {

/**
 * @brief Binomial coefficient K^n_i = n! / (i!·(n-i)!) for n ≤ ~60 (we never
 *        need more than ~20 in practice). Returned as double since it's
 *        combined with fp products below.
 * @param n Total count.
 * @param i Selection count.
 * @return Binomial coefficient as double.
 */
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

}  // namespace

void KulfanDiscipline::Initialize() {
    AddOption("n_cst_upper", "int");
    AddOption("n_cst_lower", "int");
    AddOption("npanel", "int");
    AddOption("te_upper", "double");
    AddOption("te_lower", "double");
}

void KulfanDiscipline::SetOptions(const google::protobuf::Struct& options) {
    const auto& fields = options.fields();

    if (fields.count("n_cst_upper")) {
        n_cst_upper_ =
            static_cast<int>(fields.at("n_cst_upper").number_value());
    }
    if (fields.count("n_cst_lower")) {
        n_cst_lower_ =
            static_cast<int>(fields.at("n_cst_lower").number_value());
    }
    if (fields.count("npanel")) {
        npanel_ = static_cast<int>(fields.at("npanel").number_value());
    }
    if (fields.count("te_upper")) {
        te_upper_ = fields.at("te_upper").number_value();
    }
    if (fields.count("te_lower")) {
        te_lower_ = fields.at("te_lower").number_value();
    }

    if (npanel_ < 4 || npanel_ % 2 != 0) {
        throw std::invalid_argument(
            "KulfanDiscipline: npanel must be an even integer >= 4");
    }
    if (n_cst_upper_ < 3 || n_cst_lower_ < 3) {
        /* smoothness output is the second-difference vector of length n-2,
           which requires at least 3 CPs per surface. */
        throw std::invalid_argument(
            "KulfanDiscipline: n_cst_upper / n_cst_lower must be >= 3");
    }

    N_per_side_ = npanel_ / 2;
    npoint_ = npanel_ + 1;

    BuildEvalStations();
    BuildBasisMatrices();
}

void KulfanDiscipline::Configure() {}

void KulfanDiscipline::Setup() {
    AddInput("cst_upper", {n_cst_upper_}, "");
    AddInput("cst_lower", {n_cst_lower_}, "");
    AddOutput("x", {2 * npoint_}, "m");
    /* Second-difference smoothness:  D_i = A_{i+2} − 2·A_{i+1} + A_i
       for i = 0..n_cst-3. S''(ψ) is the Bernstein expansion of D scaled by
       n(n−1), so bounding |D| directly caps the leading curvature term and
       is linear in cst (constant tri-banded partial). */
    AddOutput("smoothness_upper", {n_cst_upper_ - 2}, "");
    AddOutput("smoothness_lower", {n_cst_lower_ - 2}, "");
}

void KulfanDiscipline::SetupPartials() {
    DeclarePartials("x", "cst_upper");
    DeclarePartials("x", "cst_lower");
    DeclarePartials("smoothness_upper", "cst_upper");
    DeclarePartials("smoothness_lower", "cst_lower");
}

void KulfanDiscipline::BuildEvalStations() {
    /* Same TE-bunching distribution used by NacaPoints — gives smooth LE/TE
       clustering so the CST output plugs straight into Sirocco's paneling. */
    const int N = N_per_side_;
    const double te = 1.5;
    psi_eval_.resize(N + 1);

    VectorXd f = VectorXd::LinSpaced(N + 1, 0.0, 1.0);
    for (int i = 0; i <= N; ++i) {
        const double fi = f(i);
        double psi = 1.0 - (te + 1.0) * fi * std::pow(1.0 - fi, te) -
                     std::pow(1.0 - fi, te + 1.0);
        if (psi < 0.0) psi = 0.0;  /* guard fp roundoff */
        if (psi > 1.0) psi = 1.0;
        psi_eval_(i) = psi;
    }
}

void KulfanDiscipline::BuildBasisMatrices() {
    const int N = N_per_side_;

    /* Φ[q, i] = C(ψ_q) · K^(n-1)_i · ψ_q^i · (1−ψ_q)^(n-1-i)
       where n is the number of CST coefficients for that surface
       (Bernstein order = n - 1). */
    auto build = [&](int n_cst, MatrixXd& Phi) {
        const int order = n_cst - 1;
        Phi.resize(N + 1, n_cst);
        VectorXd binom(n_cst);
        for (int i = 0; i < n_cst; ++i) binom(i) = Binomial(order, i);
        for (int q = 0; q <= N; ++q) {
            const double psi = psi_eval_(q);
            const double one_minus_psi = 1.0 - psi;
            const double C = std::pow(psi, kN1) * std::pow(one_minus_psi, kN2);
            for (int i = 0; i < n_cst; ++i) {
                Phi(q, i) = C * binom(i) * std::pow(psi, i) *
                            std::pow(one_minus_psi, order - i);
            }
        }
    };

    build(n_cst_upper_, Phi_u_);
    build(n_cst_lower_, Phi_l_);
}

void KulfanDiscipline::Compute(const philote::Variables& inputs,
                               philote::Variables& outputs) {
    const auto& upper_var = inputs.at("cst_upper");
    const auto& lower_var = inputs.at("cst_lower");

    VectorXd cst_u(n_cst_upper_);
    for (int i = 0; i < n_cst_upper_; ++i) cst_u(i) = upper_var(i);
    VectorXd cst_l(n_cst_lower_);
    for (int i = 0; i < n_cst_lower_; ++i) cst_l(i) = lower_var(i);

    const int N = N_per_side_;
    const VectorXd y_u_shape = Phi_u_ * cst_u;  /* C·S_u  (no TE term) */
    const VectorXd y_l_shape = Phi_l_ * cst_l;  /* C·S_l */

    auto& x_var = outputs.at("x");

    /* Sirocco order: [0..N] = lower surface from TE to LE (reversed),
       index N = LE (shared), [N+1..2N] = upper surface from LE+1 to TE. */
    for (int q = 0; q <= N; ++q) {
        const double psi = psi_eval_(q);
        const double y_u = y_u_shape(q) + psi * te_upper_;
        const double y_l = -y_l_shape(q) + psi * te_lower_;

        /* Lower (reversed): output index = N - q */
        const int out_l = N - q;
        x_var(out_l) = psi;
        x_var(npoint_ + out_l) = y_l;

        /* Upper: output index = N + q. Skip q == 0 since the LE point is the
           shared lower endpoint at index N (already written above). */
        if (q > 0) {
            const int out_u = N + q;
            x_var(out_u) = psi;
            x_var(npoint_ + out_u) = y_u;
        }
    }

    /* ─── Smoothness outputs (second differences of cst) ────────────────── */
    auto& sm_u = outputs.at("smoothness_upper");
    for (int i = 0; i < n_cst_upper_ - 2; ++i) {
        sm_u(i) = cst_u(i + 2) - 2.0 * cst_u(i + 1) + cst_u(i);
    }
    auto& sm_l = outputs.at("smoothness_lower");
    for (int i = 0; i < n_cst_lower_ - 2; ++i) {
        sm_l(i) = cst_l(i + 2) - 2.0 * cst_l(i + 1) + cst_l(i);
    }
}

void KulfanDiscipline::ComputePartials(const philote::Variables& /*inputs*/,
                                       philote::Partials& partials) {
    const int N = N_per_side_;

    auto& d_u = partials[{"x", "cst_upper"}];
    auto& d_l = partials[{"x", "cst_lower"}];

    /* Flat row-major layout: partial entry (row, col) lives at
       row * n_cst + col, with row indexing the flat output vector
       (x-block then y-block). */
    const int n_rows = 2 * npoint_;
    for (int row = 0; row < n_rows; ++row) {
        for (int i = 0; i < n_cst_upper_; ++i)
            d_u(row * n_cst_upper_ + i) = 0.0;
        for (int i = 0; i < n_cst_lower_; ++i)
            d_l(row * n_cst_lower_ + i) = 0.0;
    }

    /* The surface coordinates depend on the CST coefficients only through
       the y-block (the x-coordinate of every node is the fixed ψ_eval).
       y_lower(ψ_q) = -Φ_l · cst_lower + ψ·te_lower
         → ∂y_l/∂cst_lower_i = -Φ_l(q, i),  ∂y_l/∂cst_upper = 0
       y_upper(ψ_q) =  Φ_u · cst_upper + ψ·te_upper
         → ∂y_u/∂cst_upper_i =  Φ_u(q, i),  ∂y_u/∂cst_lower = 0 */
    for (int q = 0; q <= N; ++q) {
        const int out_l = N - q;
        const int yrow_l = npoint_ + out_l;
        for (int i = 0; i < n_cst_lower_; ++i) {
            d_l(yrow_l * n_cst_lower_ + i) = -Phi_l_(q, i);
        }
        if (q > 0) {
            const int out_u = N + q;
            const int yrow_u = npoint_ + out_u;
            for (int i = 0; i < n_cst_upper_; ++i) {
                d_u(yrow_u * n_cst_upper_ + i) = Phi_u_(q, i);
            }
        }
    }

    /* ─── Smoothness partials (constant tri-banded [1, −2, 1]) ─────────────
       D_i = A_{i+2} − 2 A_{i+1} + A_i, so the Jacobian has row-major shape
       [n_cst − 2, n_cst] with row i having +1 at col i, −2 at col i+1,
       +1 at col i+2. */
    auto fill_smoothness = [](philote::Variable& dst, int n_d, int n_cst) {
        for (int i = 0; i < n_d * n_cst; ++i) dst(i) = 0.0;
        for (int i = 0; i < n_d; ++i) {
            dst(i * n_cst + i) = 1.0;
            dst(i * n_cst + i + 1) = -2.0;
            dst(i * n_cst + i + 2) = 1.0;
        }
    };
    fill_smoothness(partials[{"smoothness_upper", "cst_upper"}],
                    n_cst_upper_ - 2, n_cst_upper_);
    fill_smoothness(partials[{"smoothness_lower", "cst_lower"}],
                    n_cst_lower_ - 2, n_cst_lower_);
}

}  // namespace scirocco
