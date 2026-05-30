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

// Binomial coefficient K^n_i = n! / (i!·(n-i)!) for n ≤ ~60 (we never need
// more than ~20 in practice). Returned as double since it's combined with
// fp products below.
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
    AddOption("n_kappa_samples", "int");
    AddOption("psi_kappa_min", "double");
    AddOption("psi_kappa_max", "double");
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
    if (fields.count("n_kappa_samples")) {
        n_kappa_samples_ =
            static_cast<int>(fields.at("n_kappa_samples").number_value());
    }
    if (fields.count("psi_kappa_min")) {
        psi_kappa_min_ = fields.at("psi_kappa_min").number_value();
    }
    if (fields.count("psi_kappa_max")) {
        psi_kappa_max_ = fields.at("psi_kappa_max").number_value();
    }

    if (npanel_ < 4 || npanel_ % 2 != 0) {
        throw std::invalid_argument(
            "KulfanDiscipline: npanel must be an even integer >= 4");
    }
    if (n_cst_upper_ < 2 || n_cst_lower_ < 2) {
        throw std::invalid_argument(
            "KulfanDiscipline: n_cst_upper / n_cst_lower must be >= 2");
    }
    if (n_kappa_samples_ < 2) {
        throw std::invalid_argument(
            "KulfanDiscipline: n_kappa_samples must be >= 2");
    }
    if (psi_kappa_min_ <= 0.0 || psi_kappa_max_ <= psi_kappa_min_ ||
        psi_kappa_max_ > 1.0) {
        throw std::invalid_argument(
            "KulfanDiscipline: require 0 < psi_kappa_min < psi_kappa_max <= 1");
    }

    N_per_side_ = npanel_ / 2;
    npoint_ = npanel_ + 1;

    BuildEvalStations();
    BuildBasisMatrices();
    BuildKappaBasisMatrices();
}

void KulfanDiscipline::Configure() {}

void KulfanDiscipline::Setup() {
    AddInput("cst_upper", {n_cst_upper_}, "");
    AddInput("cst_lower", {n_cst_lower_}, "");
    AddOutput("x", {2 * npoint_}, "m");
    AddOutput("kappa_upper", {n_kappa_samples_}, "1/m");
    AddOutput("kappa_lower", {n_kappa_samples_}, "1/m");
}

void KulfanDiscipline::SetupPartials() {
    DeclarePartials("x", "cst_upper");
    DeclarePartials("x", "cst_lower");
    // kappa on a given surface is independent of the other surface's CPs.
    DeclarePartials("kappa_upper", "cst_upper");
    DeclarePartials("kappa_lower", "cst_lower");
}

void KulfanDiscipline::BuildEvalStations() {
    // Same TE-bunching distribution used by NacaPoints — gives smooth LE/TE
    // clustering so the CST output plugs straight into Sirocco's paneling.
    const int N = N_per_side_;
    const double te = 1.5;
    psi_eval_.resize(N + 1);

    VectorXd f = VectorXd::LinSpaced(N + 1, 0.0, 1.0);
    for (int i = 0; i <= N; ++i) {
        const double fi = f(i);
        double psi = 1.0 - (te + 1.0) * fi * std::pow(1.0 - fi, te) -
                     std::pow(1.0 - fi, te + 1.0);
        if (psi < 0.0) psi = 0.0;  // guard fp roundoff
        if (psi > 1.0) psi = 1.0;
        psi_eval_(i) = psi;
    }
}

void KulfanDiscipline::BuildBasisMatrices() {
    const int N = N_per_side_;

    // Φ[q, i] = C(ψ_q) · K^(n-1)_i · ψ_q^i · (1−ψ_q)^(n-1-i)
    // where n is the number of CST coefficients for that surface
    // (Bernstein order = n - 1).
    auto build = [&](int n_cst, MatrixXd& Phi) {
        const int order = n_cst - 1;
        Phi.resize(N + 1, n_cst);
        VectorXd binom(n_cst);
        for (int i = 0; i < n_cst; ++i) binom(i) = Binomial(order, i);
        for (int q = 0; q <= N; ++q) {
            const double psi = psi_eval_(q);
            const double one_minus_psi = 1.0 - psi;
            const double C =
                std::pow(psi, kN1) * std::pow(one_minus_psi, kN2);
            for (int i = 0; i < n_cst; ++i) {
                Phi(q, i) = C * binom(i) * std::pow(psi, i) *
                            std::pow(one_minus_psi, order - i);
            }
        }
    };

    build(n_cst_upper_, Phi_u_);
    build(n_cst_lower_, Phi_l_);
}

void KulfanDiscipline::BuildKappaBasisMatrices() {
    const int Nk = n_kappa_samples_;

    // Cosine-cluster ψ_kappa over [psi_kappa_min, psi_kappa_max] so we get
    // tight sampling near both LE-skip boundary and TE-skip boundary.
    psi_kappa_.resize(Nk);
    const double span = psi_kappa_max_ - psi_kappa_min_;
    for (int j = 0; j < Nk; ++j) {
        const double cos_arg = M_PI * j / (Nk - 1);
        psi_kappa_(j) = psi_kappa_min_ + 0.5 * span * (1.0 - std::cos(cos_arg));
    }

    // Class function and its first two derivatives at each ψ_kappa station.
    //   C(ψ)   = √ψ (1−ψ)
    //   C'(ψ)  = (1 − 3ψ) / (2√ψ)
    //   C''(ψ) = −(1 + 3ψ) / (4 ψ^(3/2))
    C_kappa_.resize(Nk);
    Cp_kappa_.resize(Nk);
    Cpp_kappa_.resize(Nk);
    for (int j = 0; j < Nk; ++j) {
        const double psi = psi_kappa_(j);
        const double sp = std::sqrt(psi);
        C_kappa_(j) = sp * (1.0 - psi);
        Cp_kappa_(j) = (1.0 - 3.0 * psi) / (2.0 * sp);
        Cpp_kappa_(j) = -(1.0 + 3.0 * psi) / (4.0 * psi * sp);
    }

    // Bernstein-basis values and their first two derivatives, by direct
    // differentiation of B_i(ψ) = K_i ψ^i (1−ψ)^(n−i):
    //   B'_i(ψ)  = K_i [ i ψ^(i−1) (1−ψ)^(n−i)
    //                    − (n−i) ψ^i (1−ψ)^(n−i−1) ]
    //   B''_i(ψ) = K_i [ i(i−1) ψ^(i−2) (1−ψ)^(n−i)
    //                    − 2 i (n−i) ψ^(i−1) (1−ψ)^(n−i−1)
    //                    + (n−i)(n−i−1) ψ^i (1−ψ)^(n−i−2) ]
    // ψ stays bounded away from 0 and 1 by psi_kappa_min / psi_kappa_max,
    // so all negative-exponent powers stay finite.
    auto build_kappa_basis = [&](int n_cst, MatrixXd& M0, MatrixXd& M1,
                                 MatrixXd& M2) {
        const int order = n_cst - 1;
        VectorXd binom(n_cst);
        for (int i = 0; i < n_cst; ++i) binom(i) = Binomial(order, i);

        M0.resize(Nk, n_cst);
        M1.resize(Nk, n_cst);
        M2.resize(Nk, n_cst);
        for (int j = 0; j < Nk; ++j) {
            const double psi = psi_kappa_(j);
            const double one_m = 1.0 - psi;
            for (int i = 0; i < n_cst; ++i) {
                const int p = i;          // power on ψ
                const int q = order - i;  // power on (1−ψ)

                const double pip = std::pow(psi, p);
                const double opq = std::pow(one_m, q);
                M0(j, i) = binom(i) * pip * opq;

                // B'
                double dval = 0.0;
                if (p >= 1) {
                    dval += p * std::pow(psi, p - 1) * opq;
                }
                if (q >= 1) {
                    dval -= q * pip * std::pow(one_m, q - 1);
                }
                M1(j, i) = binom(i) * dval;

                // B''
                double d2val = 0.0;
                if (p >= 2) {
                    d2val += p * (p - 1) * std::pow(psi, p - 2) * opq;
                }
                if (p >= 1 && q >= 1) {
                    d2val -= 2.0 * p * q * std::pow(psi, p - 1) *
                             std::pow(one_m, q - 1);
                }
                if (q >= 2) {
                    d2val += q * (q - 1) * pip * std::pow(one_m, q - 2);
                }
                M2(j, i) = binom(i) * d2val;
            }
        }
    };

    build_kappa_basis(n_cst_upper_, M0_u_, M1_u_, M2_u_);
    build_kappa_basis(n_cst_lower_, M0_l_, M1_l_, M2_l_);
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
    const VectorXd y_u_shape = Phi_u_ * cst_u;  // C·S_u  (no TE term)
    const VectorXd y_l_shape = Phi_l_ * cst_l;  // C·S_l

    auto& x_var = outputs.at("x");

    // Sirocco order: [0..N] = lower surface from TE to LE (reversed),
    // index N = LE (shared), [N+1..2N] = upper surface from LE+1 to TE.
    for (int q = 0; q <= N; ++q) {
        const double psi = psi_eval_(q);
        const double y_u = y_u_shape(q) + psi * te_upper_;
        const double y_l = -y_l_shape(q) + psi * te_lower_;

        // Lower (reversed): output index = N - q
        const int out_l = N - q;
        x_var(out_l) = psi;
        x_var(npoint_ + out_l) = y_l;

        // Upper: output index = N + q. Skip q == 0 since the LE point is the
        // shared lower endpoint at index N (already written above).
        if (q > 0) {
            const int out_u = N + q;
            x_var(out_u) = psi;
            x_var(npoint_ + out_u) = y_u;
        }
    }

    // ─── Curvature outputs ────────────────────────────────────────────────
    // For each ψ_kappa station compute
    //   y(ψ)  = ±C(ψ)·S(ψ) + ψ·y_te    (te term doesn't affect y'')
    //   y'(ψ) = ±(C'·S + C·S') + y_te
    //   y''(ψ) = ±(C''·S + 2 C'·S' + C·S'')
    //   κ(ψ)   = y'' / (1 + y'²)^(3/2)
    // Sign s = +1 for upper, −1 for lower.
    const VectorXd S_u = M0_u_ * cst_u;
    const VectorXd Sp_u = M1_u_ * cst_u;
    const VectorXd Spp_u = M2_u_ * cst_u;
    const VectorXd S_l = M0_l_ * cst_l;
    const VectorXd Sp_l = M1_l_ * cst_l;
    const VectorXd Spp_l = M2_l_ * cst_l;

    auto& kappa_u = outputs.at("kappa_upper");
    auto& kappa_l = outputs.at("kappa_lower");
    for (int j = 0; j < n_kappa_samples_; ++j) {
        const double C = C_kappa_(j);
        const double Cp = Cp_kappa_(j);
        const double Cpp = Cpp_kappa_(j);

        const double yp_u = (Cp * S_u(j) + C * Sp_u(j)) + te_upper_;
        const double ypp_u = Cpp * S_u(j) + 2.0 * Cp * Sp_u(j) + C * Spp_u(j);
        const double denom_u = std::pow(1.0 + yp_u * yp_u, 1.5);
        kappa_u(j) = ypp_u / denom_u;

        const double yp_l = -(Cp * S_l(j) + C * Sp_l(j)) + te_lower_;
        const double ypp_l = -(Cpp * S_l(j) + 2.0 * Cp * Sp_l(j) +
                               C * Spp_l(j));
        const double denom_l = std::pow(1.0 + yp_l * yp_l, 1.5);
        kappa_l(j) = ypp_l / denom_l;
    }
}

void KulfanDiscipline::ComputePartials(const philote::Variables& inputs,
                                       philote::Partials& partials) {
    const int N = N_per_side_;

    auto& d_u = partials[{"x", "cst_upper"}];
    auto& d_l = partials[{"x", "cst_lower"}];

    // Flat row-major layout: partial entry (row, col) lives at
    // row * n_cst + col, with row indexing the flat output vector
    // (x-block then y-block).
    const int n_rows = 2 * npoint_;
    for (int row = 0; row < n_rows; ++row) {
        for (int i = 0; i < n_cst_upper_; ++i)
            d_u(row * n_cst_upper_ + i) = 0.0;
        for (int i = 0; i < n_cst_lower_; ++i)
            d_l(row * n_cst_lower_ + i) = 0.0;
    }

    // The surface coordinates depend on the CST coefficients only through
    // the y-block (the x-coordinate of every node is the fixed ψ_eval).
    // y_lower(ψ_q) = -Φ_l · cst_lower + ψ·te_lower
    //   → ∂y_l/∂cst_lower_i = -Φ_l(q, i),  ∂y_l/∂cst_upper = 0
    // y_upper(ψ_q) =  Φ_u · cst_upper + ψ·te_upper
    //   → ∂y_u/∂cst_upper_i =  Φ_u(q, i),  ∂y_u/∂cst_lower = 0
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

    // ─── Curvature partials ──────────────────────────────────────────────
    // κ = y'' / (1 + y'²)^(3/2). With y = ±C·S + ψ·y_te, S linear in cst:
    //   ∂y'/∂A_i  = ±(C'·M0[j,i] + C·M1[j,i])
    //   ∂y''/∂A_i = ±(C''·M0[j,i] + 2 C'·M1[j,i] + C·M2[j,i])
    //   ∂κ/∂A_i   = ∂y''/∂A_i / (1+y'²)^(3/2)
    //              − 3·y'·y''·∂y'/∂A_i / (1+y'²)^(5/2)
    const auto& upper_var = inputs.at("cst_upper");
    const auto& lower_var = inputs.at("cst_lower");
    VectorXd cst_u(n_cst_upper_);
    for (int i = 0; i < n_cst_upper_; ++i) cst_u(i) = upper_var(i);
    VectorXd cst_l(n_cst_lower_);
    for (int i = 0; i < n_cst_lower_; ++i) cst_l(i) = lower_var(i);

    const VectorXd S_u = M0_u_ * cst_u;
    const VectorXd Sp_u = M1_u_ * cst_u;
    const VectorXd Spp_u = M2_u_ * cst_u;
    const VectorXd S_l = M0_l_ * cst_l;
    const VectorXd Sp_l = M1_l_ * cst_l;
    const VectorXd Spp_l = M2_l_ * cst_l;

    auto& d_ku = partials[{"kappa_upper", "cst_upper"}];
    auto& d_kl = partials[{"kappa_lower", "cst_lower"}];

    for (int j = 0; j < n_kappa_samples_; ++j) {
        const double C = C_kappa_(j);
        const double Cp = Cp_kappa_(j);
        const double Cpp = Cpp_kappa_(j);

        // Upper surface (sign +1)
        const double yp_u = (Cp * S_u(j) + C * Sp_u(j)) + te_upper_;
        const double ypp_u = Cpp * S_u(j) + 2.0 * Cp * Sp_u(j) + C * Spp_u(j);
        const double one_plus = 1.0 + yp_u * yp_u;
        const double D32_u = std::pow(one_plus, 1.5);
        const double D52_u = D32_u * one_plus;
        const double coeff_u = -3.0 * yp_u * ypp_u / D52_u;
        for (int i = 0; i < n_cst_upper_; ++i) {
            const double dyp = Cp * M0_u_(j, i) + C * M1_u_(j, i);
            const double dypp = Cpp * M0_u_(j, i) + 2.0 * Cp * M1_u_(j, i) +
                                C * M2_u_(j, i);
            d_ku(j * n_cst_upper_ + i) = dypp / D32_u + coeff_u * dyp;
        }

        // Lower surface (sign −1)
        const double yp_l = -(Cp * S_l(j) + C * Sp_l(j)) + te_lower_;
        const double ypp_l = -(Cpp * S_l(j) + 2.0 * Cp * Sp_l(j) +
                               C * Spp_l(j));
        const double one_plus_l = 1.0 + yp_l * yp_l;
        const double D32_l = std::pow(one_plus_l, 1.5);
        const double D52_l = D32_l * one_plus_l;
        const double coeff_l = -3.0 * yp_l * ypp_l / D52_l;
        for (int i = 0; i < n_cst_lower_; ++i) {
            const double dyp = -(Cp * M0_l_(j, i) + C * M1_l_(j, i));
            const double dypp = -(Cpp * M0_l_(j, i) + 2.0 * Cp * M1_l_(j, i) +
                                  C * M2_l_(j, i));
            d_kl(j * n_cst_lower_ + i) = dypp / D32_l + coeff_l * dyp;
        }
    }
}

}  // namespace scirocco
