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

    Kulfan / Class-Shape Transformation (CST) airfoil geometry discipline
*/
#pragma once

#include <explicit.h>

#include "scirocco/scirocco.h"

namespace scirocco {

/**
 * @brief Philote-Cpp explicit discipline that builds airfoil coordinates from
 *        upper / lower Kulfan CST (Class-Shape Transformation) Bernstein
 *        coefficients.
 *
 * Each surface is
 *
 *   y_upper(ψ) =  C(ψ) · S_upper(ψ) + ψ · y_te_upper
 *   y_lower(ψ) = -C(ψ) · S_lower(ψ) + ψ · y_te_lower
 *
 * with ψ = x/c ∈ [0, 1], class function `C(ψ) = ψ^0.5·(1−ψ)^1` (round LE,
 * sharp TE), and shape function `S(ψ) = Σ A_i · K^n_i · ψ^i · (1−ψ)^(n−i)`
 * (Bernstein polynomial of order n). The `A_i` are the design variables for
 * that surface — `n+1` of them per surface.
 *
 * Because Bernstein basis terms are non-negative on [0, 1], positive
 * coefficients give a non-negative interpolant *everywhere*, not just at
 * sample points — bounding `cst_upper[i] ≥ 0`, `cst_lower[i] ≥ 0` is
 * sufficient to guarantee a physical airfoil (with this sign convention)
 * for every intermediate optimizer iterate.
 *
 * @par Options
 * - `n_cst_upper` (int, default 6): number of upper-surface CST coefficients
 * - `n_cst_lower` (int, default 6): number of lower-surface CST coefficients
 * - `npanel` (int, default 200): number of panels (must be even);
 *   `npoint = npanel + 1`.
 * - `te_upper` (double, default 0.0): fixed upper-surface TE y-offset
 * - `te_lower` (double, default 0.0): fixed lower-surface TE y-offset
 *
 * @par Inputs
 * - `cst_upper` ({n_cst_upper}): upper-surface Bernstein coefficients
 * - `cst_lower` ({n_cst_lower}): lower-surface Bernstein coefficients
 *
 * @par Outputs
 * - `x` ({2·npoint}): panel node coordinates, x-block then y-block,
 *   TE_lower → LE → TE_upper ordering.
 */
class KulfanDiscipline : public philote::ExplicitDiscipline {
   public:
    KulfanDiscipline() = default;
    ~KulfanDiscipline() noexcept = default;

    void Initialize() override;
    void SetOptions(const google::protobuf::Struct& options) override;
    void Configure() override;
    void Setup() override;
    void SetupPartials() override;
    void Compute(const philote::Variables& inputs,
                 philote::Variables& outputs) override;
    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials& partials) override;

    // Test-friendly accessors
    [[nodiscard]] int npoint() const { return npoint_; }
    [[nodiscard]] int npanel() const { return npanel_; }
    [[nodiscard]] int n_cst_upper() const { return n_cst_upper_; }
    [[nodiscard]] int n_cst_lower() const { return n_cst_lower_; }
    [[nodiscard]] const VectorXd& psi_eval() const { return psi_eval_; }

   private:
    void BuildEvalStations();
    void BuildBasisMatrices();

    int n_cst_upper_ = 6;
    int n_cst_lower_ = 6;
    int npanel_ = 200;
    double te_upper_ = 0.0;
    double te_lower_ = 0.0;

    int N_per_side_ = 0;  // npanel_ / 2
    int npoint_ = 0;      // npanel_ + 1

    static constexpr double kN1 = 0.5;  // class function LE exponent
    static constexpr double kN2 = 1.0;  // class function TE exponent

    VectorXd psi_eval_;  // [N_per_side + 1] cosine-clustered ψ on [0, 1]

    // Φ_u[q, i] = C(ψ_q) · K^(n_u-1)_i · ψ_q^i · (1−ψ_q)^(n_u-1−i)
    // so that y_u(ψ_q) = (Φ_u · cst_upper)_q + ψ_q · te_upper
    MatrixXd Phi_u_;  // [N_per_side + 1, n_cst_upper_]
    MatrixXd Phi_l_;  // [N_per_side + 1, n_cst_lower_]
};

}  // namespace scirocco
