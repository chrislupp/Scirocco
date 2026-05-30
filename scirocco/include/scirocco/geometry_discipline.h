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

    Spline-based airfoil geometry discipline (camber + perpendicular thickness)
*/
#pragma once

#include <explicit.h>

#include "scirocco/scirocco.h"

namespace scirocco {

/**
 * @brief Philote-Cpp explicit discipline that builds airfoil coordinates from
 *        a cubic-spline camber line and a thickness distribution offset
 *        perpendicular to the camber line.
 *
 * Camber is a NOT-A-KNOT cubic spline `y_c(x)` whose `n_camber_cp` control
 * points sit at fixed cosine-clustered x-stations on [0, 1].
 *
 * Thickness is treated as a parametric 2-D curve `(x(ξ), t(ξ))` with
 * `x(ξ) = ξ²`. Cubic-spline interpolation on ξ-knots `ξ_i = cos-cluster([0,1])`
 * reproduces `x(ξ) = ξ²` exactly (it is a degree-2 polynomial) and so
 * `dx/dξ(0) = 0` — the thickness curve terminates perpendicular to the chord
 * at the LE, which gives the airfoil a properly rounded leading edge.
 *
 * Half-thickness is offset along the camber-line normal at each chord station,
 * not added vertically — correct for moderate-to-high camber.
 *
 * The output uses Sirocco's standard `[2 × npoint]` ordering
 * (TE_lower → LE → TE_upper) flattened row-major into a single
 * `{2·npoint}` vector matching `SiroccoDiscipline`'s `x` input.
 *
 * @par Options
 * - `n_camber_cp` (int, default 6): number of camber control points / DVs
 * - `n_thickness_cp` (int, default 8): number of free thickness control
 *   points. The LE knot is always fixed at t=0, and the TE knot is fixed at
 *   t=0 when `closed_te` is true — neither is part of the input vector.
 * - `npanel` (int, default 200): number of panels; must be even.
 *   `npoint = npanel + 1`.
 * - `closed_te` (bool, default true): force `t = 0` at the TE knot.
 *
 * @par Inputs
 * - `camber_cp` ({n_camber_cp}): y-values of camber CPs
 * - `thickness_cp` ({n_thickness_cp}): t-values of free thickness CPs
 *
 * @par Outputs
 * - `x` ({2·npoint}): airfoil node coordinates, x-block then y-block
 */
class GeometryDiscipline : public philote::ExplicitDiscipline {
   public:
    GeometryDiscipline() = default;
    ~GeometryDiscipline() noexcept = default;

    void Initialize() override;
    void SetOptions(const google::protobuf::Struct& options) override;
    void Configure() override;
    void Setup() override;
    void SetupPartials() override;
    void Compute(const philote::Variables& inputs,
                 philote::Variables& outputs) override;
    void ComputePartials(const philote::Variables& inputs,
                         philote::Partials& partials) override;

    // Read-only accessors exposed for tests
    [[nodiscard]] int npoint() const { return npoint_; }
    [[nodiscard]] int npanel() const { return npanel_; }
    [[nodiscard]] int n_camber_cp() const { return n_camber_cp_; }
    [[nodiscard]] int n_thickness_cp() const { return n_thickness_cp_; }
    [[nodiscard]] const VectorXd& x_camber_knots() const {
        return x_camber_knots_;
    }
    [[nodiscard]] const VectorXd& xi_thickness_knots() const {
        return xi_thickness_knots_;
    }

   private:
    void BuildEvalStations();
    void BuildBasisMatrices();
    void ComputeSurface(const VectorXd& camber_cp,
                        const VectorXd& thickness_cp, MatrixXd& xpoint) const;

    int n_camber_cp_ = 6;
    int n_thickness_cp_ = 8;
    int npanel_ = 200;
    bool closed_te_ = true;

    int N_per_side_ = 0;
    int npoint_ = 0;

    VectorXd x_eval_;   // [N_per_side + 1] cosine-clustered chord x in [0, 1]
    VectorXd xi_eval_;  // sqrt(x_eval_)

    VectorXd x_camber_knots_;      // [n_camber_cp_]
    VectorXd xi_thickness_knots_;  // total thickness knots incl. fixed LE/TE

    MatrixXd B_c_;   // [N_per_side + 1, n_camber_cp_]    y_c(x_eval) = B_c·cp
    MatrixXd Bp_c_;  // [N_per_side + 1, n_camber_cp_]    dy_c/dx     = Bp_c·cp
    MatrixXd B_t_;   // [N_per_side + 1, n_thickness_cp_] t(x_eval)   = B_t·cp
};

}  // namespace scirocco
