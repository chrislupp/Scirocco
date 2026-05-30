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

    Scirocco Geometry Discipline Implementation
*/

#include "scirocco/geometry_discipline.h"

#include <Eigen/Core>
#include <cmath>
#include <stdexcept>

#include "scirocco/core/spline.h"

namespace scirocco {

void GeometryDiscipline::Initialize() {
    AddOption("n_camber_cp", "int");
    AddOption("n_thickness_cp", "int");
    AddOption("npanel", "int");
    AddOption("closed_te", "bool");
}

void GeometryDiscipline::SetOptions(const google::protobuf::Struct& options) {
    const auto& fields = options.fields();

    if (fields.count("n_camber_cp")) {
        n_camber_cp_ =
            static_cast<int>(fields.at("n_camber_cp").number_value());
    }
    if (fields.count("n_thickness_cp")) {
        n_thickness_cp_ =
            static_cast<int>(fields.at("n_thickness_cp").number_value());
    }
    if (fields.count("npanel")) {
        npanel_ = static_cast<int>(fields.at("npanel").number_value());
    }
    if (fields.count("closed_te")) {
        closed_te_ = fields.at("closed_te").bool_value();
    }

    if (npanel_ < 4 || npanel_ % 2 != 0) {
        throw std::invalid_argument(
            "GeometryDiscipline: npanel must be an even integer >= 4");
    }
    if (n_camber_cp_ < 4) {
        throw std::invalid_argument(
            "GeometryDiscipline: n_camber_cp must be >= 4");
    }
    if (n_thickness_cp_ < 2) {
        throw std::invalid_argument(
            "GeometryDiscipline: n_thickness_cp must be >= 2");
    }

    N_per_side_ = npanel_ / 2;
    npoint_ = npanel_ + 1;

    BuildEvalStations();
    BuildBasisMatrices();
}

void GeometryDiscipline::Configure() {}

void GeometryDiscipline::Setup() {
    AddInput("camber_cp", {n_camber_cp_}, "");
    AddInput("thickness_cp", {n_thickness_cp_}, "");
    AddOutput("x", {2 * npoint_}, "m");
}

void GeometryDiscipline::SetupPartials() {
    DeclarePartials("x", "camber_cp");
    DeclarePartials("x", "thickness_cp");
}

void GeometryDiscipline::BuildEvalStations() {
    // Same TE-bunching distribution used by NacaPoints — gives smooth LE/TE
    // clustering and lets symmetric configurations match the NACA generator
    // closely.
    const int N = N_per_side_;
    const double te = 1.5;
    x_eval_.resize(N + 1);
    xi_eval_.resize(N + 1);

    VectorXd f = VectorXd::LinSpaced(N + 1, 0.0, 1.0);
    for (int i = 0; i <= N; ++i) {
        const double fi = f(i);
        double xi = 1.0 - (te + 1.0) * fi * std::pow(1.0 - fi, te) -
                    std::pow(1.0 - fi, te + 1.0);
        if (xi < 0.0) xi = 0.0;  // guard fp roundoff
        if (xi > 1.0) xi = 1.0;
        x_eval_(i) = xi;
        xi_eval_(i) = std::sqrt(xi);
    }
}

void GeometryDiscipline::BuildBasisMatrices() {
    const int N = N_per_side_;

    // Camber knots: cosine-clustered on [0, 1].
    x_camber_knots_.resize(n_camber_cp_);
    for (int j = 0; j < n_camber_cp_; ++j) {
        x_camber_knots_(j) =
            0.5 * (1.0 - std::cos(M_PI * j / (n_camber_cp_ - 1)));
    }

    // Cardinal-basis matrices via one-hot CP vectors. The cubic spline is
    // linear in its data values, so col-j of B_c / Bp_c is the spline /
    // spline-derivative evaluated at x_eval with the j-th CP = 1, others = 0.
    B_c_.resize(N + 1, n_camber_cp_);
    Bp_c_.resize(N + 1, n_camber_cp_);
    for (int j = 0; j < n_camber_cp_; ++j) {
        VectorXd y_hot = VectorXd::Zero(n_camber_cp_);
        y_hot(j) = 1.0;
        CubicSpline1D spline = Spline1DCreate(x_camber_knots_, y_hot);
        for (int q = 0; q <= N; ++q) {
            B_c_(q, j) = spline.Eval(x_eval_(q));
            Bp_c_(q, j) = spline.EvalDeriv(x_eval_(q));
        }
    }

    // Thickness spline knots: one fixed LE knot (xi=0, t=0), n_thickness_cp_
    // free knots, optionally one fixed TE knot (xi=1, t=0).
    const int n_total = n_thickness_cp_ + 1 + (closed_te_ ? 1 : 0);
    xi_thickness_knots_.resize(n_total);
    for (int j = 0; j < n_total; ++j) {
        xi_thickness_knots_(j) =
            0.5 * (1.0 - std::cos(M_PI * j / (n_total - 1)));
    }

    // Free CP k corresponds to full-knot index 1+k (LE skipped).
    B_t_.resize(N + 1, n_thickness_cp_);
    for (int k = 0; k < n_thickness_cp_; ++k) {
        VectorXd t_hot = VectorXd::Zero(n_total);
        t_hot(1 + k) = 1.0;
        CubicSpline1D spline = Spline1DCreate(xi_thickness_knots_, t_hot);
        for (int q = 0; q <= N; ++q) {
            B_t_(q, k) = spline.Eval(xi_eval_(q));
        }
    }
}

void GeometryDiscipline::ComputeSurface(const VectorXd& camber_cp,
                                        const VectorXd& thickness_cp,
                                        MatrixXd& xpoint) const {
    const int N = N_per_side_;
    xpoint.resize(2, npoint_);

    const VectorXd yc = B_c_ * camber_cp;
    const VectorXd dyc = Bp_c_ * camber_cp;
    const VectorXd tt = B_t_ * thickness_cp;

    // Sirocco order: indices [0..N] = lower surface from TE to LE,
    // index N = LE (shared), indices [N+1..2N] = upper surface from LE+1 to TE.
    for (int q = 0; q <= N; ++q) {
        const double n_inv = 1.0 / std::sqrt(1.0 + dyc(q) * dyc(q));
        const double tdn = tt(q) * dyc(q) * n_inv;  // t · dy · n
        const double tn = tt(q) * n_inv;            // t · n

        const double Pu_x = x_eval_(q) - tdn;
        const double Pu_y = yc(q) + tn;
        const double Pl_x = x_eval_(q) + tdn;
        const double Pl_y = yc(q) - tn;

        // Lower (reversed): out index = N - q
        const int out_l = N - q;
        xpoint(0, out_l) = Pl_x;
        xpoint(1, out_l) = Pl_y;

        // Upper: out index = N + q. Skip q = 0 since the LE point is the
        // shared lower endpoint at index N (already written above).
        if (q > 0) {
            const int out_u = N + q;
            xpoint(0, out_u) = Pu_x;
            xpoint(1, out_u) = Pu_y;
        }
    }
}

void GeometryDiscipline::Compute(const philote::Variables& inputs,
                                 philote::Variables& outputs) {
    const auto& camber_var = inputs.at("camber_cp");
    const auto& thickness_var = inputs.at("thickness_cp");

    VectorXd camber_cp(n_camber_cp_);
    for (int i = 0; i < n_camber_cp_; ++i) camber_cp(i) = camber_var(i);
    VectorXd thickness_cp(n_thickness_cp_);
    for (int i = 0; i < n_thickness_cp_; ++i)
        thickness_cp(i) = thickness_var(i);

    MatrixXd xpoint;
    ComputeSurface(camber_cp, thickness_cp, xpoint);

    auto& x_var = outputs.at("x");
    for (int i = 0; i < npoint_; ++i) {
        x_var(i) = xpoint(0, i);
        x_var(npoint_ + i) = xpoint(1, i);
    }
}

void GeometryDiscipline::ComputePartials(const philote::Variables& inputs,
                                         philote::Partials& partials) {
    const auto& camber_var = inputs.at("camber_cp");
    const auto& thickness_var = inputs.at("thickness_cp");

    VectorXd camber_cp(n_camber_cp_);
    for (int i = 0; i < n_camber_cp_; ++i) camber_cp(i) = camber_var(i);
    VectorXd thickness_cp(n_thickness_cp_);
    for (int i = 0; i < n_thickness_cp_; ++i)
        thickness_cp(i) = thickness_var(i);

    const int N = N_per_side_;
    const VectorXd dyc = Bp_c_ * camber_cp;
    const VectorXd tt = B_t_ * thickness_cp;

    auto& d_c = partials[{"x", "camber_cp"}];
    auto& d_t = partials[{"x", "thickness_cp"}];

    // Output layout: out_idx ∈ [0, npoint). Flat row 'row' for x-coord at
    // out_idx is row = out_idx; for y-coord it is row = npoint_ + out_idx.
    // Partial entry (row, col) lives at flat index row * n_cp + col.
    const int n_rows = 2 * npoint_;
    for (int row = 0; row < n_rows; ++row) {
        for (int j = 0; j < n_camber_cp_; ++j)
            d_c(row * n_camber_cp_ + j) = 0.0;
        for (int k = 0; k < n_thickness_cp_; ++k)
            d_t(row * n_thickness_cp_ + k) = 0.0;
    }

    for (int q = 0; q <= N; ++q) {
        const double dy = dyc(q);
        const double tq = tt(q);
        const double s2 = 1.0 + dy * dy;
        const double n_inv = 1.0 / std::sqrt(s2);   // n
        const double n_inv3 = n_inv / s2;           // n^3

        // tdn = t·dy·n  →  d(tdn)/d(dy) = t · n^3,  d(tdn)/d(t) = dy·n
        // tn  = t·n     →  d(tn)/d(dy)  = -t·dy·n^3, d(tn)/d(t) = n
        const double d_tdn_ddy = tq * n_inv3;
        const double d_tdn_dt = dy * n_inv;
        const double d_tn_ddy = -tq * dy * n_inv3;
        const double d_tn_dt = n_inv;

        // Upper / lower output indices (lower reversed).
        const int out_l = N - q;
        const int out_u = N + q;  // valid only for q > 0
        const int xrow_l = out_l;
        const int yrow_l = npoint_ + out_l;
        const int xrow_u = out_u;
        const int yrow_u = npoint_ + out_u;

        // ∂/∂camber_cp_j contributions (via yc and dyc):
        //   Pu.x = x_q - tdn       → dPu.x/dcp = -d_tdn_ddy * Bp_c(q,j)
        //   Pu.y = yc + tn         → dPu.y/dcp =  B_c(q,j) + d_tn_ddy * Bp_c(q,j)
        //   Pl.x = x_q + tdn       → dPl.x/dcp =  d_tdn_ddy * Bp_c(q,j)
        //   Pl.y = yc - tn         → dPl.y/dcp =  B_c(q,j) - d_tn_ddy * Bp_c(q,j)
        for (int j = 0; j < n_camber_cp_; ++j) {
            const double Bcj = B_c_(q, j);
            const double Bpcj = Bp_c_(q, j);
            d_c(xrow_l * n_camber_cp_ + j) = d_tdn_ddy * Bpcj;
            d_c(yrow_l * n_camber_cp_ + j) = Bcj - d_tn_ddy * Bpcj;
            if (q > 0) {
                d_c(xrow_u * n_camber_cp_ + j) = -d_tdn_ddy * Bpcj;
                d_c(yrow_u * n_camber_cp_ + j) = Bcj + d_tn_ddy * Bpcj;
            }
        }

        // ∂/∂thickness_cp_k contributions (via t):
        //   dPu.x/dcp = -d_tdn_dt * B_t(q,k)
        //   dPu.y/dcp =  d_tn_dt  * B_t(q,k)
        //   dPl.x/dcp =  d_tdn_dt * B_t(q,k)
        //   dPl.y/dcp = -d_tn_dt  * B_t(q,k)
        for (int k = 0; k < n_thickness_cp_; ++k) {
            const double Btk = B_t_(q, k);
            d_t(xrow_l * n_thickness_cp_ + k) = d_tdn_dt * Btk;
            d_t(yrow_l * n_thickness_cp_ + k) = -d_tn_dt * Btk;
            if (q > 0) {
                d_t(xrow_u * n_thickness_cp_ + k) = -d_tdn_dt * Btk;
                d_t(yrow_u * n_thickness_cp_ + k) = d_tn_dt * Btk;
            }
        }
    }
}

}  // namespace scirocco
