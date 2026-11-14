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
*/

/**
 * @file coupled_linear.cpp
 * @brief Mass-defect Schur complement factorization of the coupled Jacobian.
 *
 * Notation (n = number of nodes): J = [A; L] with A the boundary-layer rows
 * (3n x 4n) and L the edge-velocity rows (n x 4n). The column transform T
 * maps w = [th, m, sa, u] to v = [th, ds, sa, ue] with
 *     dds = (dm - ds du) / ue,   due = du.
 * K = J T has the block form
 *     K = [B  G]    B: A on [th, m, sa]  (3n x 3n, sparse)
 *         [L' D]    G: A on u            (3n x n, sparse)
 * with L' = -ue_m on the mass columns and D = diag(d). Then
 *     S = D - L' B^-1 G = D + ue_m (B^-1 G)_m.
 */
#include "scirocco/core/coupled_linear.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace scirocco {

namespace {

/// Relative residual accepted from a Schur solve before falling back.
constexpr double kSchurTolerance = 1e-10;

/// Iterative-refinement passes tried before falling back to a direct LU.
constexpr int kRefinementPasses = 3;

}  // namespace

bool CoupledFactorization::Compute(const SparseMatrixd& bl_jacobian,
                                   const MatrixXd& ue_m, const VectorXd& ue,
                                   const VectorXd& ds,
                                   const VectorXd& ue_diagonal) {
    n_ = static_cast<int>(ue.size());
    if (bl_jacobian.rows() != 3 * n_ || bl_jacobian.cols() != 4 * n_ ||
        ue_m.rows() != n_ || ue_m.cols() != n_ || ds.size() != n_ ||
        ue_diagonal.size() != n_) {
        throw std::invalid_argument(
            "CoupledFactorization: inconsistent system dimensions");
    }
    bl_ = bl_jacobian;
    bl_.makeCompressed();
    ue_m_ = ue_m;
    ue_ = ue;
    ds_ = ds;
    d_ = ue_diagonal;
    direct_lu_.reset();
    use_schur_ = ComputeSchur();
    if (use_schur_) return true;
    return ComputeDirect();
}

bool CoupledFactorization::ComputeSchur() {
    const int n = n_;
    if ((ue_.array() == 0.0).any()) return false;

    /* Split A into B (th, m, sa columns) and G (u columns). */
    std::vector<Triplet> tb, tg;
    tb.reserve(bl_.nonZeros());
    tg.reserve(bl_.nonZeros() / 2);
    for (int c = 0; c < bl_.outerSize(); ++c) {
        const int j = c / 4;
        const int var = c % 4;
        for (SparseMatrixd::InnerIterator it(bl_, c); it; ++it) {
            const int r = static_cast<int>(it.row());
            const double a = it.value();
            switch (var) {
                case 0:
                    tb.emplace_back(r, 3 * j, a);
                    break;
                case 1:
                    tb.emplace_back(r, 3 * j + 1, a / ue_(j));
                    tg.emplace_back(r, j, -a * ds_(j) / ue_(j));
                    break;
                case 2:
                    tb.emplace_back(r, 3 * j + 2, a);
                    break;
                default:
                    tg.emplace_back(r, j, a);
                    break;
            }
        }
    }
    SparseMatrixd b(3 * n, 3 * n), g(3 * n, n);
    b.setFromTriplets(tb.begin(), tb.end());
    g.setFromTriplets(tg.begin(), tg.end());

    /* Row equilibration: partial pivoting in SparseLU compares entries
     * within a column, so mixed equation magnitudes matter. It does not
     * change B^-1 G. */
    row_scale_ = VectorXd::Zero(3 * n);
    for (int c = 0; c < b.outerSize(); ++c)
        for (SparseMatrixd::InnerIterator it(b, c); it; ++it)
            row_scale_(it.row()) =
                std::max(row_scale_(it.row()), std::abs(it.value()));
    for (int r = 0; r < 3 * n; ++r) {
        if (!(row_scale_(r) > 0.0) || !std::isfinite(row_scale_(r)))
            return false;
        row_scale_(r) = 1.0 / row_scale_(r);
    }
    b = row_scale_.asDiagonal() * b;
    g = row_scale_.asDiagonal() * g;

    b_lu_.compute(b);
    if (b_lu_.info() != Eigen::Success) return false;

    x_ = b_lu_.solve(MatrixXd(g));
    if (b_lu_.info() != Eigen::Success || !x_.allFinite()) return false;

    /* S = D + ue_m X_m, X_m the mass rows of X. */
    MatrixXd xm(n, n);
    for (int j = 0; j < n; ++j) xm.row(j) = x_.row(3 * j + 1);
    MatrixXd s = ue_m_ * xm;
    s.diagonal() += d_;
    if (!s.allFinite()) return false;
    s_lu_.compute(s);
    /* PartialPivLU does not report singularity; check the pivots. */
    const VectorXd pivots = s_lu_.matrixLU().diagonal().cwiseAbs();
    if (!(pivots.minCoeff() > 1e-14 * pivots.maxCoeff())) return false;
    return true;
}

bool CoupledFactorization::ComputeDirect() {
    SparseMatrixd j = Assemble();
    direct_row_scale_ = VectorXd::Zero(j.rows());
    for (int c = 0; c < j.outerSize(); ++c)
        for (SparseMatrixd::InnerIterator it(j, c); it; ++it)
            direct_row_scale_(it.row()) =
                std::max(direct_row_scale_(it.row()), std::abs(it.value()));
    for (int r = 0; r < j.rows(); ++r)
        direct_row_scale_(r) =
            direct_row_scale_(r) > 0.0 ? 1.0 / direct_row_scale_(r) : 1.0;
    j = direct_row_scale_.asDiagonal() * j;
    direct_lu_ = std::make_unique<Eigen::SparseLU<SparseMatrixd>>();
    direct_lu_->compute(j);
    if (direct_lu_->info() != Eigen::Success) {
        direct_lu_.reset();
        return false;
    }
    return true;
}

SparseMatrixd CoupledFactorization::Assemble() const {
    const int n = n_;
    std::vector<Triplet> t;
    t.reserve(bl_.nonZeros() + 2 * static_cast<size_t>(n) * n + n);
    for (int c = 0; c < bl_.outerSize(); ++c)
        for (SparseMatrixd::InnerIterator it(bl_, c); it; ++it)
            t.emplace_back(it.row(), c, it.value());
    for (int i = 0; i < n; ++i) {
        const int row = 3 * n + i;
        t.emplace_back(row, 4 * i + 3, d_(i));
        for (int j = 0; j < n; ++j) {
            const double a = ue_m_(i, j);
            if (a == 0.0) continue;
            t.emplace_back(row, 4 * j + 3, -a * ds_(j));
            t.emplace_back(row, 4 * j + 1, -a * ue_(j));
        }
    }
    SparseMatrixd j(4 * n, 4 * n);
    j.setFromTriplets(t.begin(), t.end());
    return j;
}

VectorXd CoupledFactorization::Apply(const VectorXd& v) const {
    const int n = n_;
    VectorXd out(4 * n);
    out.head(3 * n) = bl_ * v;
    VectorXd mass(n), uev(n);
    for (int j = 0; j < n; ++j) {
        uev(j) = v(4 * j + 3);
        mass(j) = ds_(j) * uev(j) + ue_(j) * v(4 * j + 1);
    }
    out.tail(n) = d_.cwiseProduct(uev) - ue_m_ * mass;
    return out;
}

VectorXd CoupledFactorization::ApplyTranspose(const VectorXd& w) const {
    const int n = n_;
    VectorXd out = bl_.transpose() * w.head(3 * n);
    const VectorXd wu = w.tail(n);
    const VectorXd q = ue_m_.transpose() * wu;
    for (int j = 0; j < n; ++j) {
        out(4 * j + 1) -= ue_(j) * q(j);
        out(4 * j + 3) += d_(j) * wu(j) - ds_(j) * q(j);
    }
    return out;
}

VectorXd CoupledFactorization::AbsApplyTranspose(const VectorXd& w) const {
    const int n = n_;
    VectorXd out = bl_.cwiseAbs().transpose() * w.head(3 * n).cwiseAbs();
    const VectorXd wu = w.tail(n).cwiseAbs();
    const VectorXd q = ue_m_.cwiseAbs().transpose() * wu;
    for (int j = 0; j < n; ++j) {
        out(4 * j + 1) += std::abs(ue_(j)) * q(j);
        out(4 * j + 3) += std::abs(d_(j)) * wu(j) + std::abs(ds_(j)) * q(j);
    }
    return out;
}

VectorXd CoupledFactorization::SolveSchur(const VectorXd& rhs) const {
    const int n = n_;
    const VectorXd y0 =
        b_lu_.solve(row_scale_.cwiseProduct(rhs.head(3 * n)).eval());
    VectorXd y0m(n);
    for (int j = 0; j < n; ++j) y0m(j) = y0(3 * j + 1);
    const VectorXd u = s_lu_.solve(rhs.tail(n) + ue_m_ * y0m);
    const VectorXd y = y0 - x_ * u;
    VectorXd v(4 * n);
    for (int j = 0; j < n; ++j) {
        v(4 * j) = y(3 * j);
        v(4 * j + 1) = (y(3 * j + 1) - ds_(j) * u(j)) / ue_(j);
        v(4 * j + 2) = y(3 * j + 2);
        v(4 * j + 3) = u(j);
    }
    return v;
}

VectorXd CoupledFactorization::SolveTransposeSchur(const VectorXd& rhs) const {
    /* J^T z = g  <=>  K^T z = T^T g, with K^T = [B^T L'^T; G^T D]. */
    const int n = n_;
    VectorXd h1(3 * n), h2(n);
    for (int j = 0; j < n; ++j) {
        h1(3 * j) = rhs(4 * j);
        h1(3 * j + 1) = rhs(4 * j + 1) / ue_(j);
        h1(3 * j + 2) = rhs(4 * j + 2);
        h2(j) = rhs(4 * j + 3) - rhs(4 * j + 1) * ds_(j) / ue_(j);
    }
    /* S^T z2 = h2 - X^T h1 */
    const VectorXd z2 = s_lu_.transpose().solve(h2 - x_.transpose() * h1);
    /* B^T z1 = h1 + E_m ue_m^T z2; B = R^-1 B~ so B^-T = R B~^-T. */
    const VectorXd q = ue_m_.transpose() * z2;
    VectorXd t = h1;
    for (int j = 0; j < n; ++j) t(3 * j + 1) += q(j);
    const VectorXd z1 = row_scale_.cwiseProduct(b_lu_.transpose().solve(t));
    VectorXd z(4 * n);
    z.head(3 * n) = z1;
    z.tail(n) = z2;
    return z;
}

bool CoupledFactorization::Accurate(const VectorXd& x, const VectorXd& rhs,
                                    bool transpose) const {
    if (!x.allFinite()) return false;
    const VectorXd r = (transpose ? ApplyTranspose(x) : Apply(x)) - rhs;
    double scale = rhs.norm();
    if (transpose) {
        scale += AbsApplyTranspose(x).norm();
    } else {
        VectorXd ax(4 * n_);
        ax.head(3 * n_) = bl_.cwiseAbs() * x.cwiseAbs();
        VectorXd mass(n_);
        for (int j = 0; j < n_; ++j)
            mass(j) = std::abs(ds_(j) * x(4 * j + 3)) +
                      std::abs(ue_(j) * x(4 * j + 1));
        ax.tail(n_) = d_.cwiseAbs().cwiseProduct(x.tail(n_).cwiseAbs()) +
                      ue_m_.cwiseAbs() * mass;
        scale += ax.norm();
    }
    return r.allFinite() &&
           r.norm() <= kSchurTolerance * std::max(scale, 1e-300);
}

VectorXd CoupledFactorization::Solve(const VectorXd& rhs) {
    if (rhs.size() != 4 * n_)
        throw std::invalid_argument("CoupledFactorization: rhs size");
    if (use_schur_) {
        /* Iterative refinement recovers accuracy when B is poorly
         * conditioned (rough intermediate Newton iterates); each pass costs
         * O(N^2). */
        VectorXd x = SolveSchur(rhs);
        for (int pass = 0;; ++pass) {
            if (Accurate(x, rhs, false)) return x;
            if (pass == kRefinementPasses) break;
            x += SolveSchur(rhs - Apply(x));
        }
        use_schur_ = false;
        if (!ComputeDirect())
            throw std::runtime_error(
                "CoupledFactorization: direct factorization failed");
    }
    if (!direct_lu_)
        throw std::runtime_error("CoupledFactorization: not factored");
    VectorXd x = direct_lu_->solve(direct_row_scale_.cwiseProduct(rhs));
    if (direct_lu_->info() != Eigen::Success || !x.allFinite())
        throw std::runtime_error("CoupledFactorization: solve failed");
    return x;
}

VectorXd CoupledFactorization::SolveTranspose(const VectorXd& rhs) {
    if (rhs.size() != 4 * n_)
        throw std::invalid_argument("CoupledFactorization: rhs size");
    if (use_schur_) {
        VectorXd x = SolveTransposeSchur(rhs);
        for (int pass = 0;; ++pass) {
            if (Accurate(x, rhs, true)) return x;
            if (pass == kRefinementPasses) break;
            x += SolveTransposeSchur(rhs - ApplyTranspose(x));
        }
        use_schur_ = false;
        if (!ComputeDirect())
            throw std::runtime_error(
                "CoupledFactorization: direct factorization failed");
    }
    if (!direct_lu_)
        throw std::runtime_error("CoupledFactorization: not factored");
    /* (R J)^T y = g with x = R y. */
    VectorXd y = direct_lu_->transpose().solve(rhs);
    if (!y.allFinite())
        throw std::runtime_error("CoupledFactorization: solve failed");
    return direct_row_scale_.cwiseProduct(y);
}

}  // namespace scirocco
