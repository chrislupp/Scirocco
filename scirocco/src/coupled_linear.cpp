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
#include <cstddef>
#include <stdexcept>
#include <vector>

#ifdef SCIROCCO_HAVE_LAPACK
/* Fortran BLAS/LAPACK (LP64). The trailing arguments are the hidden lengths
 * of the character arguments; implementations that do not expect them
 * ignore them. */
extern "C" {
void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
            const int* k, const double* alpha, const double* a, const int* lda,
            const double* b, const int* ldb, const double* beta, double* c,
            const int* ldc, std::size_t, std::size_t);
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv,
             int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a,
             const int* lda, const int* ipiv, double* b, const int* ldb,
             int* info, std::size_t);
}
#endif

namespace scirocco {

namespace {

/// Relative residual accepted from a Schur solve before falling back.
constexpr double kSchurTolerance = 1e-10;

/// Iterative-refinement passes tried before falling back to a direct LU.
constexpr int kRefinementPasses = 3;

/// Tarjan's strongly connected components; components come out with every
/// dependency before its dependents (solve order).
class Components {
   public:
    explicit Components(const std::vector<std::vector<int>>& deps)
        : deps_(deps),
          index_(deps.size(), -1),
          low_(deps.size(), 0),
          on_stack_(deps.size(), false) {
        for (int v = 0; v < static_cast<int>(deps.size()); ++v)
            if (index_[v] < 0) Visit(v);
    }
    std::vector<std::vector<int>> components;

   private:
    void Visit(int v) {
        index_[v] = low_[v] = counter_++;
        stack_.push_back(v);
        on_stack_[v] = true;
        for (int w : deps_[v]) {
            if (index_[w] < 0) {
                Visit(w);
                low_[v] = std::min(low_[v], low_[w]);
            } else if (on_stack_[w]) {
                low_[v] = std::min(low_[v], index_[w]);
            }
        }
        if (low_[v] != index_[v]) return;
        std::vector<int> component;
        int w = -1;
        do {
            w = stack_.back();
            stack_.pop_back();
            on_stack_[w] = false;
            component.push_back(w);
        } while (w != v);
        std::sort(component.begin(), component.end());
        components.push_back(std::move(component));
    }
    const std::vector<std::vector<int>>& deps_;
    std::vector<int> index_, low_;
    std::vector<bool> on_stack_;
    std::vector<int> stack_;
    int counter_ = 0;
};

}  // namespace

bool BlockTriangularSolver::Compute(const SparseMatrixd& a, int block) {
    block_ = block;
    const int nodes = static_cast<int>(a.rows()) / block;
    std::vector<std::vector<int>> deps(nodes);
    for (int c = 0; c < a.outerSize(); ++c)
        for (SparseMatrixd::InnerIterator it(a, c); it; ++it) {
            const int i = static_cast<int>(it.row()) / block;
            const int j = c / block;
            if (i != j) deps[i].push_back(j);
        }
    for (auto& d : deps) {
        std::sort(d.begin(), d.end());
        d.erase(std::unique(d.begin(), d.end()), d.end());
    }
    const Components scc(deps);

    groups_.assign(scc.components.size(), Group{});
    std::vector<int> group_of(nodes), local(nodes);
    std::vector<MatrixXd> diag(groups_.size());
    for (size_t g = 0; g < scc.components.size(); ++g) {
        const auto& members = scc.components[g];
        for (size_t k = 0; k < members.size(); ++k) {
            group_of[members[k]] = static_cast<int>(g);
            local[members[k]] = static_cast<int>(k);
            for (int b = 0; b < block; ++b)
                groups_[g].rows.push_back(block * members[k] + b);
        }
        const int size = block * static_cast<int>(members.size());
        diag[g] = MatrixXd::Zero(size, size);
    }
    for (int c = 0; c < a.outerSize(); ++c)
        for (SparseMatrixd::InnerIterator it(a, c); it; ++it) {
            const int r = static_cast<int>(it.row());
            const int gi = group_of[r / block];
            const int gj = group_of[c / block];
            if (gi == gj) {
                diag[gi](block * local[r / block] + r % block,
                         block * local[c / block] + c % block) += it.value();
            } else {
                groups_[gi].coupling.push_back({r, c, it.value()});
            }
        }
    for (size_t g = 0; g < groups_.size(); ++g) {
        const MatrixXd& d = diag[g];
        groups_[g].lu.compute(d);
        const VectorXd pivots = groups_[g].lu.matrixLU().diagonal().cwiseAbs();
        const double norm = d.cwiseAbs().maxCoeff();
        if (!(pivots.minCoeff() > 1e-13 * norm)) return false;
    }
    return true;
}

void BlockTriangularSolver::SolveInPlace(RowMatrixXd& y) const {
    for (const Group& g : groups_) {
        for (const Entry& e : g.coupling)
            y.row(e.row) -= e.value * y.row(e.col);
        const int size = static_cast<int>(g.rows.size());
        if (size == block_) {
            /* A single node: its rows are contiguous. */
            y.middleRows(g.rows[0], size) =
                g.lu.solve(MatrixXd(y.middleRows(g.rows[0], size)));
        } else {
            MatrixXd local(size, y.cols());
            for (int k = 0; k < size; ++k) local.row(k) = y.row(g.rows[k]);
            local = g.lu.solve(local);
            for (int k = 0; k < size; ++k) y.row(g.rows[k]) = local.row(k);
        }
    }
}

VectorXd BlockTriangularSolver::Solve(const VectorXd& r) const {
    RowMatrixXd y = r;
    SolveInPlace(y);
    return y;
}

VectorXd BlockTriangularSolver::SolveTranspose(const VectorXd& r) const {
    /* A^T is block upper triangular in solve order: substitute backwards,
     * then push each group's result into the rows it couples to. */
    VectorXd t = r;
    VectorXd z = VectorXd::Zero(r.size());
    for (auto g = groups_.rbegin(); g != groups_.rend(); ++g) {
        const int size = static_cast<int>(g->rows.size());
        VectorXd local(size);
        for (int k = 0; k < size; ++k) local(k) = t(g->rows[k]);
        local = g->lu.transpose().solve(local);
        for (int k = 0; k < size; ++k) z(g->rows[k]) = local(k);
        for (const Entry& e : g->coupling) t(e.col) -= e.value * z(e.row);
    }
    return z;
}

int BlockTriangularSolver::LargestGroup() const {
    size_t largest = 0;
    for (const Group& g : groups_) largest = std::max(largest, g.rows.size());
    return static_cast<int>(largest) / block_;
}

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

    if (!b_solver_.Compute(b, 3)) return false;
    x_ = RowMatrixXd(g);
    b_solver_.SolveInPlace(x_);
    if (!x_.allFinite()) return false;

    /* S = D + ue_m X_m, X_m the mass rows of X. */
    MatrixXd xm(n, n);
    for (int j = 0; j < n; ++j) xm.row(j) = x_.row(3 * j + 1);
    MatrixXd s(n, n);
#ifdef SCIROCCO_HAVE_LAPACK
    {
        const int nn = n;
        const double one = 1.0, zero = 0.0;
        dgemm_("N", "N", &nn, &nn, &nn, &one, ue_m_.data(), &nn, xm.data(), &nn,
               &zero, s.data(), &nn, 1, 1);
    }
#else
    s.noalias() = ue_m_ * xm;
#endif
    s.diagonal() += d_;
    if (!s.allFinite()) return false;
    return FactorSchurComplement(std::move(s));
}

bool CoupledFactorization::FactorSchurComplement(MatrixXd s) {
#ifdef SCIROCCO_HAVE_LAPACK
    const int n = static_cast<int>(s.rows());
    s_factors_ = std::move(s);
    s_pivots_.assign(n, 0);
    int info = 0;
    dgetrf_(&n, &n, s_factors_.data(), &n, s_pivots_.data(), &info);
    if (info != 0) return false;
    const VectorXd pivots = s_factors_.diagonal().cwiseAbs();
#else
    s_lu_.compute(s);
    /* PartialPivLU does not report singularity; check the pivots. */
    const VectorXd pivots = s_lu_.matrixLU().diagonal().cwiseAbs();
#endif
    return pivots.minCoeff() > 1e-14 * pivots.maxCoeff();
}

VectorXd CoupledFactorization::SolveSchurComplement(const VectorXd& r,
                                                    bool transpose) const {
#ifdef SCIROCCO_HAVE_LAPACK
    VectorXd x = r;
    const int n = static_cast<int>(x.size());
    const int one = 1;
    int info = 0;
    dgetrs_(transpose ? "T" : "N", &n, &one, s_factors_.data(), &n,
            s_pivots_.data(), x.data(), &n, &info, 1);
    return x;
#else
    return transpose ? VectorXd(s_lu_.transpose().solve(r))
                     : VectorXd(s_lu_.solve(r));
#endif
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
        b_solver_.Solve(row_scale_.cwiseProduct(rhs.head(3 * n)));
    VectorXd y0m(n);
    for (int j = 0; j < n; ++j) y0m(j) = y0(3 * j + 1);
    const VectorXd u = SolveSchurComplement(rhs.tail(n) + ue_m_ * y0m, false);
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
    const VectorXd z2 = SolveSchurComplement(h2 - x_.transpose() * h1, true);
    /* B^T z1 = h1 + E_m ue_m^T z2; B = R^-1 B~ so B^-T = R B~^-T. */
    const VectorXd q = ue_m_.transpose() * z2;
    VectorXd t = h1;
    for (int j = 0; j < n; ++j) t(3 * j + 1) += q(j);
    const VectorXd z1 = row_scale_.cwiseProduct(b_solver_.SolveTranspose(t));
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
