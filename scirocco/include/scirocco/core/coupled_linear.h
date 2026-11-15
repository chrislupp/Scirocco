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
 * @file coupled_linear.h
 * @brief Factorization of the coupled viscous-inviscid Jacobian.
 *
 * The coupled Jacobian J has 4N columns, [th, ds, sa, ue] at each node, and
 * 4N rows: 3N boundary-layer rows (the sparse, marching R_U) followed by N
 * edge-velocity rows
 *
 *     d_i due_i - sum_j ue_m(i, j) (ds_j due_j + ue_j dds_j).
 *
 * The edge-velocity rows couple every node to every other through the dense
 * mass-influence matrix ue_m, so a general sparse LU of J fills in an N x 2N
 * dense block. Instead, substitute the mass defect m = ds ue for ds. The
 * edge-velocity rows then read d_i due_i - sum_j ue_m(i, j) dm_j, with an
 * identity-like block in ue. Eliminating the boundary-layer unknowns
 * [th, m, sa] through the sparse block B (the boundary layer at prescribed
 * mass, i.e. inverse mode, which stays regular through separation) leaves the
 * N x N dense Schur complement
 *
 *     S = diag(d) + ue_m (B^-1 G)_m,
 *
 * where G holds the boundary-layer sensitivities to ue at fixed mass and the
 * subscript m selects the mass rows. Both J x = r and J^T x = r reduce to
 * solves with B and S. The result is exact; it is the same Newton step and
 * the same adjoint as a direct factorization of J.
 *
 * If B cannot be factored, or a solve fails its residual check, the
 * factorization falls back to a sparse LU of the assembled J.
 */
#pragma once

#include <Eigen/SparseLU>
#include <memory>
#include <vector>

#include "scirocco/config.h"

namespace scirocco {

/** Row-major dense matrix, for row-wise block substitution. */
using RowMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/**
 * @brief Direct solver for a sparse matrix that is block triangular up to a
 * permutation of its node blocks.
 *
 * The boundary-layer Jacobian at prescribed mass couples each station to the
 * station upstream, except at the stagnation point (two stations solved
 * together) and the first wake station (both trailing edges). Grouping the
 * nodes into strongly connected components of the dependency graph and
 * ordering the groups topologically makes the matrix block lower
 * triangular, so solves are a forward (or, transposed, backward)
 * substitution with small dense LUs on the diagonal groups.
 */
class BlockTriangularSolver {
   public:
    /**
     * @brief Factor @p a.
     * @param a Square sparse matrix
     * @param block Rows and columns per node
     * @return false if a diagonal group is singular
     */
    bool Compute(const SparseMatrixd& a, int block);

    /** @brief Overwrite @p y with A^-1 y (any number of columns). */
    void SolveInPlace(RowMatrixXd& y) const;

    /** @brief A^-1 r. */
    [[nodiscard]] VectorXd Solve(const VectorXd& r) const;

    /** @brief A^-T r. */
    [[nodiscard]] VectorXd SolveTranspose(const VectorXd& r) const;

   private:
    struct Entry {
        int row;
        int col;
        double value;
    };
    struct Group {
        std::vector<int> rows; /* Matrix rows (and columns) of the group */
        Eigen::PartialPivLU<MatrixXd> lu;
        std::vector<Entry> coupling; /* Rows in the group, columns earlier */
    };
    std::vector<Group> groups_; /* In solve order */
    int block_ = 1;
};

/**
 * @brief Factorization of the coupled viscous-inviscid Jacobian J through
 * the mass-defect Schur complement. See the file comment for the method.
 */
class CoupledFactorization {
   public:
    /**
     * @brief Factor J.
     * @param bl_jacobian Boundary-layer rows of J, [3N x 4N]
     * @param ue_m Mass-influence matrix, d(ue)/d(mass) [N x N]
     * @param ue Edge velocity used in the edge-velocity rows [N]
     * @param ds Displacement thickness [N]
     * @param ue_diagonal Diagonal coefficient d_i of due_i in the
     *        edge-velocity rows (1 without pseudo-time terms) [N]
     * @return true if J was factored (by either method)
     */
    bool Compute(const SparseMatrixd& bl_jacobian, const MatrixXd& ue_m,
                 const VectorXd& ue, const VectorXd& ds,
                 const VectorXd& ue_diagonal);

    /** @brief Solve J x = rhs. Throws std::runtime_error on failure. */
    VectorXd Solve(const VectorXd& rhs);

    /** @brief Solve J^T x = rhs. Throws std::runtime_error on failure. */
    VectorXd SolveTranspose(const VectorXd& rhs);

    /** @brief J v. */
    [[nodiscard]] VectorXd Apply(const VectorXd& v) const;

    /** @brief J^T w. */
    [[nodiscard]] VectorXd ApplyTranspose(const VectorXd& w) const;

    /** @brief An upper bound on |J|^T |w|, for residual scales. */
    [[nodiscard]] VectorXd AbsApplyTranspose(const VectorXd& w) const;

    /** @brief Assemble J as a sparse matrix (dense edge-velocity rows). */
    [[nodiscard]] SparseMatrixd Assemble() const;

    /** @brief True if solves go through the Schur complement. */
    [[nodiscard]] bool UsesSchur() const { return use_schur_; }

   private:
    bool ComputeSchur();
    bool ComputeDirect();
    VectorXd SolveEither(const VectorXd& rhs, bool transpose);
    [[nodiscard]] VectorXd SolveSchur(const VectorXd& rhs) const;
    [[nodiscard]] VectorXd SolveTransposeSchur(const VectorXd& rhs) const;
    [[nodiscard]] bool Accurate(const VectorXd& x, const VectorXd& rhs,
                                bool transpose) const;
    bool FactorSchurComplement(const MatrixXd& s);
    [[nodiscard]] VectorXd SolveSchurComplement(const VectorXd& r,
                                                bool transpose) const;

    int n_ = 0;
    SparseMatrixd bl_;
    MatrixXd ue_m_;
    VectorXd ue_, ds_, d_;

    bool use_schur_ = false;
    VectorXd row_scale_; /* Row equilibration of B and G */
    BlockTriangularSolver b_solver_;
    RowMatrixXd x_; /* B^-1 G [3N x N] */
    /* LU of the Schur complement S: LAPACK factors and pivots when built
       with LAPACK, Eigen's otherwise */
    MatrixXd s_factors_;
    std::vector<int> s_pivots_;
    Eigen::PartialPivLU<MatrixXd> s_lu_;

    mutable std::unique_ptr<Eigen::SparseLU<SparseMatrixd>> direct_lu_;
    VectorXd direct_row_scale_;
};

}  // namespace scirocco
