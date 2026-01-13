// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CONFIG_H
#define CPPFOIL_CONFIG_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace cppfoil {

// Eigen type aliases
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using Vector2d = Eigen::Vector2d;
using Vector3d = Eigen::Vector3d;
using Vector4d = Eigen::Vector4d;
using Matrix2d = Eigen::Matrix2d;
using Matrix3d = Eigen::Matrix3d;
using Matrix4d = Eigen::Matrix4d;

// Fixed-size types for boundary layer computations
using Vector8d = Eigen::Matrix<double, 8, 1>;
using Matrix34d = Eigen::Matrix<double, 3, 4>;
using Matrix38d = Eigen::Matrix<double, 3, 8>;
using Matrix32d = Eigen::Matrix<double, 3, 2>;
using Matrix42d = Eigen::Matrix<double, 4, 2>;
using Matrix48d = Eigen::Matrix<double, 4, 8>;
using Matrix312d = Eigen::Matrix<double, 3, 12>;

// Sparse matrix types
using SparseMatrixd = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

// Constants
constexpr double kPi = 3.14159265358979323846;
constexpr double kDegToRad = kPi / 180.0;
constexpr double kRadToDeg = 180.0 / kPi;

}  // namespace cppfoil

#endif  // CPPFOIL_CONFIG_H
