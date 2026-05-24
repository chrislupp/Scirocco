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

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#ifndef SCIROCCO_CONFIG_H
#define SCIROCCO_CONFIG_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace scirocco {

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

}  // namespace scirocco

#endif  // SCIROCCO_CONFIG_H
