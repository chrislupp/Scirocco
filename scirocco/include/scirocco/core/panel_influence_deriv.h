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
#ifndef SCIROCCO_CORE_PANEL_INFLUENCE_DERIV_H
#define SCIROCCO_CORE_PANEL_INFLUENCE_DERIV_H

#include "scirocco/config.h"
#include "scirocco/core/panel_influence.h"

namespace scirocco {

// Derivative information for PanelInfo quantities
// Derivatives are w.r.t. Xj = [x1_x, x1_z, x2_x, x2_z] and xi = [xi_x, xi_z]
struct PanelInfoDeriv {
    PanelInfo info;

    // Derivatives of local coordinates w.r.t. Xj [4] and xi [2]
    Vector4d x1_Xj, x2_Xj, y_Xj;
    Vector2d x1_xi, x2_xi, y_xi;

    // Derivatives of distances
    Vector4d r1_Xj, r2_Xj;
    Vector2d r1_xi, r2_xi;

    // Derivatives of angles
    Vector4d theta1_Xj, theta2_Xj;
    Vector2d theta1_xi, theta2_xi;

    // Derivatives of panel length
    Vector4d len_Xj;
};

// Compute PanelInfo with derivatives w.r.t. Xj and xi
PanelInfoDeriv GetPanelInfoDeriv(const Matrix2d& Xj, const Vector2d& xi);

// Linear vortex streamfunction with coordinate derivatives
// Returns a, b and their derivatives w.r.t. Xj [4] and xi [2]
struct LinVortexStreamDeriv {
    double a, b;
    Vector4d a_Xj, b_Xj;
    Vector2d a_xi, b_xi;
};
LinVortexStreamDeriv PanelLinVortexStreamDeriv(const Matrix2d& Xj,
                                               const Vector2d& xi);

// Constant source streamfunction with coordinate derivatives
struct ConstSourceStreamDeriv {
    double val;
    Vector4d val_Xj;
    Vector2d val_xi;
};
ConstSourceStreamDeriv PanelConstSourceStreamDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi);

// Linear vortex velocity with coordinate derivatives
struct LinVortexVelocityDeriv {
    Vector2d a, b;
    Eigen::Matrix<double, 2, 4> a_Xj, b_Xj;
    Matrix2d a_xi, b_xi;
};
LinVortexVelocityDeriv PanelLinVortexVelocityDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi);

// Constant source velocity with coordinate derivatives
struct ConstSourceVelocityDeriv {
    Vector2d val;
    Eigen::Matrix<double, 2, 4> val_Xj;
    Matrix2d val_xi;
};
ConstSourceVelocityDeriv PanelConstSourceVelocityDeriv(const Matrix2d& Xj,
                                                       const Vector2d& xi);

// Linear source streamfunction with coordinate derivatives
struct LinSourceStreamDeriv {
    double a, b;
    Vector4d a_Xj, b_Xj;
    Vector2d a_xi, b_xi;
};
LinSourceStreamDeriv PanelLinSourceStreamDeriv(const Matrix2d& Xj,
                                               const Vector2d& xi);

// Linear source velocity with coordinate derivatives
struct LinSourceVelocityDeriv {
    Vector2d a, b;
    Eigen::Matrix<double, 2, 4> a_Xj, b_Xj;
    Matrix2d a_xi, b_xi;
};
LinSourceVelocityDeriv PanelLinSourceVelocityDeriv(const Matrix2d& Xj,
                                                   const Vector2d& xi);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_PANEL_INFLUENCE_DERIV_H
