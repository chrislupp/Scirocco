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
#ifndef SCIROCCO_CORE_SPLINE_H
#define SCIROCCO_CORE_SPLINE_H

#include "scirocco/config.h"

namespace scirocco {

// Cubic spline for one coordinate (x or y)
// Stores coefficients c[4, N-1] for N data points
// For segment i: f(t) = c(0,i)*t^3 + c(1,i)*t^2 + c(2,i)*t + c(3,i)
// where t is local parameter from 0 to ds = s[i+1] - s[i]
struct CubicSpline1D {
    VectorXd s;  // Parameter values at data points [N]
    MatrixXd c;  // Coefficients [4 x (N-1)]

    // Evaluate spline at parameter value
    [[nodiscard]] double Eval(double sval) const;

    // Evaluate derivative at parameter value
    [[nodiscard]] double EvalDeriv(double sval) const;

    // Evaluate second derivative at parameter value
    [[nodiscard]] double EvalDeriv2(double sval) const;
};

// 2D spline structure containing x and y splines
struct Spline2D {
    CubicSpline1D X;  // Spline for x coordinate
    CubicSpline1D Y;  // Spline for y coordinate
};

// Create a 2D cubic spline from points with arc-length parameterization
// INPUT:
//   points : 2D point coordinates [2 x N]
// OUTPUT:
//   Spline2D with arc-length corrected parameterization
Spline2D Spline2DCreate(const MatrixXd& points);

// Evaluate 2D spline at given parameter values
// INPUT:
//   spline : 2D spline structure
//   s_vals : parameter values at which to evaluate
// OUTPUT:
//   Points on spline [2 x len(s_vals)]
MatrixXd SplineEval(const Spline2D& spline, const VectorXd& s_vals);

// Evaluate 2D spline tangent (not normalized) at given parameter values
// INPUT:
//   spline : 2D spline structure
//   s_vals : parameter values at which to evaluate tangent
// OUTPUT:
//   dX/ds and dY/ds values at each point [2 x len(s_vals)]
MatrixXd SplineTangent(const Spline2D& spline, const VectorXd& s_vals);

// Create a 1D cubic spline from x and y data
// INPUT:
//   x : independent variable values [N]
//   y : dependent variable values [N]
// OUTPUT:
//   CubicSpline1D with x as parameter
CubicSpline1D Spline1DCreate(const VectorXd& x, const VectorXd& y);

// Quadrature points and weights for [0,1] segment
// Used for arc-length integration
// OUTPUT:
//   xq : quadrature point coordinates
//   wq : quadrature weights
void QuadSeg(VectorXd& xq, VectorXd& wq);

}  // namespace scirocco

#endif  // SCIROCCO_CORE_SPLINE_H
