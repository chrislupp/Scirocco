// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#ifndef CPPFOIL_CORE_SPLINE_H
#define CPPFOIL_CORE_SPLINE_H

#include "cppfoil/config.h"

namespace cppfoil {

// Cubic spline for one coordinate (x or y)
// Stores coefficients c[4, N-1] for N data points
// For segment i: f(t) = c(0,i)*t^3 + c(1,i)*t^2 + c(2,i)*t + c(3,i)
// where t is local parameter from 0 to ds = s[i+1] - s[i]
struct CubicSpline1D {
    VectorXd s;       // Parameter values at data points [N]
    MatrixXd c;       // Coefficients [4 x (N-1)]

    // Evaluate spline at parameter value
    double Eval(double sval) const;

    // Evaluate derivative at parameter value
    double EvalDeriv(double sval) const;

    // Evaluate second derivative at parameter value
    double EvalDeriv2(double sval) const;
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

// Quadrature points and weights for [0,1] segment
// Used for arc-length integration
// OUTPUT:
//   xq : quadrature point coordinates
//   wq : quadrature weights
void QuadSeg(VectorXd& xq, VectorXd& wq);

}  // namespace cppfoil

#endif  // CPPFOIL_CORE_SPLINE_H
