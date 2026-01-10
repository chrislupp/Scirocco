// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Spline interpolation utilities.
// Ported from XFOIL spline.f (Mark Drela, 2000).

#ifndef XFOIL_SPLINE_H_
#define XFOIL_SPLINE_H_

#include "xfoil/config.h"

namespace xfoil {

/// @brief Solves tridiagonal system A*x = d in place.
/// @param a Diagonal (modified)
/// @param b Sub-diagonal
/// @param c Super-diagonal (modified)
/// @param d Right-hand side, replaced by solution
void trisol(VectorXd& a, const VectorXd& b, VectorXd& c, VectorXd& d);

/// @brief Calculates cubic spline coefficients with natural end conditions.
/// @param x Dependent variable array
/// @param s Independent variable array
/// @param xs Output: dx/ds derivative array
void spline(const VectorXd& x, const VectorXd& s, VectorXd& xs);

/// @brief Calculates cubic spline coefficients with specified end conditions.
/// @param x Dependent variable array
/// @param s Independent variable array
/// @param xs Output: dx/ds derivative array
/// @param xs1 First point derivative (999.0 = natural, -999.0 = zero 3rd deriv)
/// @param xs2 Last point derivative (999.0 = natural, -999.0 = zero 3rd deriv)
void splind(const VectorXd& x, const VectorXd& s, VectorXd& xs,
            double xs1, double xs2);

/// @brief Simple averaging spline (non-oscillatory).
/// @param x Dependent variable array
/// @param s Independent variable array
/// @param xs Output: dx/ds derivative array
void splina(const VectorXd& x, const VectorXd& s, VectorXd& xs);

/// @brief Binary search to find interval containing ss.
/// @param ss Value to locate
/// @param s Monotonic array
/// @return Index i such that s(i-1) <= ss < s(i)
int find_interval(double ss, const VectorXd& s);

/// @brief Evaluates spline at point ss.
/// @param ss Evaluation point
/// @param x Dependent variable array
/// @param xs Spline derivative array
/// @param s Independent variable array
/// @return Interpolated value x(ss)
double seval(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s);

/// @brief Evaluates first derivative dx/ds at point ss.
double deval(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s);

/// @brief Evaluates second derivative d2x/ds2 at point ss.
double d2val(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s);

/// @brief Calculates curvature of splined 2D curve at s = ss.
double curv(double ss, const VectorXd& x, const VectorXd& xs,
            const VectorXd& y, const VectorXd& ys, const VectorXd& s);

/// @brief Calculates inverse spline s(x) using Newton iteration.
/// @param si Initial guess, returns calculated s(xi)
/// @param xi Target x value
/// @return true if converged
bool sinvrt(double& si, double xi, const VectorXd& x,
            const VectorXd& xs, const VectorXd& s);

/// @brief Calculates arc length array s for 2D points (x, y).
void scalc(const VectorXd& x, const VectorXd& y, VectorXd& s);

/// @brief Splines x(s) with derivative discontinuities at segment joints.
void segspl(const VectorXd& x, const VectorXd& s, VectorXd& xs);

}  // namespace xfoil

#endif  // XFOIL_SPLINE_H_
