// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Spline interpolation utilities.
// Ported from XFOIL spline.f (Mark Drela, 2000).

#ifndef XFOIL_SPLINE_H_
#define XFOIL_SPLINE_H_

#include <Eigen/Dense>
#include <cmath>

#include "xfoil/config.h"

namespace xfoil {

/// @brief Solves tridiagonal system A*x = d, where A has diagonal a,
///        sub-diagonal b, super-diagonal c.
/// @param a Diagonal (destroyed on output)
/// @param b Sub-diagonal
/// @param c Super-diagonal (destroyed on output)
/// @param d Right-hand side (replaced by solution on output)
template <typename Scalar>
void trisol(VectorX<Scalar>& a, const VectorX<Scalar>& b,
            VectorX<Scalar>& c, VectorX<Scalar>& d) {
    int n = a.size();

    // Forward elimination
    for (int k = 1; k < n; ++k) {
        c(k - 1) = c(k - 1) / a(k - 1);
        d(k - 1) = d(k - 1) / a(k - 1);
        a(k) = a(k) - b(k) * c(k - 1);
        d(k) = d(k) - b(k) * d(k - 1);
    }

    // Back substitution
    d(n - 1) = d(n - 1) / a(n - 1);
    for (int k = n - 2; k >= 0; --k) {
        d(k) = d(k) - c(k) * d(k + 1);
    }
}

/// @brief Calculates cubic spline coefficients for x(s).
///        Uses natural (zero second derivative) end conditions.
/// @param x Dependent variable array (input)
/// @param s Independent variable array (input)
/// @param xs dx/ds derivative array (output)
template <typename Scalar>
void spline(const VectorX<Scalar>& x, const VectorX<Scalar>& s,
            VectorX<Scalar>& xs) {
    int n = x.size();
    xs.resize(n);

    VectorX<Scalar> a(n), b(n), c(n);

    // Interior points
    for (int i = 1; i < n - 1; ++i) {
        Scalar dsm = s(i) - s(i - 1);
        Scalar dsp = s(i + 1) - s(i);
        b(i) = dsp;
        a(i) = Scalar(2.0) * (dsm + dsp);
        c(i) = dsm;
        xs(i) = Scalar(3.0) * ((x(i + 1) - x(i)) * dsm / dsp +
                               (x(i) - x(i - 1)) * dsp / dsm);
    }

    // Zero second derivative end conditions
    a(0) = Scalar(2.0);
    c(0) = Scalar(1.0);
    xs(0) = Scalar(3.0) * (x(1) - x(0)) / (s(1) - s(0));

    b(n - 1) = Scalar(1.0);
    a(n - 1) = Scalar(2.0);
    xs(n - 1) = Scalar(3.0) * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));

    // Solve tridiagonal system
    trisol(a, b, c, xs);
}

/// @brief Calculates cubic spline coefficients with specified end conditions.
/// @param x Dependent variable array (input)
/// @param s Independent variable array (input)
/// @param xs dx/ds derivative array (output)
/// @param xs1 First point derivative (999.0 = natural, -999.0 = zero 3rd deriv)
/// @param xs2 Last point derivative (999.0 = natural, -999.0 = zero 3rd deriv)
template <typename Scalar>
void splind(const VectorX<Scalar>& x, const VectorX<Scalar>& s,
            VectorX<Scalar>& xs, Scalar xs1, Scalar xs2) {
    int n = x.size();
    xs.resize(n);

    VectorX<Scalar> a(n), b(n), c(n);

    // Interior points
    for (int i = 1; i < n - 1; ++i) {
        Scalar dsm = s(i) - s(i - 1);
        Scalar dsp = s(i + 1) - s(i);
        b(i) = dsp;
        a(i) = Scalar(2.0) * (dsm + dsp);
        c(i) = dsm;
        xs(i) = Scalar(3.0) * ((x(i + 1) - x(i)) * dsm / dsp +
                               (x(i) - x(i - 1)) * dsp / dsm);
    }

    // First point conditions
    if (xs1 == Scalar(999.0)) {
        // Zero second derivative
        a(0) = Scalar(2.0);
        c(0) = Scalar(1.0);
        xs(0) = Scalar(3.0) * (x(1) - x(0)) / (s(1) - s(0));
    } else if (xs1 == Scalar(-999.0)) {
        // Zero third derivative
        a(0) = Scalar(1.0);
        c(0) = Scalar(1.0);
        xs(0) = Scalar(2.0) * (x(1) - x(0)) / (s(1) - s(0));
    } else {
        // Specified first derivative
        a(0) = Scalar(1.0);
        c(0) = Scalar(0.0);
        xs(0) = xs1;
    }

    // Last point conditions
    if (xs2 == Scalar(999.0)) {
        b(n - 1) = Scalar(1.0);
        a(n - 1) = Scalar(2.0);
        xs(n - 1) = Scalar(3.0) * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    } else if (xs2 == Scalar(-999.0)) {
        b(n - 1) = Scalar(1.0);
        a(n - 1) = Scalar(1.0);
        xs(n - 1) = Scalar(2.0) * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    } else {
        a(n - 1) = Scalar(1.0);
        b(n - 1) = Scalar(0.0);
        xs(n - 1) = xs2;
    }

    // Special case: n=2 with both zero 3rd derivative
    if (n == 2 && xs1 == Scalar(-999.0) && xs2 == Scalar(-999.0)) {
        b(n - 1) = Scalar(1.0);
        a(n - 1) = Scalar(2.0);
        xs(n - 1) = Scalar(3.0) * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    }

    // Solve tridiagonal system
    trisol(a, b, c, xs);
}

/// @brief Simple averaging spline (non-oscillatory).
/// @param x Dependent variable array (input)
/// @param s Independent variable array (input)
/// @param xs dx/ds derivative array (output)
template <typename Scalar>
void splina(const VectorX<Scalar>& x, const VectorX<Scalar>& s,
            VectorX<Scalar>& xs) {
    int n = x.size();
    xs.resize(n);

    bool lend = true;
    Scalar xs1 = Scalar(0.0);
    Scalar xs2 = Scalar(0.0);

    for (int i = 0; i < n - 1; ++i) {
        Scalar ds = s(i + 1) - s(i);
        if (ds == Scalar(0.0)) {
            xs(i) = xs1;
            lend = true;
        } else {
            Scalar dx = x(i + 1) - x(i);
            xs2 = dx / ds;
            if (lend) {
                xs(i) = xs2;
                lend = false;
            } else {
                xs(i) = Scalar(0.5) * (xs1 + xs2);
            }
        }
        xs1 = xs2;
    }
    xs(n - 1) = xs1;
}

/// @brief Binary search to find interval containing ss.
/// @param ss Value to locate
/// @param s Monotonic array
/// @return Index i such that s(i-1) <= ss < s(i)
template <typename Scalar>
int find_interval(Scalar ss, const VectorX<Scalar>& s) {
    int n = s.size();
    int ilow = 0;
    int i = n - 1;

    while (i - ilow > 1) {
        int imid = (i + ilow) / 2;
        if (ss < s(imid)) {
            i = imid;
        } else {
            ilow = imid;
        }
    }
    return i;
}

/// @brief Evaluates spline at point ss.
/// @param ss Evaluation point
/// @param x Dependent variable array
/// @param xs Spline derivative array (from spline())
/// @param s Independent variable array
/// @return Interpolated value x(ss)
template <typename Scalar>
Scalar seval(Scalar ss, const VectorX<Scalar>& x, const VectorX<Scalar>& xs,
             const VectorX<Scalar>& s) {
    int i = find_interval(ss, s);

    Scalar ds = s(i) - s(i - 1);
    Scalar t = (ss - s(i - 1)) / ds;
    Scalar cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    Scalar cx2 = ds * xs(i) - x(i) + x(i - 1);

    return t * x(i) + (Scalar(1.0) - t) * x(i - 1) +
           (t - t * t) * ((Scalar(1.0) - t) * cx1 - t * cx2);
}

/// @brief Evaluates first derivative dx/ds at point ss.
/// @param ss Evaluation point
/// @param x Dependent variable array
/// @param xs Spline derivative array (from spline())
/// @param s Independent variable array
/// @return First derivative dx/ds(ss)
template <typename Scalar>
Scalar deval(Scalar ss, const VectorX<Scalar>& x, const VectorX<Scalar>& xs,
             const VectorX<Scalar>& s) {
    int i = find_interval(ss, s);

    Scalar ds = s(i) - s(i - 1);
    Scalar t = (ss - s(i - 1)) / ds;
    Scalar cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    Scalar cx2 = ds * xs(i) - x(i) + x(i - 1);

    Scalar result = x(i) - x(i - 1) +
                    (Scalar(1.0) - Scalar(4.0) * t + Scalar(3.0) * t * t) * cx1 +
                    t * (Scalar(3.0) * t - Scalar(2.0)) * cx2;
    return result / ds;
}

/// @brief Evaluates second derivative d2x/ds2 at point ss.
/// @param ss Evaluation point
/// @param x Dependent variable array
/// @param xs Spline derivative array (from spline())
/// @param s Independent variable array
/// @return Second derivative d2x/ds2(ss)
template <typename Scalar>
Scalar d2val(Scalar ss, const VectorX<Scalar>& x, const VectorX<Scalar>& xs,
             const VectorX<Scalar>& s) {
    int i = find_interval(ss, s);

    Scalar ds = s(i) - s(i - 1);
    Scalar t = (ss - s(i - 1)) / ds;
    Scalar cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    Scalar cx2 = ds * xs(i) - x(i) + x(i - 1);

    Scalar result = (Scalar(6.0) * t - Scalar(4.0)) * cx1 +
                    (Scalar(6.0) * t - Scalar(2.0)) * cx2;
    return result / (ds * ds);
}

/// @brief Calculates curvature of splined 2D curve at s = ss.
/// @param ss Evaluation point
/// @param x X coordinate array
/// @param xs X spline derivative array
/// @param y Y coordinate array
/// @param ys Y spline derivative array
/// @param s Arc length array
/// @return Curvature at ss
template <typename Scalar>
Scalar curv(Scalar ss, const VectorX<Scalar>& x, const VectorX<Scalar>& xs,
            const VectorX<Scalar>& y, const VectorX<Scalar>& ys,
            const VectorX<Scalar>& s) {
    int i = find_interval(ss, s);

    Scalar ds = s(i) - s(i - 1);
    Scalar t = (ss - s(i - 1)) / ds;

    Scalar cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    Scalar cx2 = ds * xs(i) - x(i) + x(i - 1);
    Scalar xd = x(i) - x(i - 1) +
                (Scalar(1.0) - Scalar(4.0) * t + Scalar(3.0) * t * t) * cx1 +
                t * (Scalar(3.0) * t - Scalar(2.0)) * cx2;
    Scalar xdd = (Scalar(6.0) * t - Scalar(4.0)) * cx1 +
                 (Scalar(6.0) * t - Scalar(2.0)) * cx2;

    Scalar cy1 = ds * ys(i - 1) - y(i) + y(i - 1);
    Scalar cy2 = ds * ys(i) - y(i) + y(i - 1);
    Scalar yd = y(i) - y(i - 1) +
                (Scalar(1.0) - Scalar(4.0) * t + Scalar(3.0) * t * t) * cy1 +
                t * (Scalar(3.0) * t - Scalar(2.0)) * cy2;
    Scalar ydd = (Scalar(6.0) * t - Scalar(4.0)) * cy1 +
                 (Scalar(6.0) * t - Scalar(2.0)) * cy2;

    using std::max;
    using std::sqrt;
    Scalar sd = sqrt(xd * xd + yd * yd);
    sd = max(sd, Scalar(0.001) * ds);

    return (xd * ydd - yd * xdd) / (sd * sd * sd);
}

/// @brief Calculates the inverse spline function s(x) using Newton iteration.
/// @param si Initial guess for s, returns calculated s(xi)
/// @param xi Target x value
/// @param x Dependent variable array
/// @param xs Spline derivative array
/// @param s Independent variable array
/// @return true if converged, false otherwise
template <typename Scalar>
bool sinvrt(Scalar& si, Scalar xi, const VectorX<Scalar>& x,
            const VectorX<Scalar>& xs, const VectorX<Scalar>& s) {
    int n = s.size();
    Scalar sisav = si;

    for (int iter = 0; iter < 10; ++iter) {
        Scalar res = seval(si, x, xs, s) - xi;
        Scalar resp = deval(si, x, xs, s);
        Scalar ds = -res / resp;
        si = si + ds;

        using std::abs;
        if (abs(ds / (s(n - 1) - s(0))) < Scalar(1.0e-5)) {
            return true;
        }
    }

    // Failed to converge
    si = sisav;
    return false;
}

/// @brief Calculates arc length array s for 2D points (x, y).
/// @param x X coordinate array (input)
/// @param y Y coordinate array (input)
/// @param s Arc length array (output)
template <typename Scalar>
void scalc(const VectorX<Scalar>& x, const VectorX<Scalar>& y,
           VectorX<Scalar>& s) {
    int n = x.size();
    s.resize(n);
    s(0) = Scalar(0.0);

    for (int i = 1; i < n; ++i) {
        using std::sqrt;
        Scalar dx = x(i) - x(i - 1);
        Scalar dy = y(i) - y(i - 1);
        s(i) = s(i - 1) + sqrt(dx * dx + dy * dy);
    }
}

/// @brief Splines x(s) with derivative discontinuities at segment joints.
///        Segment joints are defined by identical successive s values.
/// @param x Dependent variable array (input)
/// @param s Independent variable array (input)
/// @param xs dx/ds derivative array (output)
template <typename Scalar>
void segspl(const VectorX<Scalar>& x, const VectorX<Scalar>& s,
            VectorX<Scalar>& xs) {
    int n = x.size();
    xs.resize(n);

    int iseg0 = 0;
    for (int iseg = 1; iseg < n - 2; ++iseg) {
        if (s(iseg) == s(iseg + 1)) {
            // Segment boundary found
            int nseg = iseg - iseg0 + 1;
            VectorX<Scalar> x_seg = x.segment(iseg0, nseg);
            VectorX<Scalar> s_seg = s.segment(iseg0, nseg);
            VectorX<Scalar> xs_seg;
            splind(x_seg, s_seg, xs_seg, Scalar(-999.0), Scalar(-999.0));
            xs.segment(iseg0, nseg) = xs_seg;
            iseg0 = iseg + 1;
        }
    }

    // Last segment
    int nseg = n - iseg0;
    VectorX<Scalar> x_seg = x.segment(iseg0, nseg);
    VectorX<Scalar> s_seg = s.segment(iseg0, nseg);
    VectorX<Scalar> xs_seg;
    splind(x_seg, s_seg, xs_seg, Scalar(-999.0), Scalar(-999.0));
    xs.segment(iseg0, nseg) = xs_seg;
}

}  // namespace xfoil

#endif  // XFOIL_SPLINE_H_
