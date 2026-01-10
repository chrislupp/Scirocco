// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Spline interpolation utilities.

#include "xfoil/spline.h"

#include <cmath>

namespace xfoil {

void trisol(VectorXd& a, const VectorXd& b, VectorXd& c, VectorXd& d) {
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

void spline(const VectorXd& x, const VectorXd& s, VectorXd& xs) {
    int n = x.size();
    xs.resize(n);

    VectorXd a(n), b(n), c(n);

    // Interior points
    for (int i = 1; i < n - 1; ++i) {
        double dsm = s(i) - s(i - 1);
        double dsp = s(i + 1) - s(i);
        b(i) = dsp;
        a(i) = 2.0 * (dsm + dsp);
        c(i) = dsm;
        xs(i) = 3.0 * ((x(i + 1) - x(i)) * dsm / dsp +
                       (x(i) - x(i - 1)) * dsp / dsm);
    }

    // Zero second derivative end conditions
    a(0) = 2.0;
    c(0) = 1.0;
    xs(0) = 3.0 * (x(1) - x(0)) / (s(1) - s(0));

    b(n - 1) = 1.0;
    a(n - 1) = 2.0;
    xs(n - 1) = 3.0 * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));

    trisol(a, b, c, xs);
}

void splind(const VectorXd& x, const VectorXd& s, VectorXd& xs,
            double xs1, double xs2) {
    int n = x.size();
    xs.resize(n);

    VectorXd a(n), b(n), c(n);

    // Interior points
    for (int i = 1; i < n - 1; ++i) {
        double dsm = s(i) - s(i - 1);
        double dsp = s(i + 1) - s(i);
        b(i) = dsp;
        a(i) = 2.0 * (dsm + dsp);
        c(i) = dsm;
        xs(i) = 3.0 * ((x(i + 1) - x(i)) * dsm / dsp +
                       (x(i) - x(i - 1)) * dsp / dsm);
    }

    // First point conditions
    if (xs1 == 999.0) {
        a(0) = 2.0;
        c(0) = 1.0;
        xs(0) = 3.0 * (x(1) - x(0)) / (s(1) - s(0));
    } else if (xs1 == -999.0) {
        a(0) = 1.0;
        c(0) = 1.0;
        xs(0) = 2.0 * (x(1) - x(0)) / (s(1) - s(0));
    } else {
        a(0) = 1.0;
        c(0) = 0.0;
        xs(0) = xs1;
    }

    // Last point conditions
    if (xs2 == 999.0) {
        b(n - 1) = 1.0;
        a(n - 1) = 2.0;
        xs(n - 1) = 3.0 * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    } else if (xs2 == -999.0) {
        b(n - 1) = 1.0;
        a(n - 1) = 1.0;
        xs(n - 1) = 2.0 * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    } else {
        a(n - 1) = 1.0;
        b(n - 1) = 0.0;
        xs(n - 1) = xs2;
    }

    // Special case: n=2 with both zero 3rd derivative
    if (n == 2 && xs1 == -999.0 && xs2 == -999.0) {
        b(n - 1) = 1.0;
        a(n - 1) = 2.0;
        xs(n - 1) = 3.0 * (x(n - 1) - x(n - 2)) / (s(n - 1) - s(n - 2));
    }

    trisol(a, b, c, xs);
}

void splina(const VectorXd& x, const VectorXd& s, VectorXd& xs) {
    int n = x.size();
    xs.resize(n);

    bool lend = true;
    double xs1 = 0.0;
    double xs2 = 0.0;

    for (int i = 0; i < n - 1; ++i) {
        double ds = s(i + 1) - s(i);
        if (ds == 0.0) {
            xs(i) = xs1;
            lend = true;
        } else {
            double dx = x(i + 1) - x(i);
            xs2 = dx / ds;
            if (lend) {
                xs(i) = xs2;
                lend = false;
            } else {
                xs(i) = 0.5 * (xs1 + xs2);
            }
        }
        xs1 = xs2;
    }
    xs(n - 1) = xs1;
}

int find_interval(double ss, const VectorXd& s) {
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

double seval(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s) {
    int i = find_interval(ss, s);

    double ds = s(i) - s(i - 1);
    double t = (ss - s(i - 1)) / ds;
    double cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    double cx2 = ds * xs(i) - x(i) + x(i - 1);

    return t * x(i) + (1.0 - t) * x(i - 1) +
           (t - t * t) * ((1.0 - t) * cx1 - t * cx2);
}

double deval(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s) {
    int i = find_interval(ss, s);

    double ds = s(i) - s(i - 1);
    double t = (ss - s(i - 1)) / ds;
    double cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    double cx2 = ds * xs(i) - x(i) + x(i - 1);

    double result = x(i) - x(i - 1) +
                    (1.0 - 4.0 * t + 3.0 * t * t) * cx1 +
                    t * (3.0 * t - 2.0) * cx2;
    return result / ds;
}

double d2val(double ss, const VectorXd& x, const VectorXd& xs,
             const VectorXd& s) {
    int i = find_interval(ss, s);

    double ds = s(i) - s(i - 1);
    double t = (ss - s(i - 1)) / ds;
    double cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    double cx2 = ds * xs(i) - x(i) + x(i - 1);

    double result = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;
    return result / (ds * ds);
}

double curv(double ss, const VectorXd& x, const VectorXd& xs,
            const VectorXd& y, const VectorXd& ys, const VectorXd& s) {
    int i = find_interval(ss, s);

    double ds = s(i) - s(i - 1);
    double t = (ss - s(i - 1)) / ds;

    double cx1 = ds * xs(i - 1) - x(i) + x(i - 1);
    double cx2 = ds * xs(i) - x(i) + x(i - 1);
    double xd = x(i) - x(i - 1) +
                (1.0 - 4.0 * t + 3.0 * t * t) * cx1 +
                t * (3.0 * t - 2.0) * cx2;
    double xdd = (6.0 * t - 4.0) * cx1 + (6.0 * t - 2.0) * cx2;

    double cy1 = ds * ys(i - 1) - y(i) + y(i - 1);
    double cy2 = ds * ys(i) - y(i) + y(i - 1);
    double yd = y(i) - y(i - 1) +
                (1.0 - 4.0 * t + 3.0 * t * t) * cy1 +
                t * (3.0 * t - 2.0) * cy2;
    double ydd = (6.0 * t - 4.0) * cy1 + (6.0 * t - 2.0) * cy2;

    double sd = std::sqrt(xd * xd + yd * yd);
    sd = std::max(sd, 0.001 * ds);

    return (xd * ydd - yd * xdd) / (sd * sd * sd);
}

bool sinvrt(double& si, double xi, const VectorXd& x,
            const VectorXd& xs, const VectorXd& s) {
    int n = s.size();
    double sisav = si;

    for (int iter = 0; iter < 10; ++iter) {
        double res = seval(si, x, xs, s) - xi;
        double resp = deval(si, x, xs, s);
        double ds = -res / resp;
        si = si + ds;

        if (std::abs(ds / (s(n - 1) - s(0))) < 1.0e-5) {
            return true;
        }
    }

    si = sisav;
    return false;
}

void scalc(const VectorXd& x, const VectorXd& y, VectorXd& s) {
    int n = x.size();
    s.resize(n);
    s(0) = 0.0;

    for (int i = 1; i < n; ++i) {
        double dx = x(i) - x(i - 1);
        double dy = y(i) - y(i - 1);
        s(i) = s(i - 1) + std::sqrt(dx * dx + dy * dy);
    }
}

void segspl(const VectorXd& x, const VectorXd& s, VectorXd& xs) {
    int n = x.size();
    xs.resize(n);

    int iseg0 = 0;
    for (int iseg = 1; iseg < n - 2; ++iseg) {
        if (s(iseg) == s(iseg + 1)) {
            // Segment boundary found
            int nseg = iseg - iseg0 + 1;
            VectorXd x_seg = x.segment(iseg0, nseg);
            VectorXd s_seg = s.segment(iseg0, nseg);
            VectorXd xs_seg;
            splind(x_seg, s_seg, xs_seg, -999.0, -999.0);
            xs.segment(iseg0, nseg) = xs_seg;
            iseg0 = iseg + 1;
        }
    }

    // Last segment
    int nseg = n - iseg0;
    VectorXd x_seg = x.segment(iseg0, nseg);
    VectorXd s_seg = s.segment(iseg0, nseg);
    VectorXd xs_seg;
    splind(x_seg, s_seg, xs_seg, -999.0, -999.0);
    xs.segment(iseg0, nseg) = xs_seg;
}

}  // namespace xfoil
