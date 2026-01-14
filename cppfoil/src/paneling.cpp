// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/paneling.h"

#include <cmath>
#include <algorithm>

namespace cppfoil {

SpacingResult SplineCurvature(const MatrixXd& points, int N,
                               double Ufac, double TEfac,
                               const VectorXd& stgt) {
    // Get x range for TE detection
    double xmin = points.row(0).minCoeff();
    double xmax = points.row(0).maxCoeff();

    // Spline the input points
    Spline2D pp = Spline2DCreate(points);

    // Create fine sampling for curvature computation
    const int nfine = 501;
    double smax = pp.X.s(pp.X.s.size() - 1);
    VectorXd s_fine = VectorXd::LinSpaced(nfine, 0.0, smax);
    MatrixXd xy_fine = SplineEval(pp, s_fine);

    // Re-spline fine points
    Spline2D pp_fine = Spline2DCreate(xy_fine);

    SpacingResult result;

    if (stgt.size() == 0) {
        // Automatic curvature-based spacing
        VectorXd s = pp_fine.X.s;
        VectorXd sk = VectorXd::Zero(nfine);

        VectorXd xq, wq;
        QuadSeg(xq, wq);

        for (int i = 0; i < nfine - 1; ++i) {
            double ds = s(i + 1) - s(i);

            // Integrate curvature over segment
            double skint = 0.0;
            for (int q = 0; q < 5; ++q) {
                double sq = s(i) + xq(q) * ds;
                double xss = pp_fine.X.EvalDeriv2(sq);
                double yss = pp_fine.Y.EvalDeriv2(sq);
                skint += wq(q) * std::sqrt(xss * xss + yss * yss);
            }
            skint = 0.01 * Ufac + 0.5 * skint * ds;

            // Force TE resolution
            double xmid = 0.5 * (xy_fine(0, i) + xy_fine(0, i + 1));
            double xx = (xmid - xmin) / (xmax - xmin);  // Close to 1 means at TE
            skint += TEfac * 0.5 * std::exp(-100.0 * (1.0 - xx));

            sk(i + 1) = sk(i) + skint;
        }

        // Offset by fraction of average to avoid problems with zero curvature
        double sk_avg = sk.sum() / nfine;
        sk = sk.array() + 2.0 * sk_avg;

        // Interpolate to get arc-length values at uniformly spaced sk
        // Use cubic spline interpolation (matching mfoil's interp1d cubic)
        VectorXd sk_target = VectorXd::LinSpaced(N, sk(0), sk(nfine - 1));

        // Create cubic spline: sk -> s (using our 1D spline)
        CubicSpline1D sk_to_s = Spline1DCreate(sk, s);

        // Evaluate at target sk values
        VectorXd s_new(N);
        for (int i = 0; i < N; ++i) {
            s_new(i) = sk_to_s.Eval(sk_target(i));
        }

        result.S = s_new;
    } else {
        result.S = stgt;
    }

    // Evaluate points and tangents at final arc-lengths
    result.X = SplineEval(pp_fine, result.S);
    result.XS = SplineTangent(pp_fine, result.S);

    return result;
}

void MakePanels(const Geom& geom, int npanel, Panel& panel,
                const VectorXd& stgt) {
    int N = npanel + 1;  // Number of nodes

    // Use curvature-based spacing (must match mfoil exactly)
    double Ufac = 2.0;   // Uniformity factor (higher = more uniform)
    double TEfac = 0.1;  // TE resolution factor (higher = more TE resolution)

    SpacingResult spacing = SplineCurvature(geom.xpoint, N, Ufac, TEfac, stgt);

    panel.N = N;
    panel.x = spacing.X;
    panel.s = spacing.S;
    panel.t = spacing.XS;

    // Normalize tangents
    for (int i = 0; i < N; ++i) {
        double tnorm = panel.t.col(i).norm();
        if (tnorm > 1e-10) {
            panel.t.col(i) /= tnorm;
        }
    }
}

}  // namespace cppfoil
