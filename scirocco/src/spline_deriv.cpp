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
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include <array>
#include <stdexcept>

#include "scirocco/core/spline.h"

namespace scirocco {
namespace {
// Each row contains derivatives with respect to all interleaved node
// coordinates.
struct SplineLin {
    MatrixXd c;
    std::array<MatrixXd, 4> dc;
};

SplineLin LinearizeSpline(const VectorXd& s, const MatrixXd& ds,
                          const VectorXd& y, const MatrixXd& dy) {
    const int n = s.size(), nd = ds.cols();
    VectorXd h = s.tail(n - 1) - s.head(n - 1);
    MatrixXd dh = ds.bottomRows(n - 1) - ds.topRows(n - 1);
    VectorXd secant = (y.tail(n - 1) - y.head(n - 1)).cwiseQuotient(h);
    MatrixXd dsecant = dy.bottomRows(n - 1) - dy.topRows(n - 1);
    for (int i = 0; i < n - 1; ++i)
        dsecant.row(i) = (dsecant.row(i) - secant(i) * dh.row(i)) / h(i);

    VectorXd z = VectorXd::Zero(n);
    MatrixXd dz = MatrixXd::Zero(n, nd);
    if (n == 3) {
        const double span = h.sum();
        z.setConstant(2 * (secant(1) - secant(0)) / span);
        for (int i = 0; i < n; ++i)
            dz.row(i) = (2 * (dsecant.row(1) - dsecant.row(0)) -
                         z(i) * dh.colwise().sum()) /
                        span;
    } else if (n > 3) {
        // Second-derivative system for a not-a-knot cubic spline.
        MatrixXd a = MatrixXd::Zero(n, n);
        VectorXd b = VectorXd::Zero(n);
        for (int i = 1; i < n - 1; ++i) {
            a(i, i - 1) = h(i - 1);
            a(i, i) = 2 * (h(i - 1) + h(i));
            a(i, i + 1) = h(i);
            b(i) = 6 * (secant(i) - secant(i - 1));
        }
        a(0, 0) = h(1);
        a(0, 1) = -h(0) - h(1);
        a(0, 2) = h(0);
        a(n - 1, n - 3) = h(n - 2);
        a(n - 1, n - 2) = -h(n - 3) - h(n - 2);
        a(n - 1, n - 1) = h(n - 3);
        auto lu = a.partialPivLu();
        z = lu.solve(b);
        // A dz = db - (dA) z. Reuse the real factorization for every direction.
        MatrixXd rhs = MatrixXd::Zero(n, nd);
        for (int i = 1; i < n - 1; ++i)
            rhs.row(i) = 6 * (dsecant.row(i) - dsecant.row(i - 1)) -
                         (z(i - 1) + 2 * z(i)) * dh.row(i - 1) -
                         (2 * z(i) + z(i + 1)) * dh.row(i);
        rhs.row(0) = -(z(0) - z(1)) * dh.row(1) - (z(2) - z(1)) * dh.row(0);
        rhs.row(n - 1) = -(z(n - 3) - z(n - 2)) * dh.row(n - 2) -
                         (z(n - 1) - z(n - 2)) * dh.row(n - 3);
        dz = lu.solve(rhs);
    }
    SplineLin out;
    out.c.resize(4, n - 1);
    for (auto& d : out.dc) d.resize(n - 1, nd);
    for (int i = 0; i < n - 1; ++i) {
        out.c(0, i) = (z(i + 1) - z(i)) / (6 * h(i));
        out.dc[0].row(i) =
            ((dz.row(i + 1) - dz.row(i)) / 6 - out.c(0, i) * dh.row(i)) / h(i);
        out.c(1, i) = z(i) / 2;
        out.dc[1].row(i) = dz.row(i) / 2;
        out.c(2, i) = secant(i) - h(i) * (2 * z(i) + z(i + 1)) / 6;
        out.dc[2].row(i) =
            dsecant.row(i) - (dh.row(i) * (2 * z(i) + z(i + 1)) +
                              h(i) * (2 * dz.row(i) + dz.row(i + 1))) /
                                 6;
        out.c(3, i) = y(i);
        out.dc[3].row(i) = dy.row(i);
    }
    return out;
}
}  // namespace

MatrixXd SplineArcLengthDeriv(const MatrixXd& points) {
    const int n = points.cols(), nd = 2 * n;
    if (n < 2) throw std::invalid_argument("Need at least 2 spline points");
    VectorXd s = VectorXd::Zero(n);
    MatrixXd ds = MatrixXd::Zero(n, nd);
    std::array<MatrixXd, 2> dy{MatrixXd::Zero(n, nd), MatrixXd::Zero(n, nd)};
    for (int i = 0; i < n; ++i) {
        dy[0](i, 2 * i) = 1;
        dy[1](i, 2 * i + 1) = 1;
        if (i == 0) continue;
        Vector2d edge = points.col(i) - points.col(i - 1);
        const double length = edge.norm();
        if (!(length > 0)) throw std::invalid_argument("Repeated spline point");
        s(i) = s(i - 1) + length;
        ds.row(i) = ds.row(i - 1);
        ds.block<1, 2>(i, 2 * i) += edge.transpose() / length;
        ds.block<1, 2>(i, 2 * (i - 1)) -= edge.transpose() / length;
    }
    VectorXd q, w;
    QuadSeg(q, w);
    // Python mfoil.spline2d performs exactly ten arc-length reintegrations.
    for (int pass = 0; pass < 10; ++pass) {
        std::array<SplineLin, 2> p{
            LinearizeSpline(s, ds, points.row(0), dy[0]),
            LinearizeSpline(s, ds, points.row(1), dy[1])};
        VectorXd next = VectorXd::Zero(n);
        MatrixXd dnext = MatrixXd::Zero(n, nd);
        for (int i = 0; i < n - 1; ++i) {
            const double h = s(i + 1) - s(i);
            const Eigen::RowVectorXd dh = ds.row(i + 1) - ds.row(i);
            double integral = 0;
            Eigen::RowVectorXd dintegral = Eigen::RowVectorXd::Zero(nd);
            for (int j = 0; j < q.size(); ++j) {
                const double t = q(j) * h;
                Vector2d v;
                MatrixXd dv(2, nd);
                for (int d = 0; d < 2; ++d) {
                    const auto& c = p[d].c;
                    v(d) = (3 * c(0, i) * t + 2 * c(1, i)) * t + c(2, i);
                    dv.row(d) = 3 * t * t * p[d].dc[0].row(i) +
                                2 * t * p[d].dc[1].row(i) + p[d].dc[2].row(i) +
                                (6 * c(0, i) * t + 2 * c(1, i)) * q(j) * dh;
                }
                const double speed = v.norm();
                if (!(speed > 0))
                    throw std::runtime_error("Zero spline tangent");
                integral += w(j) * speed;
                dintegral += w(j) * (v.transpose() * dv) / speed;
            }
            next(i + 1) = next(i) + h * integral;
            dnext.row(i + 1) = dnext.row(i) + dh * integral + h * dintegral;
        }
        s = next;
        ds = dnext;
    }
    return ds;
}
}  // namespace scirocco
