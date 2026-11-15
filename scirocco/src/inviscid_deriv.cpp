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
#include "scirocco/core/inviscid_deriv.h"

#include <cmath>
#include <vector>

#include "scirocco/core/geometry.h"
#include "scirocco/core/panel_influence_deriv.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

namespace {

/* d(tcp)/d(foil.x) and d(tdp)/d(foil.x) for one TE-adjacent coordinate.
 * The TE source and vortex panel strengths scale with tcp = |t x p| and
 * tdp = t . p, where t is the bisector of the two TE panel tangents and p the
 * unit TE gap vector; both depend only on nodes 0, 1, N-2 and N-1. */
struct TEFactorDeriv {
    int col;    /* column in the [2N] foil.x layout */
    double tcp; /* d(tcp)/d(foil.x[col]) */
    double tdp; /* d(tdp)/d(foil.x[col]) */
};

std::vector<TEFactorDeriv> TEFactorDerivs(const MatrixXd& x) {
    const int N = static_cast<int>(x.cols());
    std::vector<TEFactorDeriv> out;
    const Vector2d s = x.col(N - 1) - x.col(0);
    const double sn = s.norm();
    if (sn <= 1e-10) return out;

    const Vector2d v1 = x.col(0) - x.col(1);
    const Vector2d v2 = x.col(N - 1) - x.col(N - 2);
    const double v1n = v1.norm(), v2n = v2.norm();
    if (v1n < 1e-10 || v2n < 1e-10) return out;
    const Vector2d t1 = v1 / v1n, t2 = v2 / v2n;
    const Vector2d traw = 0.5 * (t1 + t2);
    const double tnorm = traw.norm();
    if (tnorm < 1e-10) return out;
    const Vector2d tvec = traw / tnorm;
    const Vector2d pvec = s / sn;

    const Matrix2d dt1_dv1 = (Matrix2d::Identity() - t1 * t1.transpose()) / v1n;
    const Matrix2d dt2_dv2 = (Matrix2d::Identity() - t2 * t2.transpose()) / v2n;
    const Matrix2d dt_dtraw =
        (Matrix2d::Identity() - tvec * tvec.transpose()) / tnorm;
    const Matrix2d dp_ds =
        (Matrix2d::Identity() - pvec * pvec.transpose()) / sn;
    const double cross_val = tvec(0) * pvec(1) - tvec(1) * pvec(0);
    const double sign_c = (cross_val >= 0) ? 1.0 : -1.0;

    const int te_nodes[] = {0, 1, N - 2, N - 1};
    for (int kk : te_nodes) {
        if (kk < 0 || kk >= N) continue;
        for (int dd = 0; dd < 2; ++dd) {
            Vector2d dv1 = Vector2d::Zero();
            if (kk == 0)
                dv1(dd) = 1.0;
            else if (kk == 1)
                dv1(dd) = -1.0;
            Vector2d dv2 = Vector2d::Zero();
            if (kk == N - 1)
                dv2(dd) = 1.0;
            else if (kk == N - 2)
                dv2(dd) = -1.0;
            const Vector2d dtraw_d = 0.5 * (dt1_dv1 * dv1 + dt2_dv2 * dv2);
            const Vector2d dt_d = dt_dtraw * dtraw_d;

            Vector2d ds_d = Vector2d::Zero();
            if (kk == N - 1)
                ds_d(dd) = 1.0;
            else if (kk == 0)
                ds_d(dd) = -1.0;
            const Vector2d dp_d = dp_ds * ds_d;

            const double dcross = dt_d(0) * pvec(1) + tvec(0) * dp_d(1) -
                                  dt_d(1) * pvec(0) - tvec(1) * dp_d(0);
            out.push_back({2 * kk + dd, sign_c * dcross,
                           dt_d.dot(pvec) + tvec.dot(dp_d)});
        }
    }
    return out;
}

}  // namespace

MatrixXd AICMatVecJacobianX(const FoilState& state, const VectorXd& v) {
    /* Row i of d(A*v)/d(foil.x). Perturbing foil node k affects:
     *   (a) row k: the evaluation point xi moves, for every panel;
     *   (b) panels k-1 -> k and k -> k+1: a panel endpoint moves;
     *   (c) the TE source/vortex panel (nodes 0 and N-1);
     *   (d) the TE factors tcp/tdp (nodes 0, 1, N-2, N-1).
     * The Kutta row N has no coordinate dependence. With no TE gap, row N-1
     * is the gamma-extrapolation condition and has none either. */
    const int N = state.foil.N;
    const MatrixXd& x = state.foil.x;
    MatrixXd jac = MatrixXd::Zero(N + 1, 2 * N);

    const TEInfo te = GetTEInfo(x);
    const bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);
    Matrix2d XjTE;
    XjTE.col(0) = x.col(N - 1);
    XjTE.col(1) = x.col(0);
    const double sigma_factor = 0.5 * te.tcp * (v(N - 1) - v(0));
    const double vortex_factor = 0.5 * te.tdp * (v(0) - v(N - 1));

    for (int i = 0; i < N; ++i) {
        if (nogap && i == N - 1) continue;
        const Vector2d xi = x.col(i);

        /* (1) Linear vortex panels: R(i) += a*v(j) + b*v(j+1). */
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = x.col(j);
            Xj.col(1) = x.col(j + 1);
            const auto deriv = PanelLinVortexStreamDeriv(Xj, xi);
            const double vj = v(j), vjp = v(j + 1);
            jac(i, 2 * i) += deriv.a_xi(0) * vj + deriv.b_xi(0) * vjp;
            jac(i, 2 * i + 1) += deriv.a_xi(1) * vj + deriv.b_xi(1) * vjp;
            const Vector4d dXj = deriv.a_Xj * vj + deriv.b_Xj * vjp;
            jac(i, 2 * j) += dXj(0);
            jac(i, 2 * j + 1) += dXj(1);
            jac(i, 2 * (j + 1)) += dXj(2);
            jac(i, 2 * (j + 1) + 1) += dXj(3);
        }

        /* (2) TE source panel, endpoints [N-1, 0]:
         * R(i) += a_src * 0.5*tcp * (v(N-1) - v(0)). */
        {
            const auto deriv = PanelConstSourceStreamDeriv(XjTE, xi);
            jac(i, 2 * i) += deriv.val_xi(0) * sigma_factor;
            jac(i, 2 * i + 1) += deriv.val_xi(1) * sigma_factor;
            jac(i, 2 * (N - 1)) += deriv.val_Xj(0) * sigma_factor;
            jac(i, 2 * (N - 1) + 1) += deriv.val_Xj(1) * sigma_factor;
            jac(i, 0) += deriv.val_Xj(2) * sigma_factor;
            jac(i, 1) += deriv.val_Xj(3) * sigma_factor;
        }

        /* (3) TE vortex panel, a uniform vortex:
         * R(i) += (a+b) * 0.5*tdp * (v(0) - v(N-1)). */
        {
            const auto deriv = PanelLinVortexStreamDeriv(XjTE, xi);
            jac(i, 2 * i) += (deriv.a_xi(0) + deriv.b_xi(0)) * vortex_factor;
            jac(i, 2 * i + 1) +=
                (deriv.a_xi(1) + deriv.b_xi(1)) * vortex_factor;
            const Vector4d dab_dXj = deriv.a_Xj + deriv.b_Xj;
            jac(i, 2 * (N - 1)) += dab_dXj(0) * vortex_factor;
            jac(i, 2 * (N - 1) + 1) += dab_dXj(1) * vortex_factor;
            jac(i, 0) += dab_dXj(2) * vortex_factor;
            jac(i, 1) += dab_dXj(3) * vortex_factor;
        }
    }

    /* (4) TE factor derivatives. These enter every row i < N, including the
     * no-gap row, through the TE panel influences at node i. */
    const std::vector<TEFactorDeriv> dte = TEFactorDerivs(x);
    if (!dte.empty()) {
        for (int i = 0; i < N; ++i) {
            const Vector2d xi = x.col(i);
            const double a_src = PanelConstSourceStream(XjTE, xi);
            double a_v, b_v;
            PanelLinVortexStream(XjTE, xi, a_v, b_v);
            const double src = a_src * 0.5 * (v(N - 1) - v(0));
            const double vor = (a_v + b_v) * 0.5 * (v(0) - v(N - 1));
            for (const TEFactorDeriv& d : dte)
                jac(i, d.col) += d.tcp * src + d.tdp * vor;
        }
    }
    return jac;
}

MatrixXd AICResidualJacobianX(const FoilState& state, double alpha) {
    /* R = A*gamma - rhs, rhs(i) = -xi_z*cos(alpha) + xi_x*sin(alpha). */
    const int N = state.foil.N;
    MatrixXd jac = AICMatVecJacobianX(state, state.isol.gam);
    const TEInfo te = GetTEInfo(state.foil.x);
    const bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);
    const double cosa = Cosd(alpha);
    const double sina = Sind(alpha);
    for (int i = 0; i < N; ++i) {
        if (nogap && i == N - 1) continue;
        jac(i, 2 * i) -= sina;
        jac(i, 2 * i + 1) += cosa;
    }
    return jac;
}

VectorXd AICResidualSensVJP(const FoilState& state, const VectorXd& lambda,
                            double alpha) {
    const int N = state.foil.N;
    /* Only rows 0..N-1 depend on the coordinates. */
    return AICResidualJacobianX(state, alpha).topRows(N).transpose() *
           lambda.head(N);
}

VectorXd AICMatVecVJP(const FoilState& state, const VectorXd& lambda,
                      const VectorXd& v) {
    const int N = state.foil.N;
    return AICMatVecJacobianX(state, v).topRows(N).transpose() * lambda.head(N);
}

}  // namespace scirocco
