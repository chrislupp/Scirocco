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
*/
#include "scirocco/core/viscous_geom_deriv.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/inviscid_deriv.h"
#include "scirocco/core/panel_influence_deriv.h"
#include "scirocco/core/spline.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

MatrixXd DFoilSDFoilX(const Panel& foil) {
    if (foil.arc_length_mode == ArcLengthMode::Spline)
        return SplineArcLengthDeriv(foil.x);
    if (foil.arc_length_mode == ArcLengthMode::Cumulative)
        return DFoilSDFoilXCumulative(foil);
    return MatrixXd::Zero(foil.N, 2 * foil.N);
}

bool WakeTracedAtOtherAlpha(const FoilState& state) {
    return std::isfinite(state.wake.tracing_alpha) &&
           state.oper.alpha != state.wake.tracing_alpha;
}

/* The state as BuildWake saw it: alpha and gamma at the tracing alpha. */
FoilState WakeTracingState(const FoilState& state) {
    FoilState traced = state;
    if (!WakeTracedAtOtherAlpha(state)) return traced;
    traced.oper.alpha = state.wake.tracing_alpha;
    traced.isol.gam = traced.isol.gamref.col(0) * Cosd(traced.oper.alpha) +
                      traced.isol.gamref.col(1) * Sind(traced.oper.alpha);
    return traced;
}

MatrixXd DGammaDFoilXForWake(const FoilState& state) {
    if (WakeTracedAtOtherAlpha(state)) {
        FoilState traced = WakeTracingState(state);
        return DGammaDFoilX(traced, traced.oper.alpha);
    }
    return DGammaDFoilX(state, state.oper.alpha);
}

MatrixXd DWakeXDFoilX(const FoilState& state) {
    return DWakeXDFoilX(state, DFoilSDFoilX(state.foil));
}
MatrixXd DWakeTDFoilX(const FoilState& state) {
    return DWakeTDFoilX(state, DFoilSDFoilX(state.foil));
}
MatrixXd DWakeXDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    return DWakeXDFoilX(state, dfoil_s_dx, DGammaDFoilXForWake(state));
}
MatrixXd DWakeTDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    return DWakeTDFoilX(state, dfoil_s_dx, DGammaDFoilXForWake(state));
}

/**
 * @brief d(foil.s)/d(foil.x) for the cumulative-panel-norm definition of foil.s
 * (i.e. matches scirocco's RecomputeArcLength).
 *
 * foil.s(0) = 0
 * foil.s(i) = foil.s(i-1) + |x[:,i] - x[:,i-1]|   for i >= 1
 *
 * Let p(i) = |x[:,i] - x[:,i-1]| be the i-th panel length. Then
 *   d(p(i))/d(x[d,k]) = (x[d,i] - x[d,i-1])/p(i)  if k=i
 *                     = -(x[d,i] - x[d,i-1])/p(i) if k=i-1
 *                     = 0                          otherwise
 * and d(foil.s(i))/d(x[d,k]) = sum_{j<=i} d(p(j))/d(x[d,k]).
 *
 * @param foil Panel object containing foil geometry.
 * @return MatrixXd of shape [N x 2N].
 */
MatrixXd DFoilSDFoilXCumulative(const Panel& foil) {
    int N = foil.N;
    MatrixXd out = MatrixXd::Zero(N, 2 * N);
    /* Compute panel-length unit vectors
     * dp_j/dx[d,j] = u[d,j],  dp_j/dx[d,j-1] = -u[d,j]
     * where u[:,j] = (x[:,j] - x[:,j-1]) / |x[:,j] - x[:,j-1]|
     * Then foil.s(i) = sum_{j=1..i} p(j), so
     *   d(foil.s(i))/d(x[d,k])
     *     = u[d,k]   if 1 <= k <= i        (k contributes via p(k) and p(k+1))
     *     wait — let's be careful:
     *   contribution to d(foil.s(i))/d(x[d,k]) from p(j):
     *     = +u[d,j] if k==j and j<=i
     *     = -u[d,j] if k==j-1 and j<=i, i.e. k<=i-1 and j=k+1
     *
     * So for foil.s(i):
     *   if k==0: only p(1) contributes (-u[d,1])
     *   if 1<=k<=i-1: p(k) contributes (+u[d,k]) and p(k+1) contributes
     *   (-u[d,k+1]) if k==i: p(i) contributes (+u[d,i]) if k>i: zero
     * ----------------------------------------------------------------------
     * Compute u[:,j] (the unit tangent of panel j, j ∈ [1,N-1])
     */
    MatrixXd u(2, N); /* u(:,0) unused */
    u.setZero();
    for (int j = 1; j < N; ++j) {
        Vector2d dxv = foil.x.col(j) - foil.x.col(j - 1);
        double pj = dxv.norm();
        if (pj < 1e-14) pj = 1e-14;
        u.col(j) = dxv / pj;
    }

    for (int i = 1; i < N; ++i) {
        for (int j = 1; j <= i; ++j) {
            /* p(j) contributes +u[d,j] to d(foil.s(i))/d(x[d,j])
             * and -u[d,j] to d(foil.s(i))/d(x[d,j-1])
             */
            for (int d = 0; d < 2; ++d) {
                out(i, 2 * j + d) += u(d, j);
                out(i, 2 * (j - 1) + d) -= u(d, j);
            }
        }
    }
    /* row 0 (foil.s(0)=0) is all zeros by construction */
    return out;
}

/**
 * @brief GetTEInfoDeriv — derivatives of TE-info quantities w.r.t.
 * foil.x[k]_dim
 *
 * Returns nonzero only for k ∈ {0, 1, N-2, N-1}. Computes simultaneously:
 *   dtcp, dtdp, dhTE, ddtdx
 *
 * Definitions (matching geometry.cpp:170-210):
 *   t1 = normalize(x[:,0] - x[:,1])
 *   t2 = normalize(x[:,N-1] - x[:,N-2])
 *   tbar = 0.5*(t1+t2);   t = normalize(tbar)
 *   s = x[:,N-1] - x[:,0];   p = s/|s|
 *   hTE = -s(0)*t(1) + s(1)*t(0)
 *   dtdx = t1(0)*t2(1) - t2(0)*t1(1)
 *   tcp = |t(0)*p(1) - t(1)*p(0)|
 *   tdp = t(0)*p(0) + t(1)*p(1)
 *
 * @param x Foil node coordinate matrix [2 x N].
 * @param k Node index (nonzero output only for k ∈ {0, 1, N-2, N-1}).
 * @param dim Coordinate dimension: 0 for x, 1 for y.
 * @param dtcp Output: d(tcp)/d(foil.x[k]_dim).
 * @param dtdp Output: d(tdp)/d(foil.x[k]_dim).
 * @param dhTE Output: d(hTE)/d(foil.x[k]_dim).
 * @param ddtdx Output: d(dtdx)/d(foil.x[k]_dim).
 */
void GetTEInfoDeriv(const MatrixXd& x, int k, int dim, double& dtcp,
                    double& dtdp, double& dhTE, double& ddtdx) {
    int N = static_cast<int>(x.cols());
    dtcp = 0.0;
    dtdp = 0.0;
    dhTE = 0.0;
    ddtdx = 0.0;
    if (k != 0 && k != 1 && k != N - 2 && k != N - 1) return;

    Vector2d v1 = x.col(0) - x.col(1);
    Vector2d v2 = x.col(N - 1) - x.col(N - 2);
    double v1n = v1.norm();
    double v2n = v2.norm();
    if (v1n < 1e-12 || v2n < 1e-12) return;
    Vector2d t1 = v1 / v1n;
    Vector2d t2 = v2 / v2n;

    Vector2d tbar = 0.5 * (t1 + t2);
    double tbn = tbar.norm();
    if (tbn < 1e-12) return;
    Vector2d t = tbar / tbn;

    Vector2d s = x.col(N - 1) - x.col(0);
    double sn = s.norm();
    if (sn < 1e-12) return;
    Vector2d p = s / sn;

    /* d(v1)/d(x[k]_dim), d(v2)/d(x[k]_dim), d(s)/d(x[k]_dim) */
    Vector2d dv1 = Vector2d::Zero();
    Vector2d dv2 = Vector2d::Zero();
    Vector2d ds = Vector2d::Zero();
    if (k == 0)
        dv1(dim) = 1.0;
    else if (k == 1)
        dv1(dim) = -1.0;
    if (k == N - 1)
        dv2(dim) = 1.0;
    else if (k == N - 2)
        dv2(dim) = -1.0;
    if (k == N - 1)
        ds(dim) = 1.0;
    else if (k == 0)
        ds(dim) = -1.0;

    auto normalize_jac = [](const Vector2d& u, double un) -> Matrix2d {
        return (Matrix2d::Identity() - u * u.transpose()) / un;
    };

    Vector2d dt1 = normalize_jac(t1, v1n) * dv1;
    Vector2d dt2 = normalize_jac(t2, v2n) * dv2;
    Vector2d dtbar = 0.5 * (dt1 + dt2);
    Vector2d dt = normalize_jac(t, tbn) * dtbar;
    Vector2d dp = normalize_jac(p, sn) * ds;

    /* dhTE = -ds(0)*t(1) - s(0)*dt(1) + ds(1)*t(0) + s(1)*dt(0) */
    dhTE = -ds(0) * t(1) - s(0) * dt(1) + ds(1) * t(0) + s(1) * dt(0);

    /* ddtdx = dt1(0)*t2(1) + t1(0)*dt2(1) - dt2(0)*t1(1) - t2(0)*dt1(1) */
    ddtdx = dt1(0) * t2(1) + t1(0) * dt2(1) - dt2(0) * t1(1) - t2(0) * dt1(1);

    /* tcp = |t(0)*p(1) - t(1)*p(0)| */
    double cross = t(0) * p(1) - t(1) * p(0);
    double sgn = (cross >= 0) ? 1.0 : -1.0;
    double dcross = dt(0) * p(1) + t(0) * dp(1) - dt(1) * p(0) - t(1) * dp(0);
    dtcp = sgn * dcross;

    /* tdp = t·p */
    dtdp = dt.dot(p) + t.dot(dp);
}

/* (DGammaAtNodeDFoilX helper removed — was used by the old gamma-based
 * DSstagDFoilX chain. Production now uses ue-based xi which doesn't need
 * a per-node gamma derivative.)
 */

/**
 * @brief DSstagDFoilX — ue-based xi convention (production state after commit
 * 4051c2d).
 *
 *   sstag = w1 * S[i0] + w2 * S[i1] where weights come from StagPointMove:
 *     w1 = ue[i1] / (ue[i0] + ue[i1])
 *     w2 = ue[i0] / (ue[i0] + ue[i1])
 *
 *   At fixed U (the adjoint's partial), ue at Istag is held fixed → only
 *   the foil.s chain contributes:
 *     ∂sstag/∂foil.x = w1 * ∂S[i0]/∂foil.x + w2 * ∂S[i1]/∂foil.x
 *
 *   The previous implementation used gamma-based weights and added a
 *   ∂sstag/∂gam·∂gam/∂foil.x chain — correct for gamma-based xi only.
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return VectorXd of length 2N.
 */
VectorXd DSstagDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int twoN = 2 * N;
    int i0 = state.isol.Istag[0];
    int i1 = state.isol.Istag[1];
    double ue0 = state.glob.U(3, i0);
    double ue1 = state.glob.U(3, i1);
    double den = ue0 + ue1;
    double w1 = ue1 / den;
    double w2 = ue0 / den;

    VectorXd dsstag = VectorXd::Zero(twoN);
    dsstag += w1 * dfoil_s_dx.row(i0).transpose();
    dsstag += w2 * dfoil_s_dx.row(i1).transpose();
    return dsstag;
}

/**
 * @brief DWakeSDFoilX
 *
 * wake.s[i] = foil.s[N-1] + sv[i]
 *   sv[0] = 0
 *   sv[i>0] from SpaceGeom(ds1, wakelen*chord, Nw) — depends on
 *     ds1 = 0.5 * ((foil.s[1] - foil.s[0]) + (foil.s[N-1] - foil.s[N-2]))
 *     and on chord (assumed constant in foil.x; if chord changes, add chain).
 *
 * Therefore d(wake.s[i])/d(foil.x) =
 *   d(foil.s[N-1])/d(foil.x) + (d(sv[i])/d(ds1)) * d(ds1)/d(foil.x)
 *
 * We obtain d(sv[i])/d(ds1) numerically via one extra SpaceGeom-style call
 * (analytic derivative through the Newton iteration is straightforward —
 * see below).
 *
 * Geometric spacing: sv[i] = sum_{j=1..i} ds_{j-1} * r^{j-1}, where r solves
 *   L = ds1 * (r^N - 1)/(r - 1),    N = Np - 1
 * Implicit differentiation:
 *   0 = (r^N - 1)/(r-1) * ddS + ds1 * [N*r^(N-1)/(r-1) - (r^N-1)/(r-1)^2] dr
 * => dr/dds1 = -((r^N-1)/(r-1)) / { ds1 * [N*r^(N-1)/(r-1) - (r^N-1)/(r-1)^2] }
 * Then dsv[i]/dds1 = sum_{j=0..i-1} (r^j + ds1 * j*r^(j-1) * dr/dds1).
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nw x 2N].
 */
namespace {
/* Compute (sv, dsv_dds1) given Np, ds1, L (matching wake.cpp SpaceGeom). */
void SpaceGeomWithDeriv(double ds1, double L, int Np, VectorXd& sv,
                        VectorXd& dsv_dds1) {
    if (Np <= 1) {
        sv = VectorXd::Zero(1);
        dsv_dds1 = VectorXd::Zero(1);
        return;
    }
    int N = Np - 1;
    double d = L / ds1;

    /* Initial guess for r (same as wake.cpp) */
    double a = N * (N - 1.0) * (N - 2.0) / 6.0;
    double b = N * (N - 1.0) / 2.0;
    double c = N - d;
    double disc = std::max(b * b - 4.0 * a * c, 0.0);
    double r = 1.0 + (-b + std::sqrt(disc)) / (2.0 * a);
    for (int k = 0; k < 20; ++k) {
        double rN = std::pow(r, N);
        double R = rN - 1.0 - d * (r - 1.0);
        double R_r = N * std::pow(r, N - 1) - d;
        double dr = -R / R_r;
        if (std::abs(dr) < 1e-14) break;
        r += dr;
    }

    /* dr/dds1: chord and L are constants → d/dds1 of  rN - 1 - (L/ds1)*(r-1) =
     * 0
     *   N r^{N-1} dr + L/ds1^2 * (r-1) - L/ds1 * dr = 0
     *   dr * (N r^{N-1} - L/ds1) = -L * (r-1) / ds1^2
     *   dr_dds1 = -L * (r-1) / (ds1^2 * (N r^{N-1} - L/ds1))
     */
    double dr_dds1 =
        -L * (r - 1.0) / (ds1 * ds1 * (N * std::pow(r, N - 1) - L / ds1));

    sv = VectorXd::Zero(Np);
    dsv_dds1 = VectorXd::Zero(Np);
    double dx = ds1;
    double ddx_dds1 = 1.0;
    for (int i = 1; i < Np; ++i) {
        sv(i) = sv(i - 1) + dx;
        dsv_dds1(i) = dsv_dds1(i - 1) + ddx_dds1;
        /* Next-interval update: dx <- dx * r */
        ddx_dds1 = ddx_dds1 * r + dx * dr_dds1;
        dx *= r;
    }
}
} /* namespace */

MatrixXd DWakeSDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int twoN = 2 * N;
    double wakelen = state.geom.wakelen;
    double chord = state.geom.chord;

    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;

    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    /* d(ds1)/d(foil.x) = 0.5 * ( d(foil.s[1])/d(foil.x) -
     * d(foil.s[0])/d(foil.x)
     *                          + d(foil.s[N-1])/d(foil.x)
     *                          - d(foil.s[N-2])/d(foil.x) )
     */
    VectorXd dds1 =
        0.5 *
        (dfoil_s_dx.row(1).transpose() - dfoil_s_dx.row(0).transpose() +
         dfoil_s_dx.row(N - 1).transpose() - dfoil_s_dx.row(N - 2).transpose());

    /* d(foil.s[N-1])/d(foil.x) */
    VectorXd dSTE = dfoil_s_dx.row(N - 1).transpose();

    MatrixXd dwake_s = MatrixXd::Zero(Nw, twoN);
    for (int i = 0; i < Nw; ++i) {
        /* d(wake.s[i])/d(foil.x) = dSTE + dsv_dds1[i] * dds1 */
        dwake_s.row(i) = (dSTE + dsv_dds1(i) * dds1).transpose();
    }
    return dwake_s;
}

/**
 * @brief DXiDFoilX
 *   xi[i]    = |foil.s[i] - sstag|      for i ∈ [0, N)
 *   xi[N+j]  = wake.s[j] - sstag        for j ∈ [0, Nw)
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nsys x 2N].
 */
MatrixXd DXiDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    int twoN = 2 * N;

    VectorXd dsstag = DSstagDFoilX(state, dfoil_s_dx);  /* [2N] */
    MatrixXd dwake_s = DWakeSDFoilX(state, dfoil_s_dx); /* [Nw x 2N] */

    MatrixXd dxi = MatrixXd::Zero(Nsys, twoN);
    double sstag = state.isol.sstag;

    for (int i = 0; i < N; ++i) {
        double diff = state.foil.s(i) - sstag;
        double sgn = (diff >= 0.0) ? 1.0 : -1.0;
        /* d(xi[i])/d(foil.x) = sgn * (d(foil.s[i])/d(foil.x) -
         * d(sstag)/d(foil.x))
         */
        dxi.row(i) = sgn * (dfoil_s_dx.row(i) - dsstag.transpose());
    }
    for (int j = 0; j < Nw; ++j) {
        dxi.row(N + j) = dwake_s.row(j) - dsstag.transpose();
    }
    return dxi;
}

/**
 * @brief DWgapDFoilX
 *   wgap[j] depends on hTE, dtdx, xib = (xi[N+j] - xi[N])/Lw,  Lw = flen*hTE
 *
 *   ξib = (xi[N+j] - xi[N]) / Lw
 *       = ((wake.s[j] - sstag) - (foil.s[N-1] - sstag)) / Lw   if foil.s[N-1]
 *       >= sstag
 *
 *   But wait — xi[N]  = |foil.s[N-1] - sstag|. Assuming sstag < foil.s[N-1]
 *   on the upper-surface side, xi[N] = foil.s[N-1] - sstag, and
 *   xi[N+j] - xi[N] = wake.s[j] - foil.s[N-1] = sv[j]. So sstag drops out
 *   and the chain is purely through ds1 and chord.
 *
 *   wgap[j] = hTE * (1 + (2 + flen*dtdx)*xib) * (1-xib)^2   if xib<=1 else 0
 *
 *   Let A = 1 + (2 + flen*dtdx)*xib,  B = (1-xib)^2
 *   wgap = hTE * A * B
 *   d(wgap)/d(.) = dhTE * A*B
 *                + hTE * dA * B
 *                + hTE * A * dB
 *   dA = flen*ddtdx * xib + (2 + flen*dtdx) * dxib
 *   dB = -2*(1-xib) * dxib
 *   dxib = ( d(xi[N+j]-xi[N])/d(.) * Lw - (xi[N+j]-xi[N]) * dLw ) / Lw^2
 *        = d(sv[j])/d(.)/Lw - xib * dLw/Lw
 *   dLw = flen * dhTE
 *
 * @param state FoilState containing the current solution and geometry.
 * @param dfoil_s_dx d(foil.s)/d(foil.x) matrix of shape [N x 2N].
 * @return MatrixXd of shape [Nw x 2N].
 */
MatrixXd DWgapDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int twoN = 2 * N;
    double wakelen = state.geom.wakelen;
    double chord = state.geom.chord;
    const double flen = 2.5;

    /* Recompute TE info (matches SetWakeGap) */
    TEInfo te = GetTEInfo(state.foil.x);
    double hTE = te.hTE;
    double dtdx_raw = te.dtdx;
    double clip_lo = -3.0 / flen;
    double clip_hi = 3.0 / flen;
    bool clipped = (dtdx_raw < clip_lo) || (dtdx_raw > clip_hi);
    double dtdx = std::min(std::max(dtdx_raw, clip_lo), clip_hi);
    double Lw = std::max(flen * hTE, 1e-10);

    /* d(hTE), d(dtdx_raw) for each (k, dim) over TE-affecting nodes
     * We'll need 4 nodes × 2 dims = 8 entries → store sparsely.
     */
    auto te_deriv_col = [&](int k, int dim, double& dhTE, double& ddtdx) {
        double dtcp, dtdp, dhTE_loc, ddtdx_loc;
        GetTEInfoDeriv(state.foil.x, k, dim, dtcp, dtdp, dhTE_loc, ddtdx_loc);
        dhTE = dhTE_loc;
        ddtdx = clipped ? 0.0 : ddtdx_loc;
    };

    /* d(sv[j])/d(foil.x) via SpaceGeom analytical sensitivity */
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    VectorXd dds1 =
        0.5 *
        (dfoil_s_dx.row(1).transpose() - dfoil_s_dx.row(0).transpose() +
         dfoil_s_dx.row(N - 1).transpose() - dfoil_s_dx.row(N - 2).transpose());

    MatrixXd dwgap = MatrixXd::Zero(Nw, twoN);

    for (int j = 0; j < Nw; ++j) {
        double xib = sv(j) / Lw;
        if (xib > 1.0) {
            /* wgap[j] = 0  outside the support; derivative is 0 too. */
            continue;
        }
        double A = 1.0 + (2.0 + flen * dtdx) * xib;
        double B = (1.0 - xib) * (1.0 - xib);

        /* For each (k, dim), assemble d(wgap[j])/d(x[k]_dim) */
        for (int k = 0; k < N; ++k) {
            for (int dim = 0; dim < 2; ++dim) {
                int col = 2 * k + dim;

                double dhTE = 0.0;
                double ddtdx = 0.0;
                if (k == 0 || k == 1 || k == N - 2 || k == N - 1) {
                    te_deriv_col(k, dim, dhTE, ddtdx);
                }
                double dLw = flen * dhTE;
                double dsv = dsv_dds1(j) * dds1(col);
                /* xib = sv/Lw → dxib = dsv/Lw - sv/Lw^2 * dLw = dsv/Lw -
                 * xib*dLw/Lw
                 */
                double dxib = dsv / Lw - xib * dLw / Lw;
                double dA = flen * ddtdx * xib + (2.0 + flen * dtdx) * dxib;
                double dB = -2.0 * (1.0 - xib) * dxib;
                double dwj = dhTE * A * B + hTE * dA * B + hTE * A * dB;
                dwgap(j, col) = dwj;
            }
        }
    }
    return dwgap;
}

/**
 * @brief Helper: full forward Jacobian d(gamma)/d(foil.x) [(N+1) x 2N].
 *
 * gamma satisfies A(foil.x) * gamma = rhs(alpha, foil.x), so
 *   d(gamma)/d(foil.x) = -A^{-1} * dR/d(foil.x) where R = A*gamma - rhs.
 * dR/d(foil.x) is assembled once (O(N^2) panel derivatives), then one LU
 * solve with 2N right-hand sides gives all of d(gamma)/d(foil.x).
 *
 * @param state FoilState containing the current solution and geometry.
 * @param alpha Angle of attack in degrees.
 * @return MatrixXd of shape [N+1 x 2N].
 */
MatrixXd DGammaDFoilX(const FoilState& state, double alpha) {
    return -SolveAIC(state, AICResidualJacobianX(state, alpha));
}

/**
 * @brief Helper: InviscidVelocity at xi with full forward sensitivity matrix.
 * Returns V, ∂V/∂xi [2x2], and the combined Jacobian
 *   ∂V/∂(foil.x_panels) at fixed xi  +  V_G * d(gamma)/d(foil.x)
 * of shape [2 x 2N].
 */
/* ============================================================================
 */
VelDeriv InviscidVelocityWithDerivExt(const FoilState& state,
                                      const Vector2d& xi,
                                      const VectorXd& vortex,
                                      const MatrixXd& dvortex_dx, double alpha);

VelDeriv InviscidVelocityWithDeriv(const FoilState& state, const Vector2d& xi,
                                   const MatrixXd& dgam_dx, double alpha) {
    return InviscidVelocityWithDerivExt(state, xi, state.isol.gam, dgam_dx,
                                        alpha);
}

VelDeriv InviscidVelocityWithDerivExt(const FoilState& state,
                                      const Vector2d& xi,
                                      const VectorXd& vortex,
                                      const MatrixXd& dvortex_dx,
                                      double alpha) {
    int N = state.foil.N;
    int twoN = 2 * N;
    const MatrixXd& X = state.foil.x;
    const VectorXd& G = vortex;
    const MatrixXd& dgam_dx = dvortex_dx;
    double Vinf = state.oper.Vinf;

    VelDeriv out;
    out.V.setZero();
    out.V_xi.setZero();
    out.V_xfoil = MatrixXd::Zero(2, twoN);
    out.V_G = MatrixXd::Zero(2, N + 1);

    bool has_gam_chain = (dgam_dx.cols() == twoN);

    /* Helper lambda: write a 2-vector contribution into V_xfoil columns
     * [2*k + 0, 2*k + 1] from a 2x2 a_Xj block (cols 0,1 = x1=Xj[k0], cols 2,3
     * = x2=Xj[k1]).
     */
    auto add_Xj_block = [&](const Eigen::Matrix<double, 2, 4>& M_Xj, int k0,
                            int k1) {
        out.V_xfoil.col(2 * k0 + 0) += M_Xj.col(0);
        out.V_xfoil.col(2 * k0 + 1) += M_Xj.col(1);
        out.V_xfoil.col(2 * k1 + 0) += M_Xj.col(2);
        out.V_xfoil.col(2 * k1 + 1) += M_Xj.col(3);
    };

    /* ---- Linear vortex panel contributions (j → j+1) ---- */
    for (int j = 0; j < N - 1; ++j) {
        Matrix2d Xj;
        Xj.col(0) = X.col(j);
        Xj.col(1) = X.col(j + 1);

        LinVortexVelocityDeriv d = PanelLinVortexVelocityDeriv(Xj, xi);
        out.V += d.a * G(j) + d.b * G(j + 1);
        out.V_xi += d.a_xi * G(j) + d.b_xi * G(j + 1);
        out.V_G.col(j) += d.a;
        out.V_G.col(j + 1) += d.b;
        /* Direct panel-geom contribution (G held fixed) */
        add_Xj_block(G(j) * d.a_Xj + G(j + 1) * d.b_Xj, j, j + 1);
        /* Gamma chain: V_G contribution = a * dgam_dx.row(j) + b *
         * dgam_dx.row(j+1)
         */
        if (has_gam_chain) {
            out.V_xfoil += d.a * dgam_dx.row(j) + d.b * dgam_dx.row(j + 1);
        }
    }

    /* ---- TE info ---- */
    TEInfo te = GetTEInfo(X);
    double tcp = te.tcp;
    double tdp = te.tdp;

    /* ---- TE constant source panel (Xj = [X[N-1], X[0]]) ---- */
    Matrix2d XjTE;
    XjTE.col(0) = X.col(N - 1);
    XjTE.col(1) = X.col(0);
    ConstSourceVelocityDeriv ds = PanelConstSourceVelocityDeriv(XjTE, xi);
    Vector2d v_src = ds.val;
    double sigma_te = 0.5 * tcp * (G(N - 1) - G(0));
    out.V += v_src * sigma_te;
    out.V_xi += ds.val_xi * sigma_te;
    out.V_G.col(N - 1) += 0.5 * tcp * v_src;
    out.V_G.col(0) -= 0.5 * tcp * v_src;
    /* Direct geom of THIS panel: d(v_src)/dXj * sigma_te + v_src *
     * d(sigma_te)/dXj_te_nodes d(sigma_te)/dXj of this panel =
     * 0.5*(G(N-1)-G(0)) * d(tcp)/dXj_te_via_GetTEInfo BUT tcp depends on the
     * global TE info — handled below as a separate global-TE contribution.
     */
    add_Xj_block(sigma_te * ds.val_Xj, N - 1, 0);
    /* Gamma chain: v_src * 0.5*tcp * (dgam[N-1] - dgam[0]) */
    if (has_gam_chain) {
        out.V_xfoil +=
            v_src * (0.5 * tcp * (dgam_dx.row(N - 1) - dgam_dx.row(0)));
    }

    /* ---- TE vortex panel ---- */
    LinVortexVelocityDeriv dv = PanelLinVortexVelocityDeriv(XjTE, xi);
    Vector2d f_vor = dv.a + dv.b;
    double tdp_factor = 0.5 * tdp * (G(0) - G(N - 1));
    out.V += f_vor * tdp_factor;
    out.V_xi += (dv.a_xi + dv.b_xi) * tdp_factor;
    out.V_G.col(0) += 0.5 * tdp * f_vor;
    out.V_G.col(N - 1) -= 0.5 * tdp * f_vor;
    /* Direct geom on Xj of this panel (G held fixed) */
    add_Xj_block(tdp_factor * (dv.a_Xj + dv.b_Xj), N - 1, 0);
    /* Gamma chain: f_vor * 0.5*tdp * (dgam[0] - dgam[N-1]) */
    if (has_gam_chain) {
        out.V_xfoil +=
            f_vor * (0.5 * tdp * (dgam_dx.row(0) - dgam_dx.row(N - 1)));
    }

    /* ---- Global TE info dependence on foil.x[{0,1,N-2,N-1}] ----
     * v_src * sigma_te uses tcp → d(tcp)/d(foil.x[k]_dim) feeds in for k ∈
     * {0,1,N-2,N-1} f_vor * 0.5*tdp*(G(0)-G(N-1)) uses tdp similarly.
     */
    {
        std::vector<int> te_nodes = {0, 1, N - 2, N - 1};
        for (int k : te_nodes) {
            for (int dim = 0; dim < 2; ++dim) {
                double dtcp, dtdp, dhTE, ddtdx;
                GetTEInfoDeriv(X, k, dim, dtcp, dtdp, dhTE, ddtdx);
                int col = 2 * k + dim;
                /* d(sigma_te)/d(tcp) * dtcp = 0.5*(G(N-1)-G(0))*dtcp */
                double dsig = 0.5 * (G(N - 1) - G(0)) * dtcp;
                out.V_xfoil.col(col) += v_src * dsig;
                /* d(tdp_factor)/d(tdp) * dtdp = 0.5*(G(0)-G(N-1))*dtdp */
                double dtdp_factor = 0.5 * (G(0) - G(N - 1)) * dtdp;
                out.V_xfoil.col(col) += f_vor * dtdp_factor;
            }
        }
    }

    /* ---- Freestream ---- */
    out.V(0) += Vinf * Cosd(alpha);
    out.V(1) += Vinf * Sind(alpha);
    /* No foil.x or xi dependence. */

    return out;
}

/**
 * @brief DWakeXDFoilX — full forward Jacobian of wake.x w.r.t. foil.x via
 * linearised predictor-corrector sweep.
 *
 * Returns [2*Nw x 2N], row 2*j+d = d(wake.x(d, j))/d(foil.x).
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2*Nw x 2N].
 */
MatrixXd DWakeXDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx,
                      const MatrixXd& dgam_dx_wake) {
    if (WakeTracedAtOtherAlpha(state)) {
        return DWakeXDFoilX(WakeTracingState(state), dfoil_s_dx, dgam_dx_wake);
    }

    int N = state.foil.N;
    int Nw = state.wake.N;
    int twoN = 2 * N;
    double alpha = state.oper.alpha;
    double chord = state.geom.chord;
    double wakelen = state.geom.wakelen;

    /* d(gamma)/d(foil.x) at the tracing alpha (== alpha here). */
    const MatrixXd& dgam_dx = dgam_dx_wake;

    /* SpaceGeom sensitivities matching scirocco's wake.cpp. */
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    VectorXd dds1_dx =
        0.5 *
        (dfoil_s_dx.row(1).transpose() - dfoil_s_dx.row(0).transpose() +
         dfoil_s_dx.row(N - 1).transpose() - dfoil_s_dx.row(N - 2).transpose());

    /* ---- Initial wake.x[:, 0] ----
     * x_te = 0.5*(foil.x[:,0] + foil.x[:,N-1])
     * n_te = foil.x[:,N-1] - foil.x[:,0]
     * t_te_raw = (n_te(1), -n_te(0)), then sign-flipped so t_te[0] > 0
     * wake.x[:, 0] = x_te + 1e-5 * chord * t_te
     */
    Vector2d xy1 = state.foil.x.col(0);
    Vector2d xyN = state.foil.x.col(N - 1);
    Vector2d n_te = xyN - xy1;
    Vector2d t_te_raw(n_te(1), -n_te(0));
    double sign = (t_te_raw(0) < 0) ? -1.0 : 1.0;

    MatrixXd dxyw_i = MatrixXd::Zero(2, twoN);
    /* Average term 0.5*(foil.x[:,0] + foil.x[:,N-1]) */
    for (int dim = 0; dim < 2; ++dim) {
        dxyw_i(dim, 2 * 0 + dim) += 0.5;
        dxyw_i(dim, 2 * (N - 1) + dim) += 0.5;
    }
    /* Tangent term: 1e-5 * chord * sign * (n_te[1], -n_te[0])
     * d(n_te[0])/d(foil.x[d',k]) = δ(d',0) * (δ(k,N-1) - δ(k,0))
     * d(n_te[1])/d(foil.x[d',k]) = δ(d',1) * (δ(k,N-1) - δ(k,0))
     * So d(t_te[0])/d(foil.x[1, N-1]) = +sign*1, d(t_te[0])/d(foil.x[1, 0]) =
     * -sign*1
     *    d(t_te[1])/d(foil.x[0, N-1]) = -sign*1, d(t_te[1])/d(foil.x[0, 0]) =
     *    +sign*1
     */
    double factor = 1e-5 * chord * sign;
    dxyw_i(0, 2 * (N - 1) + 1) += factor;
    dxyw_i(0, 2 * 0 + 1) -= factor;
    dxyw_i(1, 2 * (N - 1) + 0) -= factor;
    dxyw_i(1, 2 * 0 + 0) += factor;

    /* Layout: row d*Nw + j = d(wake.x(d, j))/d(foil.x) (matches the
     * Python extractor's np.reshape(dwake_x[2, Nw, 2N], (2*Nw, 2N)) in C
     * order).
     */
    MatrixXd dwake_x = MatrixXd::Zero(2 * Nw, twoN);
    dwake_x.row(0 * Nw + 0) = dxyw_i.row(0);
    dwake_x.row(1 * Nw + 0) = dxyw_i.row(1);

    Vector2d xyw_i = state.wake.x.col(0);

    /* ---- Predictor-corrector sweep, linearised ---- */
    for (int i = 0; i < Nw - 1; ++i) {
        double ds_i = sv(i + 1) - sv(i);
        VectorXd dds_i_dx =
            (dsv_dds1(i + 1) - dsv_dds1(i)) * dds1_dx; /* [2N] */

        /* V0 at xyw_i */
        VelDeriv vd0 = InviscidVelocityWithDeriv(state, xyw_i, dgam_dx, alpha);
        Vector2d V0 = vd0.V;
        double V0n = V0.norm();
        if (V0n < 1e-12) V0n = 1e-12;
        Vector2d v1 = V0 / V0n;
        Matrix2d norm_jac1 = (Matrix2d::Identity() - v1 * v1.transpose()) / V0n;
        MatrixXd dV0_dx = vd0.V_xfoil + vd0.V_xi * dxyw_i; /* [2 x 2N] */
        MatrixXd dv1_dx = norm_jac1 * dV0_dx;

        /* Predictor step */
        Vector2d x_pred = xyw_i + ds_i * v1;
        MatrixXd dx_pred =
            dxyw_i + v1 * dds_i_dx.transpose() + ds_i * dv1_dx; /* [2 x 2N] */

        /* V1 at x_pred */
        VelDeriv vd1 = InviscidVelocityWithDeriv(state, x_pred, dgam_dx, alpha);
        Vector2d V1 = vd1.V;
        double V1n = V1.norm();
        if (V1n < 1e-12) V1n = 1e-12;
        Vector2d v2 = V1 / V1n;
        Matrix2d norm_jac2 = (Matrix2d::Identity() - v2 * v2.transpose()) / V1n;
        MatrixXd dV1_dx = vd1.V_xfoil + vd1.V_xi * dx_pred;
        MatrixXd dv2_dx = norm_jac2 * dV1_dx;

        /* Corrector */
        Vector2d xyw_next = xyw_i + ds_i * 0.5 * (v1 + v2);
        MatrixXd dxyw_next = dxyw_i + (0.5 * (v1 + v2)) * dds_i_dx.transpose() +
                             ds_i * 0.5 * (dv1_dx + dv2_dx);

        dwake_x.row(0 * Nw + (i + 1)) = dxyw_next.row(0);
        dwake_x.row(1 * Nw + (i + 1)) = dxyw_next.row(1);

        xyw_i = xyw_next;
        dxyw_i = dxyw_next;
    }

    return dwake_x;
}

/**
 * @brief DWakeTDFoilX — derivative of wake.t (unit tangents) w.r.t. foil.x.
 *
 * Mirrors DWakeXDFoilX's iterative sweep. wake.t[j] is set during BuildWake
 * by:
 *   for i = 0..Nw-2:
 *     v1 = V(wake.x[i]) / |V|;    wake.t[i] = v1   (overwrites prior v2)
 *     v2 = V(x_pred) / |V|;       wake.t[i+1] = v2 (overwritten next iter
 *                                                   except for i = Nw-2)
 * So:
 *   wake.t[j]  (j < Nw-1) = v1 in iteration i=j → derivative = dv1_dx_at_i=j
 *   wake.t[Nw-1]          = v2 in iteration i=Nw-2 → derivative = dv2_dx
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2*Nw x 2N].
 */
MatrixXd DWakeTDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx,
                      const MatrixXd& dgam_dx_wake) {
    if (WakeTracedAtOtherAlpha(state)) {
        return DWakeTDFoilX(WakeTracingState(state), dfoil_s_dx, dgam_dx_wake);
    }

    int N = state.foil.N;
    int Nw = state.wake.N;
    int twoN = 2 * N;
    double alpha = state.oper.alpha;
    double chord = state.geom.chord;
    double wakelen = state.geom.wakelen;

    const MatrixXd& dgam_dx = dgam_dx_wake;

    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    VectorXd dds1_dx =
        0.5 *
        (dfoil_s_dx.row(1).transpose() - dfoil_s_dx.row(0).transpose() +
         dfoil_s_dx.row(N - 1).transpose() - dfoil_s_dx.row(N - 2).transpose());

    /* Initial wake.x[:, 0] sensitivities (same as DWakeXDFoilX) */
    Vector2d xy1 = state.foil.x.col(0);
    Vector2d xyN = state.foil.x.col(N - 1);
    Vector2d n_te = xyN - xy1;
    Vector2d t_te_raw(n_te(1), -n_te(0));
    double sign = (t_te_raw(0) < 0) ? -1.0 : 1.0;

    MatrixXd dxyw_i = MatrixXd::Zero(2, twoN);
    for (int dim = 0; dim < 2; ++dim) {
        dxyw_i(dim, 2 * 0 + dim) += 0.5;
        dxyw_i(dim, 2 * (N - 1) + dim) += 0.5;
    }
    double factor = 1e-5 * chord * sign;
    dxyw_i(0, 2 * (N - 1) + 1) += factor;
    dxyw_i(0, 2 * 0 + 1) -= factor;
    dxyw_i(1, 2 * (N - 1) + 0) -= factor;
    dxyw_i(1, 2 * 0 + 0) += factor;

    MatrixXd dwake_t = MatrixXd::Zero(2 * Nw, twoN);
    Vector2d xyw_i = state.wake.x.col(0);

    for (int i = 0; i < Nw - 1; ++i) {
        double ds_i = sv(i + 1) - sv(i);
        VectorXd dds_i_dx = (dsv_dds1(i + 1) - dsv_dds1(i)) * dds1_dx;

        /* v1 at xyw_i  → wake.t[i] */
        VelDeriv vd0 = InviscidVelocityWithDeriv(state, xyw_i, dgam_dx, alpha);
        Vector2d V0 = vd0.V;
        double V0n = V0.norm();
        if (V0n < 1e-12) V0n = 1e-12;
        Vector2d v1 = V0 / V0n;
        Matrix2d norm_jac1 = (Matrix2d::Identity() - v1 * v1.transpose()) / V0n;
        MatrixXd dV0_dx = vd0.V_xfoil + vd0.V_xi * dxyw_i;
        MatrixXd dv1_dx = norm_jac1 * dV0_dx;

        /* Store wake.t[i] = v1 */
        dwake_t.row(0 * Nw + i) = dv1_dx.row(0);
        dwake_t.row(1 * Nw + i) = dv1_dx.row(1);

        /* Predictor step */
        Vector2d x_pred = xyw_i + ds_i * v1;
        MatrixXd dx_pred = dxyw_i + v1 * dds_i_dx.transpose() + ds_i * dv1_dx;

        /* v2 at x_pred → wake.t[i+1] (overwritten except for last iter) */
        VelDeriv vd1 = InviscidVelocityWithDeriv(state, x_pred, dgam_dx, alpha);
        Vector2d V1 = vd1.V;
        double V1n = V1.norm();
        if (V1n < 1e-12) V1n = 1e-12;
        Vector2d v2 = V1 / V1n;
        Matrix2d norm_jac2 = (Matrix2d::Identity() - v2 * v2.transpose()) / V1n;
        MatrixXd dV1_dx = vd1.V_xfoil + vd1.V_xi * dx_pred;
        MatrixXd dv2_dx = norm_jac2 * dV1_dx;

        /* Tentatively store wake.t[i+1] = v2 (overwritten unless last) */
        dwake_t.row(0 * Nw + (i + 1)) = dv2_dx.row(0);
        dwake_t.row(1 * Nw + (i + 1)) = dv2_dx.row(1);

        /* Update for next iter (corrector → xyw_{i+1}) */
        Vector2d xyw_next = xyw_i + ds_i * 0.5 * (v1 + v2);
        MatrixXd dxyw_next = dxyw_i + (0.5 * (v1 + v2)) * dds_i_dx.transpose() +
                             ds_i * 0.5 * (dv1_dx + dv2_dx);
        xyw_i = xyw_next;
        dxyw_i = dxyw_next;
    }

    return dwake_t;
}

WakeSweepLinearization LinearizeWakeSweep(const FoilState& state,
                                          const MatrixXd& dfoil_s_dx) {
    if (WakeTracedAtOtherAlpha(state)) {
        WakeSweepLinearization lin =
            LinearizeWakeSweep(WakeTracingState(state), dfoil_s_dx);
        lin.traced_at_other_alpha = true;
        return lin;
    }

    /* Mirrors DWakeXDFoilX / DWakeTDFoilX: same primal sweep, but stores the
     * local derivatives instead of propagating 2N tangent directions. */
    const int N = state.foil.N;
    const int Nw = state.wake.N;
    const int twoN = 2 * N;
    const double alpha = state.oper.alpha;
    const double chord = state.geom.chord;

    WakeSweepLinearization lin;
    lin.N = N;
    lin.Nw = Nw;

    const double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                              (state.foil.s(N - 1) - state.foil.s(N - 2)));
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, state.geom.wakelen * chord, Nw, sv, dsv_dds1);
    lin.dds1_dx =
        0.5 *
        (dfoil_s_dx.row(1).transpose() - dfoil_s_dx.row(0).transpose() +
         dfoil_s_dx.row(N - 1).transpose() - dfoil_s_dx.row(N - 2).transpose());

    /* wake.x[:, 0] = TE midpoint + 1e-5 * chord * (rotated TE gap), see
     * DWakeXDFoilX for the derivation. */
    const Vector2d n_te = state.foil.x.col(N - 1) - state.foil.x.col(0);
    const double sign = (n_te(1) < 0) ? -1.0 : 1.0; /* t_te_raw(0) = n_te(1) */
    lin.D0 = MatrixXd::Zero(2, twoN);
    for (int dim = 0; dim < 2; ++dim) {
        lin.D0(dim, 2 * 0 + dim) += 0.5;
        lin.D0(dim, 2 * (N - 1) + dim) += 0.5;
    }
    const double factor = 1e-5 * chord * sign;
    lin.D0(0, 2 * (N - 1) + 1) += factor;
    lin.D0(0, 2 * 0 + 1) -= factor;
    lin.D0(1, 2 * (N - 1) + 0) -= factor;
    lin.D0(1, 2 * 0 + 0) += factor;

    const int nsteps = std::max(Nw - 1, 0);
    lin.ds.resize(nsteps);
    lin.dds_dds1.resize(nsteps);
    lin.v1.resize(nsteps);
    lin.v2.resize(nsteps);
    lin.J1.resize(nsteps);
    lin.J2.resize(nsteps);
    lin.Vxi0.resize(nsteps);
    lin.Vxi1.resize(nsteps);
    lin.P0.resize(nsteps);
    lin.P1.resize(nsteps);
    lin.VG0.resize(nsteps);
    lin.VG1.resize(nsteps);

    const MatrixXd no_gamma_chain; /* empty: V_xfoil is the fixed-gamma part */
    Vector2d xyw_i = state.wake.x.col(0);
    for (int i = 0; i < nsteps; ++i) {
        const double ds_i = sv(i + 1) - sv(i);
        lin.ds(i) = ds_i;
        lin.dds_dds1(i) = dsv_dds1(i + 1) - dsv_dds1(i);

        VelDeriv vd0 = InviscidVelocityWithDerivExt(
            state, xyw_i, state.isol.gam, no_gamma_chain, alpha);
        double V0n = vd0.V.norm();
        if (V0n < 1e-12) V0n = 1e-12;
        const Vector2d v1 = vd0.V / V0n;
        lin.v1[i] = v1;
        lin.J1[i] = (Matrix2d::Identity() - v1 * v1.transpose()) / V0n;
        lin.Vxi0[i] = vd0.V_xi;
        lin.P0[i] = std::move(vd0.V_xfoil);
        lin.VG0[i] = std::move(vd0.V_G);

        const Vector2d x_pred = xyw_i + ds_i * v1;
        VelDeriv vd1 = InviscidVelocityWithDerivExt(
            state, x_pred, state.isol.gam, no_gamma_chain, alpha);
        double V1n = vd1.V.norm();
        if (V1n < 1e-12) V1n = 1e-12;
        const Vector2d v2 = vd1.V / V1n;
        lin.v2[i] = v2;
        lin.J2[i] = (Matrix2d::Identity() - v2 * v2.transpose()) / V1n;
        lin.Vxi1[i] = vd1.V_xi;
        lin.P1[i] = std::move(vd1.V_xfoil);
        lin.VG1[i] = std::move(vd1.V_G);

        xyw_i = xyw_i + ds_i * 0.5 * (v1 + v2);
    }
    return lin;
}

void WakeSweepVJP(const WakeSweepLinearization& lin, const MatrixXd& xbar,
                  const MatrixXd& tbar, MatrixXd& foil_bar, MatrixXd& gam_bar) {
    const int N = lin.N;
    const int Nw = lin.Nw;
    const long K = xbar.cols();
    if (xbar.rows() != 2 * Nw || tbar.rows() != 2 * Nw || tbar.cols() != K ||
        foil_bar.rows() != 2 * N || foil_bar.cols() != K ||
        gam_bar.rows() != N + 1 || gam_bar.cols() != K) {
        throw std::invalid_argument("WakeSweepVJP: inconsistent dimensions");
    }

    /* a.block(2*j) holds the adjoint of wake.x[:, j], starting from the
     * direct seeds; later points feed earlier ones as the sweep runs back. */
    MatrixXd a(2 * Nw, K);
    for (int j = 0; j < Nw; ++j) {
        a.row(2 * j) = xbar.row(j);
        a.row(2 * j + 1) = xbar.row(Nw + j);
    }
    auto t_seed = [&](int j) {
        MatrixXd t(2, K);
        t.row(0) = tbar.row(j);
        t.row(1) = tbar.row(Nw + j);
        return t;
    };

    Eigen::RowVectorXd ds1_bar = Eigen::RowVectorXd::Zero(K);
    for (int i = Nw - 2; i >= 0; --i) {
        const double ds_i = lin.ds(i);
        const Vector2d& v1 = lin.v1[i];
        const Vector2d& v2 = lin.v2[i];
        const MatrixXd a_next = a.middleRows(2 * (i + 1), 2);

        /* Corrector: x_{i+1} = x_i + ds_i * 0.5 * (v1 + v2). */
        MatrixXd v1bar = (0.5 * ds_i) * a_next;
        MatrixXd v2bar = (0.5 * ds_i) * a_next;
        Eigen::RowVectorXd dsbar = (0.5 * (v1 + v2)).transpose() * a_next;
        MatrixXd a_i = a.middleRows(2 * i, 2) + a_next;

        /* wake.t[i] = v1 of step i; wake.t[Nw-1] = v2 of the last step. */
        v1bar += t_seed(i);
        if (i == Nw - 2) v2bar += t_seed(Nw - 1);

        /* v2 = V(x_pred)/|V(x_pred)| */
        const MatrixXd V1bar = lin.J2[i].transpose() * v2bar;
        foil_bar.noalias() += lin.P1[i].transpose() * V1bar;
        gam_bar.noalias() += lin.VG1[i].transpose() * V1bar;
        const MatrixXd xpbar = lin.Vxi1[i].transpose() * V1bar;

        /* Predictor: x_pred = x_i + ds_i * v1. */
        a_i += xpbar;
        v1bar += ds_i * xpbar;
        dsbar += v1.transpose() * xpbar;

        /* v1 = V(x_i)/|V(x_i)| */
        const MatrixXd V0bar = lin.J1[i].transpose() * v1bar;
        foil_bar.noalias() += lin.P0[i].transpose() * V0bar;
        gam_bar.noalias() += lin.VG0[i].transpose() * V0bar;
        a_i += lin.Vxi0[i].transpose() * V0bar;

        a.middleRows(2 * i, 2) = a_i;
        /* ds_i depends on foil.x through the first wake spacing ds1. */
        ds1_bar += lin.dds_dds1(i) * dsbar;
    }

    foil_bar.noalias() += lin.D0.transpose() * a.middleRows(0, 2);
    foil_bar.noalias() += lin.dds1_dx * ds1_bar;
}

/* ============================================================================
 * Alpha-direction analytic chains. These mirror the foil.x chains but with
 * scalar derivatives (one column instead of 2N). Used by the analytic
 * DRVDAlpha to satisfy the "no FD" mandate.
 * ============================================================================
 */

/**
 * @brief d(gamma)/d(alpha) [N+1] — closed-form from gamref[:,0]·cos α +
 * gamref[:,1]·sin α.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length N+1.
 */
VectorXd DGammaDAlpha(const FoilState& state) {
    int N = state.foil.N;
    double alpha = state.oper.alpha;
    double ca = Cosd(alpha);
    double sa = Sind(alpha);
    /* gam[i] = gamref[i,0]·cos α + gamref[i,1]·sin α
     * dgam[i]/dα = (-gamref[i,0]·sin α + gamref[i,1]·cos α)·deg2rad
     */
    VectorXd dgam(N + 1);
    for (int i = 0; i < N; ++i) {
        dgam(i) =
            (-state.isol.gamref(i, 0) * sa + state.isol.gamref(i, 1) * ca) *
            kDegToRad;
    }
    dgam(N) = 0.0; /* streamfunction component */
    return dgam;
}

/**
 * @brief InviscidVelocity α-derivative: ∂V/∂α [2] at xi, with the chain through
 * V_G · dgam/dα (gamma part) plus freestream V_∞ ·(-sin α, cos α)·deg2rad.
 *
 * @param state FoilState containing the current solution and geometry.
 * @param xi Evaluation point in 2D space.
 * @param dgam_dalpha d(gamma)/d(alpha) vector of length N+1.
 * @return Vector2d of ∂V/∂α.
 */
Vector2d InviscidVelocityDAlpha(const FoilState& state, const Vector2d& xi,
                                const VectorXd& dgam_dalpha) {
    /* ∂V/∂α = V_G · dgam/dα + V_∞ · (-sin α, cos α) · deg2rad.
     * V_G is the per-gamma sensitivity returned by InviscidVelocity with
     * gradient=true. Use the existing InviscidVelocityWithDerivExt with
     * vortex strengths = dgam_dalpha and freestream contribution added.
     */
    int N = state.foil.N;
    double alpha = state.oper.alpha;
    double Vinf = state.oper.Vinf;
    /* Use InviscidVelocityWithDerivExt with vortex strengths = dgam_dalpha.
     * It returns V including freestream. Subtract Vinf and add the
     * explicit α-derivative of the freestream.
     */
    MatrixXd dgam_dummy = MatrixXd::Zero(N + 1, 1); /* we use the vortex slot */
    VelDeriv vd =
        InviscidVelocityWithDerivExt(state, xi, dgam_dalpha, dgam_dummy, alpha);
    Vector2d Vinf_vec(Vinf * Cosd(alpha), Vinf * Sind(alpha));
    Vector2d V_G_dgam = vd.V - Vinf_vec; /* gamma-induced part only */
    /* Freestream α-derivative */
    Vector2d dVinf_dalpha(-Vinf * Sind(alpha) * kDegToRad,
                          Vinf * Cosd(alpha) * kDegToRad);
    return V_G_dgam + dVinf_dalpha;
}

/**
 * @brief d(wake.x)/d(alpha) [2 x Nw] — mirrors DWakeXDFoilX's
 * predictor-corrector sweep with dgam_dα (length N+1) and direct freestream
 * alpha contribution in the velocity.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2 x Nw].
 */
MatrixXd DWakeXDAlpha(const FoilState& state) {
    if (std::isfinite(state.wake.tracing_alpha) &&
        state.oper.alpha != state.wake.tracing_alpha) {
        FoilState traced = state;
        traced.oper.alpha = state.wake.tracing_alpha;
        traced.isol.gam = traced.isol.gamref.col(0) * Cosd(traced.oper.alpha) +
                          traced.isol.gamref.col(1) * Sind(traced.oper.alpha);
        return DWakeXDAlpha(traced);
    }

    /* Predictor-corrector chain in scalar α direction. Mirrors DWakeXDFoilX
     * but with dgam_dα (length N+1) instead of dgam_dx (matrix), and with
     * direct α-dependence in freestream V handled via InviscidVelocityDAlpha.
     */
    int N = state.foil.N;
    int Nw = state.wake.N;
    double alpha = state.oper.alpha;
    double chord = state.geom.chord;
    double wakelen = state.geom.wakelen;

    VectorXd dgam_dalpha = DGammaDAlpha(state);

    /* foil.s doesn't depend on α, so ds1 has no α-derivative, hence sv too. */
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    /* Initial wake.x[:, 0] = x_te + 1e-5·chord·t_te. No α dependence. */
    Vector2d dxyw_i = Vector2d::Zero();

    MatrixXd dwake_x = MatrixXd::Zero(2, Nw);
    Vector2d xyw_i = state.wake.x.col(0);

    for (int i = 0; i < Nw - 1; ++i) {
        double ds_i = sv(i + 1) - sv(i);
        /* dsv/dα = 0 (ds1 doesn't depend on α) */

        /* V0 at xyw_i. ∂V/∂α = V_G · dgam_dalpha + freestream_α.
         * Also need ∂V/∂xi for the chain through dxyw_i.
         * Get V_xi from InviscidVelocityWithDeriv (any direction gives same
         * V_xi).
         */
        VelDeriv vd_for_xi = InviscidVelocityWithDeriv(
            state, xyw_i, MatrixXd::Zero(N + 1, 1), alpha);
        Vector2d V0 = vd_for_xi.V;
        Matrix2d V0_xi = vd_for_xi.V_xi;
        Vector2d V0_alpha = InviscidVelocityDAlpha(state, xyw_i, dgam_dalpha);

        double V0n = V0.norm();
        if (V0n < 1e-12) V0n = 1e-12;
        Vector2d v1 = V0 / V0n;
        Matrix2d norm_jac1 = (Matrix2d::Identity() - v1 * v1.transpose()) / V0n;
        Vector2d dV0 = V0_alpha + V0_xi * dxyw_i;
        Vector2d dv1 = norm_jac1 * dV0;

        /* Predictor step */
        Vector2d x_pred = xyw_i + ds_i * v1;
        Vector2d dx_pred = dxyw_i + ds_i * dv1; /* ds_i not α-dependent */

        /* V1 at x_pred */
        VelDeriv vd1_xi = InviscidVelocityWithDeriv(
            state, x_pred, MatrixXd::Zero(N + 1, 1), alpha);
        Vector2d V1 = vd1_xi.V;
        Matrix2d V1_xi = vd1_xi.V_xi;
        Vector2d V1_alpha = InviscidVelocityDAlpha(state, x_pred, dgam_dalpha);

        double V1n = V1.norm();
        if (V1n < 1e-12) V1n = 1e-12;
        Vector2d v2 = V1 / V1n;
        Matrix2d norm_jac2 = (Matrix2d::Identity() - v2 * v2.transpose()) / V1n;
        Vector2d dV1 = V1_alpha + V1_xi * dx_pred;
        Vector2d dv2 = norm_jac2 * dV1;

        /* Corrector */
        Vector2d xyw_next = xyw_i + ds_i * 0.5 * (v1 + v2);
        Vector2d dxyw_next = dxyw_i + ds_i * 0.5 * (dv1 + dv2);

        dwake_x.col(i + 1) = dxyw_next;

        xyw_i = xyw_next;
        dxyw_i = dxyw_next;
    }

    return dwake_x;
}

/**
 * @brief d(wake.t)/d(alpha) [2 x Nw] — same iterative sweep as DWakeXDAlpha but
 * stores dv1 (for j<Nw-1) and dv2 (for j=Nw-1) as the wake.t α-derivative.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return MatrixXd of shape [2 x Nw].
 */
MatrixXd DWakeTDAlpha(const FoilState& state) {
    if (std::isfinite(state.wake.tracing_alpha) &&
        state.oper.alpha != state.wake.tracing_alpha) {
        FoilState traced = state;
        traced.oper.alpha = state.wake.tracing_alpha;
        traced.isol.gam = traced.isol.gamref.col(0) * Cosd(traced.oper.alpha) +
                          traced.isol.gamref.col(1) * Sind(traced.oper.alpha);
        return DWakeTDAlpha(traced);
    }

    /* Same iterative sweep as DWakeXDAlpha but stores dv1 (for j<Nw-1)
     * and dv2 (for j=Nw-1) as the wake.t α-derivative.
     */
    int N = state.foil.N;
    int Nw = state.wake.N;
    double alpha = state.oper.alpha;
    double chord = state.geom.chord;
    double wakelen = state.geom.wakelen;

    VectorXd dgam_dalpha = DGammaDAlpha(state);

    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    Vector2d dxyw_i = Vector2d::Zero();
    MatrixXd dwake_t = MatrixXd::Zero(2, Nw);
    Vector2d xyw_i = state.wake.x.col(0);

    for (int i = 0; i < Nw - 1; ++i) {
        double ds_i = sv(i + 1) - sv(i);

        VelDeriv vd_xi = InviscidVelocityWithDeriv(
            state, xyw_i, MatrixXd::Zero(N + 1, 1), alpha);
        Vector2d V0 = vd_xi.V;
        Vector2d V0_alpha = InviscidVelocityDAlpha(state, xyw_i, dgam_dalpha);
        double V0n = V0.norm();
        if (V0n < 1e-12) V0n = 1e-12;
        Vector2d v1 = V0 / V0n;
        Matrix2d norm_jac1 = (Matrix2d::Identity() - v1 * v1.transpose()) / V0n;
        Vector2d dV0 = V0_alpha + vd_xi.V_xi * dxyw_i;
        Vector2d dv1 = norm_jac1 * dV0;
        dwake_t.col(i) = dv1;

        Vector2d x_pred = xyw_i + ds_i * v1;
        Vector2d dx_pred = dxyw_i + ds_i * dv1;

        VelDeriv vd1_xi = InviscidVelocityWithDeriv(
            state, x_pred, MatrixXd::Zero(N + 1, 1), alpha);
        Vector2d V1 = vd1_xi.V;
        Vector2d V1_alpha = InviscidVelocityDAlpha(state, x_pred, dgam_dalpha);
        double V1n = V1.norm();
        if (V1n < 1e-12) V1n = 1e-12;
        Vector2d v2 = V1 / V1n;
        Matrix2d norm_jac2 = (Matrix2d::Identity() - v2 * v2.transpose()) / V1n;
        Vector2d dV1 = V1_alpha + vd1_xi.V_xi * dx_pred;
        Vector2d dv2 = norm_jac2 * dV1;
        dwake_t.col(i + 1) = dv2; /* overwritten on next iter unless last */

        Vector2d xyw_next = xyw_i + ds_i * 0.5 * (v1 + v2);
        Vector2d dxyw_next = dxyw_i + ds_i * 0.5 * (dv1 + dv2);
        xyw_i = xyw_next;
        dxyw_i = dxyw_next;
    }
    return dwake_t;
}

/**
 * @brief d(wake.s)/d(alpha) [Nw] — cumulative panel norms along wake.x
 * perturbation.
 *
 * wake.s[i] depends on α via:
 *   wake.s[i] = foil.s[N-1] + sv[i] where sv depends on ds1 (no α).
 * BUT BuildWake actually re-uses foil.s as a base; wake.s itself is set
 * from sv. Per BuildWake the cumulative panel norms aren't used — wake.s
 * is sw = foil.s[N-1] + sv. So dwake.s/dα = 0!
 *
 * (The cumulative-norm reinterpretation is done by recompute_geom in the
 * CS extract script for the FOIL only. wake.s is taken directly from sv.)
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of zeros, length Nw.
 */
VectorXd DWakeSDAlpha(const FoilState& state) {
    /* wake.s[i] depends on α via:
     *   wake.s[i] = foil.s[N-1] + sv[i] where sv depends on ds1 (no α).
     * BUT BuildWake actually re-uses foil.s as a base; wake.s itself is set
     * from sv. Per BuildWake the cumulative panel norms aren't used — wake.s
     * is sw = foil.s[N-1] + sv. So dwake.s/dα = 0!
     *
     * (The cumulative-norm reinterpretation is done by recompute_geom in the
     * CS extract script for the FOIL only. wake.s is taken directly from sv.)
     */
    int Nw = state.wake.N;
    return VectorXd::Zero(Nw);
}

/**
 * @brief d(sstag)/d(alpha) — scalar. Chain through dgamma at the stagnation
 * panel.
 *
 * Production state is ue-based xi (after the post-Newton-swap removal in
 * commit 4051c2d). sstag = w1·foil.s[Istag[0]] + w2·foil.s[Istag[1]]
 * with weights depending on ue[Istag[0]] and ue[Istag[1]] (via
 * StagPointMove). At fixed U (which is what ∂R/∂α partial requires),
 * ue at Istag is fixed → sstag is fixed → ∂sstag/∂α = 0.
 *
 * The previous implementation used sstag_g·dgam/dα (gamma-based chain),
 * valid only when production xi was gamma-based.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return 0.0 (sstag is α-independent at fixed U).
 */
double DSstagDAlpha(const FoilState& state) {
    /* Production state is ue-based xi (after the post-Newton-swap removal in
     * commit 4051c2d). sstag = w1·foil.s[Istag[0]] + w2·foil.s[Istag[1]]
     * with weights depending on ue[Istag[0]] and ue[Istag[1]] (via
     * StagPointMove). At fixed U (which is what ∂R/∂α partial requires),
     * ue at Istag is fixed → sstag is fixed → ∂sstag/∂α = 0.
     *
     * The previous implementation used sstag_g·dgam/dα (gamma-based chain),
     * valid only when production xi was gamma-based.
     */
    (void)state;
    return 0.0;
}

/**
 * @brief d(xi)/d(alpha) [Nsys] — chain through sstag (foil-side) and wake.s
 * (wake-side).
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of length Nsys.
 */
VectorXd DXiDAlpha(const FoilState& state) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    double dsstag = DSstagDAlpha(state);
    VectorXd dwake_s = DWakeSDAlpha(state); /* all zero, but keep for clarity */
    VectorXd dxi(Nsys);
    /* foil-side: xi[i] = |foil.s[i] - sstag| → sign(foil.s[i]-sstag) ·
     *   (d(foil.s)/dα - dsstag/dα) = -sign · dsstag/dα (foil.s α-indep)
     */
    for (int i = 0; i < N; ++i) {
        double diff = state.foil.s(i) - state.isol.sstag;
        double sgn = (diff >= 0) ? 1.0 : -1.0;
        dxi(i) = sgn * (0.0 - dsstag);
    }
    /* wake-side: xi[N+j] = wake.s[j] - sstag → dwake.s[j]/dα - dsstag/dα */
    for (int j = 0; j < Nw; ++j) {
        dxi(N + j) = dwake_s(j) - dsstag;
    }
    return dxi;
}

/**
 * @brief d(wgap)/d(alpha) [Nw] — chain through wake.s and xi (wgap's xib).
 *
 * wgap depends on hTE (no α), dtdx (no α), foil-TE coords (no α),
 * and xi[N+i] - xi[N] (wake-only arc-length from TE).
 * wake.s is α-independent (DWakeSDAlpha returns 0), and the wake-portion
 * of xi simplifies to wake.s - sstag (with sstag α-dep canceling out
 * between consecutive wake stations) → dwgap/dα = 0.
 *
 * @param state FoilState containing the current solution and geometry.
 * @return VectorXd of zeros, length Nw.
 */
VectorXd DWgapDAlpha(const FoilState& state) {
    /* wgap depends on hTE (no α), dtdx (no α), foil-TE coords (no α),
     * and xi[N+i] - xi[N] (wake-only arc-length from TE).
     * wake.s is α-independent (DWakeSDAlpha returns 0), and the wake-portion
     * of xi simplifies to wake.s - sstag (with sstag α-dep canceling out
     * between consecutive wake stations) → dwgap/dα = 0.
     */
    int Nw = state.wake.N;
    return VectorXd::Zero(Nw);
}

} /* namespace scirocco */
