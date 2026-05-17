/*
    cppfoil: C++ Airfoil Analysis Library

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
#include "cppfoil/core/viscous_geom_deriv.h"

#include <cmath>
#include <stdexcept>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid_deriv.h"

namespace cppfoil {

// ============================================================================
// d(foil.s)/d(foil.x) for the cumulative-panel-norm definition of foil.s
// (i.e. matches cppfoil's RecomputeArcLength).
//
// foil.s(0) = 0
// foil.s(i) = foil.s(i-1) + |x[:,i] - x[:,i-1]|   for i >= 1
//
// Let p(i) = |x[:,i] - x[:,i-1]| be the i-th panel length. Then
//   d(p(i))/d(x[d,k]) = (x[d,i] - x[d,i-1])/p(i)  if k=i
//                     = -(x[d,i] - x[d,i-1])/p(i) if k=i-1
//                     = 0                          otherwise
// and d(foil.s(i))/d(x[d,k]) = sum_{j<=i} d(p(j))/d(x[d,k]).
// ============================================================================
MatrixXd DFoilSDFoilXCumulative(const Panel& foil) {
    int N = foil.N;
    MatrixXd out = MatrixXd::Zero(N, 2 * N);
    // Compute panel-length unit vectors
    // dp_j/dx[d,j] = u[d,j],  dp_j/dx[d,j-1] = -u[d,j]
    // where u[:,j] = (x[:,j] - x[:,j-1]) / |x[:,j] - x[:,j-1]|
    // Then foil.s(i) = sum_{j=1..i} p(j), so
    //   d(foil.s(i))/d(x[d,k])
    //     = u[d,k]   if 1 <= k <= i        (k contributes via p(k) and p(k+1))
    //     wait — let's be careful:
    //   contribution to d(foil.s(i))/d(x[d,k]) from p(j):
    //     = +u[d,j] if k==j and j<=i
    //     = -u[d,j] if k==j-1 and j<=i, i.e. k<=i-1 and j=k+1
    //
    // So for foil.s(i):
    //   if k==0: only p(1) contributes (-u[d,1])
    //   if 1<=k<=i-1: p(k) contributes (+u[d,k]) and p(k+1) contributes (-u[d,k+1])
    //   if k==i: p(i) contributes (+u[d,i])
    //   if k>i: zero
    // ----------------------------------------------------------------------
    // Compute u[:,j] (the unit tangent of panel j, j ∈ [1,N-1])
    MatrixXd u(2, N);  // u(:,0) unused
    u.setZero();
    for (int j = 1; j < N; ++j) {
        Vector2d dxv = foil.x.col(j) - foil.x.col(j - 1);
        double pj = dxv.norm();
        if (pj < 1e-14) pj = 1e-14;
        u.col(j) = dxv / pj;
    }

    for (int i = 1; i < N; ++i) {
        for (int j = 1; j <= i; ++j) {
            // p(j) contributes +u[d,j] to d(foil.s(i))/d(x[d,j])
            // and -u[d,j] to d(foil.s(i))/d(x[d,j-1])
            for (int d = 0; d < 2; ++d) {
                out(i, 2 * j + d) += u(d, j);
                out(i, 2 * (j - 1) + d) -= u(d, j);
            }
        }
    }
    // row 0 (foil.s(0)=0) is all zeros by construction
    return out;
}

// ============================================================================
// GetTEInfoDeriv — derivatives of TE-info quantities w.r.t. foil.x[k]_dim
//
// Returns nonzero only for k ∈ {0, 1, N-2, N-1}. Computes simultaneously:
//   dtcp, dtdp, dhTE, ddtdx
//
// Definitions (matching geometry.cpp:170-210):
//   t1 = normalize(x[:,0] - x[:,1])
//   t2 = normalize(x[:,N-1] - x[:,N-2])
//   tbar = 0.5*(t1+t2);   t = normalize(tbar)
//   s = x[:,N-1] - x[:,0];   p = s/|s|
//   hTE = -s(0)*t(1) + s(1)*t(0)
//   dtdx = t1(0)*t2(1) - t2(0)*t1(1)
//   tcp = |t(0)*p(1) - t(1)*p(0)|
//   tdp = t(0)*p(0) + t(1)*p(1)
// ============================================================================
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

    // d(v1)/d(x[k]_dim), d(v2)/d(x[k]_dim), d(s)/d(x[k]_dim)
    Vector2d dv1 = Vector2d::Zero();
    Vector2d dv2 = Vector2d::Zero();
    Vector2d ds  = Vector2d::Zero();
    if (k == 0)         dv1(dim) =  1.0;
    else if (k == 1)    dv1(dim) = -1.0;
    if (k == N - 1)     dv2(dim) =  1.0;
    else if (k == N - 2) dv2(dim) = -1.0;
    if (k == N - 1)     ds(dim)  =  1.0;
    else if (k == 0)    ds(dim)  = -1.0;

    auto normalize_jac = [](const Vector2d& u, double un) -> Matrix2d {
        return (Matrix2d::Identity() - u * u.transpose()) / un;
    };

    Vector2d dt1 = normalize_jac(t1, v1n) * dv1;
    Vector2d dt2 = normalize_jac(t2, v2n) * dv2;
    Vector2d dtbar = 0.5 * (dt1 + dt2);
    Vector2d dt = normalize_jac(t, tbn) * dtbar;
    Vector2d dp = normalize_jac(p, sn) * ds;

    // dhTE = -ds(0)*t(1) - s(0)*dt(1) + ds(1)*t(0) + s(1)*dt(0)
    dhTE = -ds(0) * t(1) - s(0) * dt(1) + ds(1) * t(0) + s(1) * dt(0);

    // ddtdx = dt1(0)*t2(1) + t1(0)*dt2(1) - dt2(0)*t1(1) - t2(0)*dt1(1)
    ddtdx = dt1(0) * t2(1) + t1(0) * dt2(1) - dt2(0) * t1(1) - t2(0) * dt1(1);

    // tcp = |t(0)*p(1) - t(1)*p(0)|
    double cross = t(0) * p(1) - t(1) * p(0);
    double sgn = (cross >= 0) ? 1.0 : -1.0;
    double dcross = dt(0) * p(1) + t(0) * dp(1) - dt(1) * p(0) - t(1) * dp(0);
    dtcp = sgn * dcross;

    // tdp = t·p
    dtdp = dt.dot(p) + t.dot(dp);
}

// ============================================================================
// Helper: forward sensitivity of inviscid gamma at one node via reverse mode
//
// d(gamma_i)/d(foil.x) = -(AIC^T)^{-1} e_i^T · dR_AIC/d(foil.x)
//                      = -AICResidualSensVJP(state, λ_i, alpha),
// where λ_i = (AIC^T)^{-1} e_i (e_i ∈ R^{N+1}, the i-th basis vector).
// ============================================================================
static VectorXd DGammaAtNodeDFoilX(const FoilState& state, int gamma_node,
                                   double alpha) {
    const MatrixXd& AIC = state.isol.AIC;
    int dim = static_cast<int>(AIC.rows());
    if (gamma_node < 0 || gamma_node >= dim) {
        throw std::out_of_range("DGammaAtNodeDFoilX: gamma_node out of range");
    }
    VectorXd e = VectorXd::Zero(dim);
    e(gamma_node) = 1.0;
    MatrixXd AIC_T = AIC.transpose();
    auto lu = AIC_T.partialPivLu();
    VectorXd lambda = lu.solve(e);
    VectorXd vjp = AICResidualSensVJP(state, lambda, alpha);
    return -vjp;
}

// ============================================================================
// DSstagDFoilX
//   sstag = w1 * S[i0] + w2 * S[i1],
//   where i0=Istag[0], i1=Istag[1], g0=gam(i0), g1=gam(i1),
//         den = g1 - g0,  w1 = g1/den,  w2 = -g0/den.
//
//   ∂sstag/∂S[i0] = w1
//   ∂sstag/∂S[i1] = w2
//   ∂sstag/∂g0    = sstag_g(0)   (already cached)
//   ∂sstag/∂g1    = sstag_g(1)
// ============================================================================
VectorXd DSstagDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int twoN = 2 * N;
    int i0 = state.isol.Istag[0];
    int i1 = state.isol.Istag[1];
    double g0 = state.isol.gam(i0);
    double g1 = state.isol.gam(i1);
    double den = g1 - g0;
    double w1 = g1 / den;
    double w2 = -g0 / den;
    double dsstag_dg0 = state.isol.sstag_g(0);
    double dsstag_dg1 = state.isol.sstag_g(1);
    double alpha = state.oper.alpha;

    // d(g0)/d(foil.x), d(g1)/d(foil.x)
    VectorXd dg0 = DGammaAtNodeDFoilX(state, i0, alpha);
    VectorXd dg1 = DGammaAtNodeDFoilX(state, i1, alpha);

    VectorXd dsstag = VectorXd::Zero(twoN);
    // Chain through foil.s
    dsstag += w1 * dfoil_s_dx.row(i0).transpose();
    dsstag += w2 * dfoil_s_dx.row(i1).transpose();
    // Chain through gamma
    dsstag += dsstag_dg0 * dg0;
    dsstag += dsstag_dg1 * dg1;
    return dsstag;
}

// ============================================================================
// DWakeSDFoilX
//
// wake.s[i] = foil.s[N-1] + sv[i]
//   sv[0] = 0
//   sv[i>0] from SpaceGeom(ds1, wakelen*chord, Nw) — depends on
//     ds1 = 0.5 * ((foil.s[1] - foil.s[0]) + (foil.s[N-1] - foil.s[N-2]))
//     and on chord (assumed constant in foil.x; if chord changes, add chain).
//
// Therefore d(wake.s[i])/d(foil.x) =
//   d(foil.s[N-1])/d(foil.x) + (d(sv[i])/d(ds1)) * d(ds1)/d(foil.x)
//
// We obtain d(sv[i])/d(ds1) numerically via one extra SpaceGeom-style call
// (analytic derivative through the Newton iteration is straightforward —
// see below).
//
// Geometric spacing: sv[i] = sum_{j=1..i} ds_{j-1} * r^{j-1}, where r solves
//   L = ds1 * (r^N - 1)/(r - 1),    N = Np - 1
// Implicit differentiation:
//   0 = (r^N - 1)/(r-1) * ddS + ds1 * [N*r^(N-1)/(r-1) - (r^N-1)/(r-1)^2] dr
// => dr/dds1 = -((r^N-1)/(r-1)) / { ds1 * [N*r^(N-1)/(r-1) - (r^N-1)/(r-1)^2] }
// Then dsv[i]/dds1 = sum_{j=0..i-1} (r^j + ds1 * j*r^(j-1) * dr/dds1).
// ============================================================================
namespace {
// Compute (sv, dsv_dds1) given Np, ds1, L (matching wake.cpp SpaceGeom).
void SpaceGeomWithDeriv(double ds1, double L, int Np, VectorXd& sv,
                        VectorXd& dsv_dds1) {
    if (Np <= 1) {
        sv = VectorXd::Zero(1);
        dsv_dds1 = VectorXd::Zero(1);
        return;
    }
    int N = Np - 1;
    double d = L / ds1;

    // Initial guess for r (same as wake.cpp)
    double a = N * (N - 1.0) * (N - 2.0) / 6.0;
    double b = N * (N - 1.0) / 2.0;
    double c = N - d;
    double disc = std::max(b * b - 4.0 * a * c, 0.0);
    double r = 1.0 + (-b + std::sqrt(disc)) / (2.0 * a);
    for (int k = 0; k < 10; ++k) {
        double rN = std::pow(r, N);
        double R = rN - 1.0 - d * (r - 1.0);
        double R_r = N * std::pow(r, N - 1) - d;
        double dr = -R / R_r;
        if (std::abs(dr) < 1e-12) break;
        r += dr;
    }

    // dr/dds1: chord and L are constants → d/dds1 of  rN - 1 - (L/ds1)*(r-1) = 0
    //   N r^{N-1} dr + L/ds1^2 * (r-1) - L/ds1 * dr = 0
    //   dr * (N r^{N-1} - L/ds1) = -L * (r-1) / ds1^2
    //   dr_dds1 = -L * (r-1) / (ds1^2 * (N r^{N-1} - L/ds1))
    double dr_dds1 = -L * (r - 1.0) /
                     (ds1 * ds1 * (N * std::pow(r, N - 1) - L / ds1));

    sv = VectorXd::Zero(Np);
    dsv_dds1 = VectorXd::Zero(Np);
    double dx = ds1;
    double ddx_dds1 = 1.0;
    for (int i = 1; i < Np; ++i) {
        sv(i) = sv(i - 1) + dx;
        dsv_dds1(i) = dsv_dds1(i - 1) + ddx_dds1;
        // Next-interval update: dx <- dx * r
        ddx_dds1 = ddx_dds1 * r + dx * dr_dds1;
        dx *= r;
    }
}
}  // namespace

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

    // d(ds1)/d(foil.x) = 0.5 * ( d(foil.s[1])/d(foil.x) - d(foil.s[0])/d(foil.x)
    //                          + d(foil.s[N-1])/d(foil.x)
    //                          - d(foil.s[N-2])/d(foil.x) )
    VectorXd dds1 = 0.5 * (dfoil_s_dx.row(1).transpose() -
                           dfoil_s_dx.row(0).transpose() +
                           dfoil_s_dx.row(N - 1).transpose() -
                           dfoil_s_dx.row(N - 2).transpose());

    // d(foil.s[N-1])/d(foil.x)
    VectorXd dSTE = dfoil_s_dx.row(N - 1).transpose();

    MatrixXd dwake_s = MatrixXd::Zero(Nw, twoN);
    for (int i = 0; i < Nw; ++i) {
        // d(wake.s[i])/d(foil.x) = dSTE + dsv_dds1[i] * dds1
        dwake_s.row(i) = (dSTE + dsv_dds1(i) * dds1).transpose();
    }
    return dwake_s;
}

// ============================================================================
// DXiDFoilX
//   xi[i]    = |foil.s[i] - sstag|      for i ∈ [0, N)
//   xi[N+j]  = wake.s[j] - sstag        for j ∈ [0, Nw)
// ============================================================================
MatrixXd DXiDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    int twoN = 2 * N;

    VectorXd dsstag = DSstagDFoilX(state, dfoil_s_dx);  // [2N]
    MatrixXd dwake_s = DWakeSDFoilX(state, dfoil_s_dx); // [Nw x 2N]

    MatrixXd dxi = MatrixXd::Zero(Nsys, twoN);
    double sstag = state.isol.sstag;

    for (int i = 0; i < N; ++i) {
        double diff = state.foil.s(i) - sstag;
        double sgn = (diff >= 0.0) ? 1.0 : -1.0;
        // d(xi[i])/d(foil.x) = sgn * (d(foil.s[i])/d(foil.x) - d(sstag)/d(foil.x))
        dxi.row(i) =
            sgn * (dfoil_s_dx.row(i) - dsstag.transpose());
    }
    for (int j = 0; j < Nw; ++j) {
        dxi.row(N + j) = dwake_s.row(j) - dsstag.transpose();
    }
    return dxi;
}

// ============================================================================
// DWgapDFoilX
//   wgap[j] depends on hTE, dtdx, xib = (xi[N+j] - xi[N])/Lw,  Lw = flen*hTE
//
//   ξib = (xi[N+j] - xi[N]) / Lw
//       = ((wake.s[j] - sstag) - (foil.s[N-1] - sstag)) / Lw   if foil.s[N-1] >= sstag
//
//   But wait — xi[N]  = |foil.s[N-1] - sstag|. Assuming sstag < foil.s[N-1]
//   on the upper-surface side, xi[N] = foil.s[N-1] - sstag, and
//   xi[N+j] - xi[N] = wake.s[j] - foil.s[N-1] = sv[j]. So sstag drops out
//   and the chain is purely through ds1 and chord.
//
//   wgap[j] = hTE * (1 + (2 + flen*dtdx)*xib) * (1-xib)^2   if xib<=1 else 0
//
//   Let A = 1 + (2 + flen*dtdx)*xib,  B = (1-xib)^2
//   wgap = hTE * A * B
//   d(wgap)/d(.) = dhTE * A*B
//                + hTE * dA * B
//                + hTE * A * dB
//   dA = flen*ddtdx * xib + (2 + flen*dtdx) * dxib
//   dB = -2*(1-xib) * dxib
//   dxib = ( d(xi[N+j]-xi[N])/d(.) * Lw - (xi[N+j]-xi[N]) * dLw ) / Lw^2
//        = d(sv[j])/d(.)/Lw - xib * dLw/Lw
//   dLw = flen * dhTE
// ============================================================================
MatrixXd DWgapDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int twoN = 2 * N;
    double wakelen = state.geom.wakelen;
    double chord = state.geom.chord;
    const double flen = 2.5;

    // Recompute TE info (matches SetWakeGap)
    TEInfo te = GetTEInfo(state.foil.x);
    double hTE = te.hTE;
    double dtdx_raw = te.dtdx;
    double clip_lo = -3.0 / flen;
    double clip_hi = 3.0 / flen;
    bool clipped = (dtdx_raw < clip_lo) || (dtdx_raw > clip_hi);
    double dtdx = std::min(std::max(dtdx_raw, clip_lo), clip_hi);
    double Lw = std::max(flen * hTE, 1e-10);

    // d(hTE), d(dtdx_raw) for each (k, dim) over TE-affecting nodes
    // We'll need 4 nodes × 2 dims = 8 entries → store sparsely.
    auto te_deriv_col = [&](int k, int dim, double& dhTE, double& ddtdx) {
        double dtcp, dtdp, dhTE_loc, ddtdx_loc;
        GetTEInfoDeriv(state.foil.x, k, dim, dtcp, dtdp, dhTE_loc, ddtdx_loc);
        dhTE = dhTE_loc;
        ddtdx = clipped ? 0.0 : ddtdx_loc;
    };

    // d(sv[j])/d(foil.x) via SpaceGeom analytical sensitivity
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));
    double L = wakelen * chord;
    VectorXd sv, dsv_dds1;
    SpaceGeomWithDeriv(ds1, L, Nw, sv, dsv_dds1);

    VectorXd dds1 = 0.5 * (dfoil_s_dx.row(1).transpose() -
                           dfoil_s_dx.row(0).transpose() +
                           dfoil_s_dx.row(N - 1).transpose() -
                           dfoil_s_dx.row(N - 2).transpose());

    MatrixXd dwgap = MatrixXd::Zero(Nw, twoN);

    for (int j = 0; j < Nw; ++j) {
        double xib = sv(j) / Lw;
        if (xib > 1.0) {
            // wgap[j] = 0  outside the support; derivative is 0 too.
            continue;
        }
        double A = 1.0 + (2.0 + flen * dtdx) * xib;
        double B = (1.0 - xib) * (1.0 - xib);

        // For each (k, dim), assemble d(wgap[j])/d(x[k]_dim)
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
                // xib = sv/Lw → dxib = dsv/Lw - sv/Lw^2 * dLw = dsv/Lw - xib*dLw/Lw
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

// ============================================================================
// DWakeXDFoilX — STUB
//
// The wake-coordinate chain requires a forward sweep through BuildWake's
// predictor-corrector iteration. Each step is:
//
//   V0  = InviscidVelocity(foil.x, gam, Vinf, alpha, xyw[:,i])
//   v1  = V0 / |V0|
//   x_pred = xyw[:,i] + ds_i * v1
//   V1  = InviscidVelocity(foil.x, gam, Vinf, alpha, x_pred)
//   v2  = V1 / |V1|
//   xyw[:,i+1] = xyw[:,i] + ds_i * 0.5 * (v1 + v2)
//
// Each iteration mixes:
//   - dxyw[:,i]/d(foil.x)         (carries forward)
//   - dInviscidVelocity/d(foil.x) at fixed eval point (panel influence chain)
//   - dInviscidVelocity/d(eval_point) at fixed foil (the eval-point chain)
//   - d(gamma)/d(foil.x) (via reverse-mode AIC)
//   - dds_i/d(foil.x) through dsv/dds1 (from DWakeSDFoilX helper)
//
// This is the chain that caused the prior ~0.216 adjoint error at TE.
// Will be implemented after the closed-form chains pass CS validation.
// ============================================================================
MatrixXd DWakeXDFoilX(const FoilState& state) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    return MatrixXd::Zero(2 * Nw, 2 * N);
}

}  // namespace cppfoil
