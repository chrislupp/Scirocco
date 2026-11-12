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
#include "scirocco/core/captured_transition.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "scirocco/core/boundary_layer.h"
#include "scirocco/core/geometry.h"
#include "scirocco/core/inviscid.h"
#include "scirocco/core/wake.h"
#include "scirocco/properties/bl_properties.h"

namespace scirocco {

namespace {

using Vector5d = Eigen::Matrix<double, 5, 1>;
using Vector10d = Eigen::Matrix<double, 10, 1>;
using Matrix52d = Eigen::Matrix<double, 5, 2>;

/* State rows: glob.U holds [th, ds, ct, ue]; n is appended as row 4. */
constexpr int TH = 0, DS = 1, CT = 2, UE = 3, AMP = 4;
constexpr std::array<int, 4> kLam = {TH, DS, AMP, UE};
constexpr std::array<int, 4> kTurb = {TH, DS, CT, UE};

Vector4d View(const Vector5d& U, const std::array<int, 4>& rows) {
    return {U(rows[0]), U(rows[1]), U(rows[2]), U(rows[3])};
}

Vector5d E5(const Vector4d& v4, const std::array<int, 4>& rows) {
    Vector5d out = Vector5d::Zero();
    for (int j = 0; j < 4; ++j) out(rows[j]) += v4(j);
    return out;
}

Vector10d E10(const Vector8d& v8, const std::array<int, 4>& rows) {
    Vector10d out = Vector10d::Zero();
    for (int j = 0; j < 4; ++j) {
        out(rows[j]) += v8(j);
        out(5 + rows[j]) += v8(4 + j);
    }
    return out;
}

Vector10d Two(const Vector5d& a, const Vector5d& b) {
    Vector10d out;
    out << a, b;
    return out;
}

Vector5d Unit(int i) {
    Vector5d e = Vector5d::Zero();
    e(i) = 1.0;
    return e;
}

/* S((n - center) / w) and d/dn */
std::pair<double, double> Switch(double n, double center, double w) {
    double t = std::tanh((n - center) / w);
    return {0.5 * (1.0 + t), 0.5 * (1.0 - t * t) / w};
}

double ClosureCenter(const Param& p) {
    return p.ncrit + p.capt_b * p.capt_ws + 2.3 * p.capt_wc;
}

std::pair<double, double> Gamma(double n, const Param& p) {
    return Switch(n, ClosureCenter(p), p.capt_wc);
}
std::pair<double, double> Growth(double n, const Param& p) {
    return Switch(n, p.ncrit, p.capt_ws);
}
std::pair<double, double> Release(double n, const Param& p) {
    return Switch(n, p.ncrit - 2.5 * p.capt_ws, p.capt_ws);
}
double RateC(const Param& p) { return 2.2 * p.capt_wc / p.capt_L; }

/* Upwind helper for 5-state linearizations (upw_U has 10 entries). */
std::pair<double, Vector10d> Up(double upw, const Vector10d& upw_U, double f1,
                                const Vector5d& f1_U, double f2,
                                const Vector5d& f2_U) {
    Vector10d f_U = -upw_U * f1 + upw_U * f2;
    f_U.head<5>() += (1.0 - upw) * f1_U;
    f_U.tail<5>() += upw * f2_U;
    return {(1.0 - upw) * f1 + upw * f2, f_U};
}

/* Hs, cf or cDi blended at one station: (1 - g) laminar + g turbulent. */
template <typename Fn>
std::pair<double, Vector5d> Blended(Fn fn, const Vector5d& U,
                                    const Param& param) {
    Param pt = param;
    pt.turb = true;
    auto [ft, ft_U] = fn(View(U, kTurb), pt);
    if (param.wake) return {ft, E5(ft_U, kTurb)};
    Param pl = param;
    pl.turb = false;
    auto [fl, fl_U] = fn(View(U, kLam), pl);
    auto [g, g_n] = Gamma(U(AMP), param);
    Vector5d f_U = (1.0 - g) * E5(fl_U, kLam) + g * E5(ft_U, kTurb);
    f_U(AMP) += (ft - fl) * g_n;
    return {(1.0 - g) * fl + g * ft, f_U};
}

struct XtLin {
    double v;
    Vector5d v_U;
    double v_x;
};

/* f x / th with a blended closure f. */
template <typename Fn>
XtLin BlendedXt(Fn fn, const Vector5d& U, double x, const Param& param) {
    auto [f, f_U] = Blended(fn, U, param);
    XtLin r;
    r.v = f * x / U(TH);
    r.v_U = f_U * x / U(TH);
    r.v_U(TH) -= r.v / U(TH);
    r.v_x = f / U(TH);
    return r;
}

/* Post-onset n rate at one station: (1 - s) d + s smax(d, C / th). */
struct RateLin {
    double q;
    Vector5d q_U;
};

RateLin Rate(double d, const Vector5d& d_U, double s, double s_n, double C,
             double th) {
    constexpr double ke = 0.1;
    double c = C / th;
    Vector5d c_U = -C / (th * th) * Unit(TH);
    double eps = ke * c;
    Vector5d eps_U = ke * c_U;
    double r = std::sqrt((d - c) * (d - c) + eps * eps);
    double m = 0.5 * (d + c + r);
    Vector5d m_U =
        0.5 * (d_U + c_U + ((d - c) * (d_U - c_U) + eps * eps_U) / r);
    RateLin out;
    out.q = (1.0 - s) * d + s * m;
    out.q_U = (1.0 - s) * d_U + s * m_U + (m - d) * s_n * Unit(AMP);
    return out;
}

}  // namespace

CapturedResidual CapturedIntervalResidual(const Param& param,
                                          const Eigen::Vector2d& x,
                                          const Matrix52d& Uin,
                                          const Eigen::Vector2d& Aux) {
    Matrix52d U = Uin;
    U.row(DS) -= Aux.transpose();
    const Vector5d U1 = U.col(0), U2 = U.col(1);
    const Vector5d Um = 0.5 * (U1 + U2);
    const Vector4d Ut1 = View(U1, kTurb), Ut2 = View(U2, kTurb);
    Param pt = param;
    pt.turb = true;
    Param pl = param;
    pl.turb = false;
    pl.wake = false;
    const Vector10d Z10 = Vector10d::Zero();
    const double th0 = U1(TH), th1 = U2(TH);
    const double ct0 = U1(CT), ct1 = U2(CT);
    const double n0 = U1(AMP), n1 = U2(AMP);

    auto [uk1, uk1_u] = GetUk(U1(UE), param);
    auto [uk2, uk2_u] = GetUk(U2(UE), param);
    double thlog = std::log(th1 / th0);
    Vector10d thlog_U = Two(-Unit(TH) / th0, Unit(TH) / th1);
    double uelog = std::log(uk2 / uk1);
    Vector10d uelog_U =
        Two(-uk1_u(3) / uk1 * Unit(UE), uk2_u(3) / uk2 * Unit(UE));
    double xlog = std::log(x(1) / x(0));
    Vector2d xlog_x(-1.0 / x(0), 1.0 / x(1));
    double dx = x(1) - x(0);
    Vector2d dx_x(-1.0, 1.0);

    auto [upw, upw8] = GetUpw(Ut1, Ut2, pt);
    const Vector10d upw_U = E10(upw8, kTurb);

    auto [H1, H1_U] = GetH(Ut1);
    auto [H2, H2_U] = GetH(Ut2);
    double H = 0.5 * (H1 + H2);
    Vector10d H_U = 0.5 * Two(E5(H1_U, kTurb), E5(H2_U, kTurb));

    auto hs_fn = [](const Vector4d& u, const Param& p) { return GetHs(u, p); };
    auto cf_fn = [](const Vector4d& u, const Param& p) { return GetCf(u, p); };
    auto cdi_fn = [](const Vector4d& u, const Param& p) {
        return GetCDi(u, p);
    };
    auto [Hs1, Hs1_U] = Blended(hs_fn, U1, param);
    auto [Hs2, Hs2_U] = Blended(hs_fn, U2, param);
    auto [Hs, Hs_U] = Up(0.5, Z10, Hs1, Hs1_U, Hs2, Hs2_U);
    double Hslog = std::log(Hs2 / Hs1);
    Vector10d Hslog_U = Two(-Hs1_U / Hs1, Hs2_U / Hs2);

    if (param.simi) {
        thlog = 0.0;
        thlog_U.setZero();
        Hslog = 0.0;
        Hslog_U.setZero();
        uelog = 1.0;
        uelog_U.setZero();
        xlog = 1.0;
        xlog_x.setZero();
        dx = 0.5 * (x(0) + x(1));
        dx_x << 0.5, 0.5;
    }

    auto [Hw1, Hw1_U] = GetHw(Ut1, Aux(0));
    auto [Hw2, Hw2_U] = GetHw(Ut2, Aux(1));
    double Hw = 0.5 * (Hw1 + Hw2);
    Vector10d Hw_U = 0.5 * Two(E5(Hw1_U, kTurb), E5(Hw2_U, kTurb));

    /* ---- shear-lag equation (turbulent closures, target released by r) */
    double salog = std::log(ct1 / ct0);
    Vector10d salog_U = Two(-Unit(CT) / ct0, Unit(CT) / ct1);
    auto [de1, de1_U] = GetDe(Ut1, pt);
    auto [de2, de2_U] = GetDe(Ut2, pt);
    auto [de, de_U] =
        Up(0.5, Z10, de1, E5(de1_U, kTurb), de2, E5(de2_U, kTurb));
    auto [Us1, Us1_U] = GetUs(Ut1, pt);
    auto [Us2, Us2_U] = GetUs(Ut2, pt);
    auto [Us, Us_U] =
        Up(0.5, Z10, Us1, E5(Us1_U, kTurb), Us2, E5(Us2_U, kTurb));
    auto [Hk1, Hk1_U] = GetHk(Ut1, pt);
    auto [Hk2, Hk2_U] = GetHk(Ut2, pt);
    auto [Hk, Hk_U] =
        Up(upw, upw_U, Hk1, E5(Hk1_U, kTurb), Hk2, E5(Hk2_U, kTurb));
    auto [Rt1, Rt1_U] = GetRet(Ut1, pt);
    auto [Rt2, Rt2_U] = GetRet(Ut2, pt);
    auto [Ret, Ret_U] =
        Up(0.5, Z10, Rt1, E5(Rt1_U, kTurb), Rt2, E5(Rt2_U, kTurb));
    auto [cft1, cft1_U] = GetCf(Ut1, pt);
    auto [cft2, cft2_U] = GetCf(Ut2, pt);
    auto [cft, cft_U] =
        Up(upw, upw_U, cft1, E5(cft1_U, kTurb), cft2, E5(cft2_U, kTurb));
    double dsa = 0.5 * (U(DS, 0) + U(DS, 1));
    /* GetUq works on 8-entry (two 4-state) linearizations: map through the
     * turbulent view, whose rows are identical in both layouts. */
    auto to8 = [](const Vector10d& v) {
        Vector8d o;
        o << v(TH), v(DS), v(CT), v(UE), v(5 + TH), v(5 + DS), v(5 + CT),
            v(5 + UE);
        return o;
    };
    Vector8d dsa8 = Vector8d::Zero();
    dsa8(1) = 0.5;
    dsa8(5) = 0.5;
    auto [uq, uq8] =
        GetUq(dsa, dsa8, cft, to8(cft_U), Hk, to8(Hk_U), Ret, to8(Ret_U), pt);
    Vector10d uq_U = E10(uq8, kTurb);

    std::array<double, 2> tg{};
    std::array<Vector5d, 2> tg_U{};
    for (int k = 0; k < 2; ++k) {
        const Vector5d& Uk = k == 0 ? U1 : U2;
        auto [ceq, ceq_U] = GetCteq(View(Uk, kTurb), pt);
        if (param.wake) {
            tg[k] = ceq;
            tg_U[k] = E5(ceq_U, kTurb);
            continue;
        }
        auto [ctr, ctr_U] = GetCttr(View(Uk, kTurb), pt);
        auto [r, r_n] = Release(Uk(AMP), param);
        tg[k] = (1.0 - r) * ctr + r * ceq;
        tg_U[k] = (1.0 - r) * E5(ctr_U, kTurb) + r * E5(ceq_U, kTurb);
        tg_U[k](AMP) += (ceq - ctr) * r_n;
    }
    auto [tgt, tgt_U] = Up(upw, upw_U, tg[0], tg_U[0], tg[1], tg_U[1]);
    auto [saa, saa_U] = Up(upw, upw_U, ct0, Unit(CT), ct1, Unit(CT));
    double rb = 1.0;
    Vector10d rb_U = Z10;
    if (!param.wake) {
        auto [ra, ra_n] = Release(n0, param);
        auto [rc, rc_n] = Release(n1, param);
        rb = 0.5 * (ra + rc);
        rb_U = 0.5 * Two(ra_n * Unit(AMP), rc_n * Unit(AMP));
    }
    double Clag = param.SlagK / param.GB / (1.0 + Us);
    Vector10d Clag_U = -Clag / (1.0 + Us) * Us_U;
    double ald = param.wake ? param.Dlr : 1.0;
    double pg = 2.0 * de * (uq * dx - uelog) * param.Cuq;
    Vector10d pg_U = 2.0 * de_U * (uq * dx - uelog) * param.Cuq +
                     2.0 * de * (uq_U * dx - uelog_U) * param.Cuq;
    double Rlag = Clag * (tgt - ald * saa) * dx - 2.0 * de * salog + rb * pg;
    Vector10d Rlag_U = Clag_U * (tgt - ald * saa) * dx +
                       Clag * (tgt_U - ald * saa_U) * dx - 2.0 * de_U * salog -
                       2.0 * de * salog_U + rb_U * pg + rb * pg_U;
    Vector2d Rlag_x =
        Clag * (tgt - ald * saa) * dx_x + rb * 2.0 * de * uq * param.Cuq * dx_x;
    if (param.simi) {
        /* The ODE degenerates at the stagnation point (dx ~ 1e-6): pin ct
         * algebraically, with the 2 de scale of the transport term. */
        double lg = std::log(ct1 / tg[1]);
        Rlag = 2.0 * de * lg;
        Vector10d lg_U = Z10;
        lg_U.tail<5>() = Unit(CT) / ct1 - tg_U[1] / tg[1];
        Rlag_U = 2.0 * de_U * lg + 2.0 * de * lg_U;
        Rlag_x.setZero();
    }

    /* ---- amplification equation */
    double Ramp;
    Vector10d Ramp_U;
    Vector2d Ramp_x = Vector2d::Zero();
    if (param.wake) {
        Ramp = n1 - n0;
        Ramp_U = Two(-Unit(AMP), Unit(AMP));
    } else if (param.simi) {
        Ramp = n0 + n1;
        Ramp_U = Two(Unit(AMP), Unit(AMP));
    } else {
        auto [d1, d1_U] = GetDamp(View(U1, kLam), pl);
        auto [d2, d2_U] = GetDamp(View(U2, kLam), pl);
        double C = RateC(param);
        auto [s1, s1_n] = Growth(n0, param);
        auto [s2, s2_n] = Growth(n1, param);
        RateLin q1 = Rate(d1, E5(d1_U, kLam), s1, s1_n, C, th0);
        RateLin q2 = Rate(d2, E5(d2_U, kLam), s2, s2_n, C, th1);
        double damp = 0.5 * (q1.q + q2.q);
        Vector10d damp_U = 0.5 * Two(q1.q_U, q2.q_U);
        Ramp = n1 - n0 - damp * dx;
        Ramp_U = Two(-Unit(AMP), Unit(AMP)) - damp_U * dx;
        Ramp_x = -damp * dx_x;
    }

    /* ---- momentum equation (blended cf) */
    auto [Ms1, Ms1_U] = GetMach2(Ut1, param);
    auto [Ms2, Ms2_U] = GetMach2(Ut2, param);
    double Ms = 0.5 * (Ms1 + Ms2);
    Vector10d Ms_U = 0.5 * Two(E5(Ms1_U, kTurb), E5(Ms2_U, kTurb));
    XtLin c1 = BlendedXt(cf_fn, U1, x(0), param);
    XtLin c2 = BlendedXt(cf_fn, U2, x(1), param);
    XtLin cm = BlendedXt(cf_fn, Um, 0.5 * (x(0) + x(1)), param);
    double cfxt = 0.25 * c1.v + 0.5 * cm.v + 0.25 * c2.v;
    Vector10d cfxt_U = 0.25 * Two(c1.v_U + cm.v_U, cm.v_U + c2.v_U);
    Vector2d cfxt_x(0.25 * (c1.v_x + cm.v_x), 0.25 * (cm.v_x + c2.v_x));
    double Rmom = thlog + (2.0 + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt;
    Vector10d Rmom_U = thlog_U + (H_U + Hw_U - Ms_U) * uelog +
                       (2.0 + H + Hw - Ms) * uelog_U - 0.5 * xlog * cfxt_U;
    Vector2d Rmom_x = -0.5 * xlog_x * cfxt - 0.5 * xlog * cfxt_x;

    /* ---- shape equation (blended Hs, cf, cDi) */
    XtLin D1 = BlendedXt(cdi_fn, U1, x(0), param);
    XtLin D2 = BlendedXt(cdi_fn, U2, x(1), param);
    auto [cDixt, cDixt_U] = Up(upw, upw_U, D1.v, D1.v_U, D2.v, D2.v_U);
    Vector2d cDixt_x((1.0 - upw) * D1.v_x, upw * D2.v_x);
    auto [cfxtu, cfxtu_U] = Up(upw, upw_U, c1.v, c1.v_U, c2.v, c2.v_U);
    Vector2d cfxtu_x((1.0 - upw) * c1.v_x, upw * c2.v_x);
    auto [Hss1, Hss1_U] = GetHss(Ut1, param);
    auto [Hss2, Hss2_U] = GetHss(Ut2, param);
    auto [Hss, Hss_U] =
        Up(0.5, Z10, Hss1, E5(Hss1_U, kTurb), Hss2, E5(Hss2_U, kTurb));
    double Rshape = Hslog + (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog +
                    xlog * (0.5 * cfxtu - cDixt);
    Vector10d Rshape_U =
        Hslog_U +
        (2.0 * Hss_U / Hs - 2.0 * Hss / (Hs * Hs) * Hs_U - H_U - Hw_U) * uelog +
        (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog_U +
        xlog * (0.5 * cfxtu_U - cDixt_U);
    Vector2d Rshape_x =
        xlog_x * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x - cDixt_x);

    CapturedResidual out;
    out.R << Rmom, Rshape, Ramp, Rlag;
    out.R_U.row(0) = Rmom_U.transpose();
    out.R_U.row(1) = Rshape_U.transpose();
    out.R_U.row(2) = Ramp_U.transpose();
    out.R_U.row(3) = Rlag_U.transpose();
    out.R_x.row(0) = Rmom_x.transpose();
    out.R_x.row(1) = Rshape_x.transpose();
    out.R_x.row(2) = Ramp_x.transpose();
    out.R_x.row(3) = Rlag_x.transpose();
    return out;
}

namespace {

struct Stagnation5 {
    Vector5d Ust;
    Eigen::Matrix<double, 5, 10> Ust_U;
    Eigen::Matrix<double, 5, 2> Ust_x;
    double xst;
};

/* 5-state stagnation extrapolation (ue quadratic, ct from the first node). */
Stagnation5 StagnationState5(const Matrix52d& U, const Vector2d& x) {
    const Vector5d U1 = U.col(0), U2 = U.col(1);
    double x1 = x(0), x2 = x(1);
    double dx = x2 - x1;
    Vector2d dx_x(-1.0, 1.0);
    double rx = x2 / x1;
    Vector2d rx_x(-rx / x1, 1.0 / x1);
    double w1 = x2 / dx;
    Vector2d w1_x = -w1 / dx * dx_x + Vector2d(0.0, 1.0) / dx;
    double w2 = -x1 / dx;
    Vector2d w2_x = -w2 / dx * dx_x + Vector2d(-1.0, 0.0) / dx;
    double wk1 = rx / dx;
    Vector2d wk1_x = rx_x / dx - wk1 / dx * dx_x;
    double wk2 = -1.0 / (rx * dx);
    Vector2d wk2_x = -wk2 * (Vector2d(rx_x / rx) + Vector2d(dx_x / dx));
    double K = wk1 * U1(UE) + wk2 * U2(UE);
    Vector2d K_x = U1(UE) * wk1_x + U2(UE) * wk2_x;
    Stagnation5 s;
    s.xst = 1e-6;
    s.Ust = U1 * w1 + U2 * w2;
    s.Ust(UE) = K * s.xst;
    Eigen::Matrix<double, 5, 5> E = Eigen::Matrix<double, 5, 5>::Identity();
    E(UE, UE) = 0.0;
    s.Ust_U << w1 * E, w2 * E;
    s.Ust_U(UE, UE) = wk1 * s.xst;
    s.Ust_U(UE, 5 + UE) = wk2 * s.xst;
    s.Ust_x = U1 * w1_x.transpose() + U2 * w2_x.transpose();
    s.Ust_x.row(UE) = K_x.transpose() * s.xst;
    /* ct: first node's value; extrapolation can go negative */
    s.Ust(CT) = U1(CT);
    s.Ust_U.row(CT).setZero();
    s.Ust_U(CT, CT) = 1.0;
    s.Ust_x.row(CT).setZero();
    return s;
}

Vector5d NodeState(const FoilState& state, int j) {
    Vector5d u;
    u << state.glob.U.col(j), state.glob.Un(j);
    return u;
}

struct Assembly {
    VectorXd R;        /* 4 Nsys */
    SparseMatrixd R_U; /* 4 Nsys x 5 Nsys */
    SparseMatrixd R_x; /* 4 Nsys x Nsys */
};

/* Global captured system, including the stagnation-movement chain in the ue
 * columns of the two stagnation-panel nodes (as BuildGlobSys does). */
Assembly BuildGlobSysCaptured(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    Assembly a;
    a.R = VectorXd::Zero(4 * Nsys);
    std::vector<Triplet> tU, tx;
    auto put = [&](int node, const Vector4d& Ri, const Eigen::MatrixXd& Ri_U,
                   const std::vector<int>& J, const Eigen::MatrixXd* Ri_x) {
        a.R.segment<4>(4 * node) += Ri;
        for (size_t j = 0; j < J.size(); ++j) {
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 5; ++c) {
                    double v = Ri_U(r, 5 * j + c);
                    if (v != 0.0)
                        tU.emplace_back(4 * node + r, 5 * J[j] + c, v);
                }
                if (Ri_x) {
                    double v = (*Ri_x)(r, j);
                    if (v != 0.0) tx.emplace_back(4 * node + r, J[j], v);
                }
            }
        }
    };

    for (int si = 0; si < 3; ++si) {
        const std::vector<int>& Is = state.vsol.Is[si];
        const int N = static_cast<int>(Is.size());
        Param param = BuildParam(state, si);
        auto xi = [&](int i) { return state.isol.xi(Is[i]); };
        auto aux = [&](int i) { return si == 2 ? state.vsol.wgap(i) : 0.0; };
        int i0 = (si < 2 && xi(0) < 1e-8 * xi(N - 1)) ? 1 : 0;
        if (si < 2) {
            Matrix52d Up2;
            Up2.col(0) = NodeState(state, Is[i0]);
            Up2.col(1) = NodeState(state, Is[i0 + 1]);
            Stagnation5 st =
                StagnationState5(Up2, Vector2d(xi(i0), xi(i0 + 1)));
            param.turb = false;
            param.simi = true;
            Matrix52d Ust2;
            Ust2.col(0) = st.Ust;
            Ust2.col(1) = st.Ust;
            CapturedResidual r1 =
                CapturedIntervalResidual(param, Vector2d(st.xst, st.xst), Ust2,
                                         Vector2d(aux(i0), aux(i0)));
            param.simi = false;
            Eigen::Matrix<double, 4, 5> R1_Ust =
                r1.R_U.leftCols<5>() + r1.R_U.rightCols<5>();
            std::vector<int> J = {Is[i0], Is[i0 + 1]};
            if (i0 == 1) {
                /* first node sits on the stagnation point */
                const std::array<int, 4> rows = {TH, DS, CT, AMP};
                Vector4d Rs;
                Eigen::MatrixXd Rs_U = Eigen::MatrixXd::Zero(4, 15);
                Eigen::MatrixXd Rs_x = Eigen::MatrixXd::Zero(4, 3);
                Vector5d U0 = NodeState(state, Is[0]);
                for (int k = 0; k < 4; ++k) {
                    Rs(k) = U0(rows[k]) - st.Ust(rows[k]);
                    Rs_U(k, rows[k]) = 1.0;
                    Rs_U.block(k, 5, 1, 10) = -st.Ust_U.row(rows[k]);
                    Rs_x.block(k, 1, 1, 2) = -st.Ust_x.row(rows[k]);
                }
                put(Is[0], Rs, Rs_U, {Is[0], J[0], J[1]}, &Rs_x);
            }
            Eigen::MatrixXd R1_U = R1_Ust * st.Ust_U;
            Eigen::MatrixXd R1_x = R1_Ust * st.Ust_x;
            put(Is[i0], r1.R, R1_U, J, &R1_x);
        } else {
            /* wake initialization: th, ds sums; theta-weighted ct; n = 0 */
            int il = state.vsol.Is[0].back(), iu = state.vsol.Is[1].back();
            int iw = Is[0];
            Vector5d Ul = NodeState(state, il), Uu = NodeState(state, iu),
                     Uw = NodeState(state, iw);
            double hTE = GetTEInfo(state.foil.x).hTE;
            double thsum = Ul(TH) + Uu(TH);
            double ctw = (Ul(CT) * Ul(TH) + Uu(CT) * Uu(TH)) / thsum;
            Vector4d R;
            R << Uw(TH) - thsum, Uw(DS) - (Ul(DS) + Uu(DS) + hTE), Uw(AMP),
                Uw(CT) - ctw;
            Eigen::MatrixXd R_U = Eigen::MatrixXd::Zero(4, 15);
            for (int b : {0, 5}) {
                R_U(0, b + TH) = -1.0;
                R_U(1, b + DS) = -1.0;
            }
            R_U(0, 10 + TH) = 1.0;
            R_U(1, 10 + DS) = 1.0;
            R_U(2, 10 + AMP) = 1.0;
            R_U(3, 10 + CT) = 1.0;
            R_U(3, CT) = -Ul(TH) / thsum;
            R_U(3, TH) = -(Ul(CT) - ctw) / thsum;
            R_U(3, 5 + CT) = -Uu(TH) / thsum;
            R_U(3, 5 + TH) = -(Uu(CT) - ctw) / thsum;
            put(iw, R, R_U, {il, iu, iw}, nullptr);
            param.turb = true;
            param.wake = true;
        }
        for (int i = i0 + 1; i < N; ++i) {
            Matrix52d Ui;
            Ui.col(0) = NodeState(state, Is[i - 1]);
            Ui.col(1) = NodeState(state, Is[i]);
            CapturedResidual ri =
                CapturedIntervalResidual(param, Vector2d(xi(i - 1), xi(i)), Ui,
                                         Vector2d(aux(i - 1), aux(i)));
            Eigen::MatrixXd RU = ri.R_U, Rx = ri.R_x;
            put(Is[i], ri.R, RU, {Is[i - 1], Is[i]}, &Rx);
        }
    }
    a.R_U.resize(4 * Nsys, 5 * Nsys);
    a.R_U.setFromTriplets(tU.begin(), tU.end());
    a.R_x.resize(4 * Nsys, Nsys);
    a.R_x.setFromTriplets(tx.begin(), tx.end());

    /* stagnation-point movement: R_ue += R_x x_st st_ue */
    VectorXd x_st(Nsys);
    const int N = state.foil.N;
    for (int i = 0; i < Nsys; ++i) {
        x_st(i) = i < N ? -state.isol.sgnue(i) : -1.0;
    }
    VectorXd R_st = a.R_x * x_st;
    std::vector<Triplet> tst;
    for (int k = 0; k < 2; ++k) {
        int col = 5 * state.isol.Istag[k] + UE;
        for (int r = 0; r < 4 * Nsys; ++r) {
            if (R_st(r) != 0.0) {
                tst.emplace_back(r, col, R_st(r) * state.isol.sstag_ue(k));
            }
        }
    }
    SparseMatrixd S(4 * Nsys, 5 * Nsys);
    S.setFromTriplets(tst.begin(), tst.end());
    a.R_U += S;
    return a;
}

/* Newton step for [th, ds, ct, ue, n] with the ue coupling rows. */
bool SolveGlobCaptured(FoilState& state, const Assembly& a,
                       Eigen::MatrixXd& dU) {
    const int Nsys = state.glob.Nsys;
    VectorXd ue = state.glob.U.row(UE).transpose();
    VectorXd ds = state.glob.U.row(DS).transpose();
    double uemax = ue.cwiseAbs().maxCoeff();
    for (int i = 0; i < Nsys; ++i) ue(i) = std::max(ue(i), 1e-10 * uemax);
    VectorXd ueinv = GetUeInv(state);
    VectorXd R(5 * Nsys);
    R.head(4 * Nsys) = a.R;
    R.tail(Nsys) = ue - (ueinv + state.vsol.ue_m * ds.cwiseProduct(ue));

    std::vector<Triplet> t;
    t.reserve(a.R_U.nonZeros() + 2 * Nsys * Nsys);
    for (int k = 0; k < a.R_U.outerSize(); ++k) {
        for (SparseMatrixd::InnerIterator it(a.R_U, k); it; ++it) {
            t.emplace_back(it.row(), it.col(), it.value());
        }
    }
    for (int i = 0; i < Nsys; ++i) {
        int row = 4 * Nsys + i;
        t.emplace_back(row, 5 * i + UE, 1.0);
        for (int j = 0; j < Nsys; ++j) {
            double m = state.vsol.ue_m(i, j);
            if (m == 0.0) continue;
            t.emplace_back(row, 5 * j + UE, -m * ds(j));
            t.emplace_back(row, 5 * j + DS, -m * ue(j));
        }
    }
    SparseMatrixd J(5 * Nsys, 5 * Nsys);
    J.setFromTriplets(t.begin(), t.end());
    Eigen::SparseLU<SparseMatrixd> lu;
    lu.compute(J);
    if (lu.info() != Eigen::Success) return false;
    VectorXd dV = -lu.solve(R);
    if (lu.info() != Eigen::Success || !dV.allFinite()) return false;
    dU.resize(5, Nsys);
    for (int i = 0; i < Nsys; ++i) dU.col(i) = dV.segment<5>(5 * i);
    return true;
}

/* Damped update keeping th, ds, ct positive and changes bounded. */
void UpdateStateCaptured(FoilState& state, const Eigen::MatrixXd& dU) {
    const int Nsys = state.glob.Nsys;
    Eigen::MatrixXd U(5, Nsys);
    U.topRows<4>() = state.glob.U;
    U.row(AMP) = state.glob.Un.transpose();
    const Param& p = state.param;
    double omega = 1.0;
    for (int k : {TH, DS}) {
        double fmin = (dU.row(k).array() / U.row(k).array()).minCoeff();
        if (fmin < -0.5) omega = std::min(omega, std::abs(0.5 / fmin));
    }
    /* ct decreases are limited only where ct matters (released or
     * blended); n changes only where the closure blend is unsaturated. */
    std::vector<bool> is_wake(Nsys, false);
    for (int j : state.vsol.Is[2]) is_wake[j] = true;
    double rmin = 0.0, dnmax = 0.0;
    for (int j = 0; j < Nsys; ++j) {
        double n = U(AMP, j);
        bool act = is_wake[j] ||
                   std::max(Release(n, p).first, Gamma(n, p).first) > 1e-3;
        if (act) rmin = std::min(rmin, dU(CT, j) / U(CT, j));
        if (Gamma(n, p).first < 0.999) {
            dnmax = std::max(dnmax, std::abs(dU(AMP, j)));
        }
    }
    if (rmin < -0.8) omega = std::min(omega, std::abs(0.8 / rmin));
    if (dnmax > 2.0) omega = std::min(omega, 2.0 / dnmax);
    double fmax = dU.row(UE).cwiseAbs().maxCoeff() / state.oper.Vinf;
    if (fmax > 0.2) omega = std::min(omega, 0.2 / fmax);

    Eigen::MatrixXd Unew = U + omega * dU;
    for (int j = 0; j < Nsys; ++j) {
        if (Unew(CT, j) < 0.2 * U(CT, j)) Unew(CT, j) = 0.2 * U(CT, j);
    }
    state.glob.U = Unew.topRows<4>();
    state.glob.Un = Unew.row(AMP).transpose();
    /* fix bad Hk, as UpdateState does */
    for (int si = 0; si < 3; ++si) {
        double Hkmin = si == 2 ? 1.00005 : 1.02;
        Param param = BuildParam(state, si);
        for (int j : state.vsol.Is[si]) {
            Vector4d Uj = state.glob.U.col(j);
            double Hk = GetHk(Uj, param).value;
            if (Hk < Hkmin) {
                state.glob.U(DS, j) += 2.0 * (Hkmin - Hk) * state.glob.U(DS, j);
            }
        }
    }
}

/* Converged cut-cell [th, ds, sa, ue] -> captured [th, ds, ct, ue] + n. */
void ConvertState(FoilState& state) {
    const int Nsys = state.glob.Nsys;
    const Param& P = state.param;
    MatrixXd U4 = state.glob.U;
    state.glob.Un = VectorXd::Zero(Nsys);
    double C = RateC(P);
    for (int si = 0; si < 3; ++si) {
        const std::vector<int>& Is = state.vsol.Is[si];
        Param param = BuildParam(state, si);
        for (size_t k = 0; k < Is.size(); ++k) {
            int j = Is[k];
            if (si < 2 && !state.vsol.turb[j]) {
                Param pt = param;
                pt.turb = true;
                state.glob.U(CT, j) = GetCttr(Vector4d(U4.col(j)), pt).value;
            }
            if (si == 2 || k == 0) continue;
            /* march n with the captured amplification equation */
            Param pl = param;
            pl.turb = false;
            pl.wake = false;
            int jm = Is[k - 1];
            Vector4d Um(U4(TH, jm), U4(DS, jm), state.glob.Un(jm), U4(UE, jm));
            Vector4d Uj(U4(TH, j), U4(DS, j), Um(2), U4(UE, j));
            double dx = state.isol.xi(j) - state.isol.xi(jm);
            double d1 = GetDamp(Um, pl).value;
            auto [sm, sm_n] = Growth(Um(2), P);
            (void)sm_n;
            double q1 = Rate(d1, Vector5d::Zero(), sm, 0.0, C, Um(0)).q;
            for (int it = 0; it < 60; ++it) {
                auto [d2, d2_U] = GetDamp(Uj, pl);
                auto [s2, s2_n] = Growth(Uj(2), P);
                RateLin q2 = Rate(d2, E5(d2_U, kLam), s2, s2_n, C, Uj(0));
                double Ra = Uj(2) - Um(2) - 0.5 * (q1 + q2.q) * dx;
                Uj(2) -= Ra / (1.0 - 0.5 * q2.q_U(AMP) * dx);
                if (std::abs(Ra) < 1e-12) break;
            }
            state.glob.Un(j) = Uj(2);
        }
    }
}

}  // namespace

void SolveCaptured(FoilState& state) {
    ConvertState(state);
    state.glob.conv = false;
    const int nNewton = std::max(state.param.niglob, 40);
    for (int it = 0; it < nNewton; ++it) {
        Assembly a = BuildGlobSysCaptured(state);
        double Rn = a.R.norm();
        if (state.param.verb >= 1) {
            std::cout << "  captured Newton " << it << ", |R| = " << Rn << '\n';
        }
        if (!std::isfinite(Rn)) break;
        if (Rn < std::max(state.param.rtol, 1e-10)) {
            state.glob.conv = true;
            break;
        }
        Eigen::MatrixXd dU;
        if (!SolveGlobCaptured(state, a, dU)) break;
        UpdateStateCaptured(state, dU);
        StagPointMove(state);
    }
    /* report transition where the closure blend crosses 1/2 */
    for (int si = 0; si < 3; ++si) {
        for (int j : state.vsol.Is[si]) {
            state.vsol.turb[j] =
                si == 2 || Gamma(state.glob.Un(j), state.param).first > 0.5;
        }
    }
}

}  // namespace scirocco
