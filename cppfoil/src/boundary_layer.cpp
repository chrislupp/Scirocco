// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/boundary_layer.h"

#include <cmath>
#include "cppfoil/properties/bl_properties.h"

namespace cppfoil {

std::pair<double, double> ThwaitesInit(double K, double nu) {
    // Thwaites correlation for stagnation point
    // K = due/dx at stagnation
    // nu = kinematic viscosity

    double th = std::sqrt(0.075 * nu / std::max(K, 1e-10));
    double H = 2.55;  // Shape parameter at stagnation
    double ds = H * th;

    return {th, ds};
}

Param BuildParam(const FoilState& state, int si) {
    Param param = state.param;

    // Set surface-specific flags
    param.wake = (si == 2);
    param.simi = false;
    param.turb = false;

    return param;
}

void StationParam(const FoilState& state, Param& param, int i) {
    // Update turbulence flag
    int N = state.foil.N;
    if (i < static_cast<int>(state.vsol.turb.size())) {
        param.turb = state.vsol.turb[i];
    }

    // Wake flag
    param.wake = (i >= N);
}

void InitBoundaryLayer(FoilState& state) {
    // Initialize boundary layer from inviscid solution
    // This is a simplified initialization

    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    state.glob.Nsys = Nsys;
    state.glob.U.resize(4, Nsys);
    state.glob.U.setZero();

    // Initialize with small values
    double th0 = 0.001 * state.geom.chord;
    double ds0 = 2.5 * th0;
    double sa0 = 0.0;

    for (int i = 0; i < Nsys; ++i) {
        state.glob.U(0, i) = th0;  // theta
        state.glob.U(1, i) = ds0;  // delta*
        state.glob.U(2, i) = sa0;  // sa (amp or ctau)

        // Edge velocity from inviscid
        if (i < N) {
            state.glob.U(3, i) = state.isol.gam(i) * state.isol.sgnue(i);
        } else if (state.wake.N > 0) {
            state.glob.U(3, i) = state.isol.uewi(i - N);
        }
    }

    // Initialize turbulence flags (all laminar initially)
    state.vsol.turb.resize(Nsys, false);
}

ResidualLin ResidualStation(const Param& param, const Vector2d& x,
                             const Matrix42d& Uin, const Vector2d& Aux) {
    // Compute residual at a boundary layer station
    // U = [th, ds, sa, ue] at two stations
    // Aux = wake gap at each station

    ResidualLin result;
    result.value = Vector3d::Zero();
    result.lin_U = Matrix38d::Zero();
    result.lin_x = Matrix32d::Zero();

    // Copy and modify for wake gap
    Matrix42d U = Uin;
    U.row(1) -= Aux.transpose();  // Remove wake gap from delta*

    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);

    double th1 = std::max(U1(0), 1e-12);
    double th2 = std::max(U2(0), 1e-12);
    double sa1 = std::max(U1(2), 1e-12);
    double sa2 = std::max(U2(2), 1e-12);

    // Compressibility-corrected velocities
    auto [uk1, uk1_U] = GetUk(U1(3), param);
    auto [uk2, uk2_U] = GetUk(U2(3), param);

    // Log changes
    double thlog = std::log(th2 / th1);
    Vector8d thlog_U = Vector8d::Zero();
    thlog_U(0) = -1.0 / th1;
    thlog_U(4) = 1.0 / th2;

    double uelog = std::log(std::max(uk2, 1e-12) / std::max(uk1, 1e-12));
    Vector8d uelog_U = Vector8d::Zero();
    uelog_U(3) = -uk1_U(3) / std::max(uk1, 1e-12);
    uelog_U(7) = uk2_U(3) / std::max(uk2, 1e-12);

    double xlog = std::log(std::max(x(1), 1e-12) / std::max(x(0), 1e-12));
    Vector2d xlog_x(-1.0 / std::max(x(0), 1e-12), 1.0 / std::max(x(1), 1e-12));

    double dx = x(1) - x(0);
    Vector2d dx_x(-1.0, 1.0);

    // Similarity station special handling
    if (param.simi) {
        thlog = 0.0;
        thlog_U.setZero();
        uelog = 1.0;
        uelog_U.setZero();
        xlog = 1.0;
        xlog_x.setZero();
        dx = 0.5 * (x(0) + x(1));
        dx_x << 0.5, 0.5;
    }

    // Upwinding factor
    auto [upw, upw_U] = GetUpw(U1, U2, param);

    // Shape parameters
    auto [H1, H1_U1] = GetH(U1);
    auto [H2, H2_U2] = GetH(U2);
    double H = 0.5 * (H1 + H2);
    Vector8d H_U;
    H_U << 0.5 * H1_U1, 0.5 * H2_U2;

    // Kinetic energy shape parameter Hs
    auto [Hs1, Hs1_U1] = GetHs(U1, param);
    auto [Hs2, Hs2_U2] = GetHs(U2, param);
    auto [Hs, Hs_U] = Upwind(0.5, Vector8d::Zero(), Hs1, Hs1_U1, Hs2, Hs2_U2);

    double Hslog = std::log(std::max(Hs2, 1e-12) / std::max(Hs1, 1e-12));
    Vector8d Hslog_U;
    Hslog_U << -Hs1_U1 / std::max(Hs1, 1e-12), Hs2_U2 / std::max(Hs2, 1e-12);

    if (param.simi) {
        Hslog = 0.0;
        Hslog_U.setZero();
    }

    // Wake shape parameter with gap
    auto [Hw1, Hw1_U1] = GetHw(U1, Aux(0));
    auto [Hw2, Hw2_U2] = GetHw(U2, Aux(1));
    double Hw = 0.5 * (Hw1 + Hw2);
    Vector8d Hw_U;
    Hw_U << 0.5 * Hw1_U1, 0.5 * Hw2_U2;

    // Squared Mach number
    auto [Ms1, Ms1_U1] = GetMach2(U1, param);
    auto [Ms2, Ms2_U2] = GetMach2(U2, param);
    auto [Ms, Ms_U] = Upwind(0.5, Vector8d::Zero(), Ms1, Ms1_U1, Ms2, Ms2_U2);

    // cf * x / theta
    auto [cfxt1, cfxt1_U1, cfxt1_x1] = GetCfxt(U1, x(0), param);
    auto [cfxt2, cfxt2_U2, cfxt2_x2] = GetCfxt(U2, x(1), param);
    double xm = 0.5 * (x(0) + x(1));
    Vector4d Um = 0.5 * (U1 + U2);
    auto [cfxtm, cfxtm_Um, cfxtm_xm] = GetCfxt(Um, xm, param);
    double cfxt = 0.25 * cfxt1 + 0.5 * cfxtm + 0.25 * cfxt2;
    Vector8d cfxt_U;
    cfxt_U << 0.25 * cfxt1_U1 + 0.25 * cfxtm_Um, 0.25 * cfxtm_Um + 0.25 * cfxt2_U2;
    Vector2d cfxt_x(0.25 * cfxt1_x1 + 0.25 * cfxtm_xm, 0.25 * cfxtm_xm + 0.25 * cfxt2_x2);

    // Momentum equation
    double Rmom = thlog + (2.0 + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt;
    Vector8d Rmom_U = thlog_U + (H_U + Hw_U - Ms_U) * uelog +
                       (2.0 + H + Hw - Ms) * uelog_U - 0.5 * xlog * cfxt_U;
    Vector2d Rmom_x = -0.5 * xlog_x * cfxt - 0.5 * xlog * cfxt_x;

    // cDi * x / theta (dissipation)
    auto [cDixt1, cDixt1_U1, cDixt1_x1] = GetCDixt(U1, x(0), param);
    auto [cDixt2, cDixt2_U2, cDixt2_x2] = GetCDixt(U2, x(1), param);
    auto [cDixt, cDixt_U] = Upwind(upw, upw_U, cDixt1, cDixt1_U1, cDixt2, cDixt2_U2);
    Vector2d cDixt_x((1.0 - upw) * cDixt1_x1, upw * cDixt2_x2);

    // Upwinded cf*x/theta
    auto [cfxtu, cfxtu_U] = Upwind(upw, upw_U, cfxt1, cfxt1_U1, cfxt2, cfxt2_U2);
    Vector2d cfxtu_x((1.0 - upw) * cfxt1_x1, upw * cfxt2_x2);

    // Density shape parameter Hss
    auto [Hss1, Hss1_U1] = GetHss(U1, param);
    auto [Hss2, Hss2_U2] = GetHss(U2, param);
    auto [Hss, Hss_U] = Upwind(0.5, Vector8d::Zero(), Hss1, Hss1_U1, Hss2, Hss2_U2);

    // Shape parameter equation
    double Rshape = Hslog + (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog +
                    xlog * (0.5 * cfxtu - cDixt);
    Vector8d Rshape_U = Hslog_U +
                        (2.0 * Hss_U / Hs - 2.0 * Hss / (Hs * Hs) * Hs_U - H_U - Hw_U) * uelog +
                        (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog_U +
                        xlog * (0.5 * cfxtu_U - cDixt_U);
    Vector2d Rshape_x = xlog_x * (0.5 * cfxtu - cDixt) +
                         xlog * (0.5 * cfxtu_x - cDixt_x);

    // Amplification or shear lag equation
    double Rlag;
    Vector8d Rlag_U;
    Vector2d Rlag_x;

    if (param.turb) {
        // Turbulent: shear lag equation
        double salog = std::log(sa2 / sa1);
        Vector8d salog_U = Vector8d::Zero();
        salog_U(2) = -1.0 / sa1;
        salog_U(6) = 1.0 / sa2;

        auto [de1, de1_U1] = GetDe(U1, param);
        auto [de2, de2_U2] = GetDe(U2, param);
        auto [de, de_U] = Upwind(0.5, Vector8d::Zero(), de1, de1_U1, de2, de2_U2);

        auto [Us1, Us1_U1] = GetUs(U1, param);
        auto [Us2, Us2_U2] = GetUs(U2, param);
        auto [Us, Us_U] = Upwind(0.5, Vector8d::Zero(), Us1, Us1_U1, Us2, Us2_U2);

        auto [cteq1, cteq1_U1] = GetCteq(U1, param);
        auto [cteq2, cteq2_U2] = GetCteq(U2, param);
        auto [cteq, cteq_U] = Upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2);

        // Root shear coefficient upwinded
        Vector4d sa1_U1 = Vector4d::Zero();
        sa1_U1(2) = 1.0;
        Vector4d sa2_U2 = Vector4d::Zero();
        sa2_U2(2) = 1.0;
        auto [saa, saa_U] = Upwind(upw, upw_U, sa1, sa1_U1, sa2, sa2_U2);

        double Klag = param.SlagK;
        double beta = param.GB;
        double Clag = Klag / beta / (1.0 + Us);
        Vector8d Clag_U = -Clag / (1.0 + Us) * Us_U;

        double ald = param.wake ? param.Dlr : 1.0;

        Rlag = Clag * (cteq - ald * saa) * dx - 2.0 * de * salog;
        Rlag_U = Clag_U * (cteq - ald * saa) * dx +
                 Clag * (cteq_U - ald * saa_U) * dx -
                 2.0 * de_U * salog - 2.0 * de * salog_U;
        Rlag_x = Clag * (cteq - ald * saa) * dx_x;
    } else {
        // Laminar: amplification factor equation
        if (param.simi) {
            Rlag = sa1 + sa2;
            Rlag_U = Vector8d::Zero();
            Rlag_U(2) = 1.0;
            Rlag_U(6) = 1.0;
            Rlag_x = Vector2d::Zero();
        } else {
            auto [damp1, damp1_U1] = GetDamp(U1, param);
            auto [damp2, damp2_U2] = GetDamp(U2, param);
            auto [damp, damp_U] = Upwind(0.5, Vector8d::Zero(), damp1, damp1_U1, damp2, damp2_U2);

            Rlag = sa2 - sa1 - damp * dx;
            Rlag_U = Vector8d::Zero();
            Rlag_U(2) = -1.0;
            Rlag_U(6) = 1.0;
            Rlag_U -= damp_U * dx;
            Rlag_x = -damp * dx_x;
        }
    }

    // Assemble results
    result.value << Rmom, Rshape, Rlag;
    result.lin_U.row(0) = Rmom_U;
    result.lin_U.row(1) = Rshape_U;
    result.lin_U.row(2) = Rlag_U;
    result.lin_x.row(0) = Rmom_x;
    result.lin_x.row(1) = Rshape_x;
    result.lin_x.row(2) = Rlag_x;

    return result;
}

StagnationResult StagnationState(const Matrix42d& U, const Vector2d& x) {
    // Extrapolates two states to stagnation point
    // Fits quadratic to edge velocity (0 at x=0), linear to other states

    StagnationResult result;

    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);
    double x1 = x(0);
    double x2 = x(1);

    double dx = x2 - x1;
    Vector2d dx_x(-1.0, 1.0);
    double rx = x2 / x1;
    Vector2d rx_x(-rx / x1, 1.0 / x1);

    // Linear extrapolation weights
    double w1 = x2 / dx;
    Vector2d w1_x = -w1 / dx * dx_x + Vector2d(0.0, 1.0) / dx;
    double w2 = -x1 / dx;
    Vector2d w2_x = -w2 / dx * dx_x + Vector2d(-1.0, 0.0) / dx;

    // Stagnation state (linear extrapolation for th, ds, sa)
    result.Ust = U1 * w1 + U2 * w2;

    // Quadratic extrapolation of edge velocity for better slope: ue = K*x
    double wk1 = rx / dx;
    Vector2d wk1_x = rx_x / dx - wk1 / dx * dx_x;
    double wk2 = -1.0 / (rx * dx);
    Vector2d wk2_x = -wk2 * (rx_x / rx + dx_x / dx);

    double K = wk1 * U1(3) + wk2 * U2(3);
    Vector8d K_U = Vector8d::Zero();
    K_U(3) = wk1;
    K_U(7) = wk2;
    Vector2d K_x = U1(3) * wk1_x + U2(3) * wk2_x;

    // Stagnation coord must be small but nonzero
    result.xst = 1e-6;
    result.Ust(3) = K * result.xst;

    // Build linearization w.r.t. U
    result.Ust_U.setZero();
    for (int i = 0; i < 3; ++i) {
        result.Ust_U(i, i) = w1;
        result.Ust_U(i, i + 4) = w2;
    }
    result.Ust_U.row(3) = K_U.transpose() * result.xst;

    // Build linearization w.r.t. x
    result.Ust_x.setZero();
    for (int i = 0; i < 3; ++i) {
        result.Ust_x.row(i) = U1(i) * w1_x + U2(i) * w2_x;
    }
    result.Ust_x.row(3) = K_x * result.xst;

    return result;
}

WakeResult WakeSys(const FoilState& state, Param& param) {
    // Construct residual system for wake initialization
    // Enforces theta_wake = theta_lower + theta_upper
    //         ds_wake = ds_lower + ds_upper + hTE
    //         ctau_wake = weighted average of lower/upper ctau

    WakeResult result;

    int il = state.vsol.Is[0].back();   // Lower surface TE index
    int iu = state.vsol.Is[1].back();   // Upper surface TE index
    int iw = state.vsol.Is[2].front();  // First wake index

    Vector4d Ul = state.glob.U.col(il);
    Vector4d Uu = state.glob.U.col(iu);
    Vector4d Uw = state.glob.U.col(iw);

    // Get TE gap
    double hTE = state.geom.hTE;

    // Get shear stress at wake start
    // Use transition ctau if not already turbulent
    param.turb = true;
    param.wake = false;

    double ctl, ctu;
    Vector4d ctl_Ul = Vector4d::Zero();
    Vector4d ctu_Uu = Vector4d::Zero();

    if (state.vsol.turb[il]) {
        ctl = Ul(2);
        ctl_Ul(2) = 1.0;
    } else {
        auto [ct, ct_U] = GetCttr(Ul, param);
        ctl = ct;
        ctl_Ul = ct_U;
    }

    if (state.vsol.turb[iu]) {
        ctu = Uu(2);
        ctu_Uu(2) = 1.0;
    } else {
        auto [ct, ct_U] = GetCttr(Uu, param);
        ctu = ct;
        ctu_Uu = ct_U;
    }

    // Theta-weighted average of ctau
    double thsum = Ul(0) + Uu(0);
    double ctw = (ctl * Ul(0) + ctu * Uu(0)) / thsum;

    Vector4d ctw_Ul = (ctl_Ul * Ul(0) + (ctl - ctw) * Vector4d(1, 0, 0, 0)) / thsum;
    Vector4d ctw_Uu = (ctu_Uu * Uu(0) + (ctu - ctw) * Vector4d(1, 0, 0, 0)) / thsum;

    // Build residual: wake state = sum of TE states (plus gap for ds)
    result.R(0) = Uw(0) - (Ul(0) + Uu(0));
    result.R(1) = Uw(1) - (Ul(1) + Uu(1) + hTE);
    result.R(2) = Uw(2) - ctw;

    // Store indices
    result.J = {il, iu, iw};

    // Build Jacobian (3x12 for lower, upper, wake 4-states each)
    result.R_U.setZero();

    // Block 0: d/dUl (columns 0-3)
    result.R_U(0, 0) = -1.0;
    result.R_U(1, 1) = -1.0;
    result.R_U.block<1, 4>(2, 0) = -ctw_Ul.transpose();

    // Block 1: d/dUu (columns 4-7)
    result.R_U(0, 4) = -1.0;
    result.R_U(1, 5) = -1.0;
    result.R_U.block<1, 4>(2, 4) = -ctw_Uu.transpose();

    // Block 2: d/dUw (columns 8-11)
    result.R_U(0, 8) = 1.0;
    result.R_U(1, 9) = 1.0;
    result.R_U(2, 10) = 1.0;

    return result;
}

}  // namespace cppfoil
