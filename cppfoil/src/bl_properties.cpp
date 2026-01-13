// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/properties/bl_properties.h"

#include <cmath>
#include <algorithm>

namespace cppfoil {

ScalarLin4 GetH(const Vector4d& U) {
    // H = delta*/theta
    double th = std::max(U(0), 1e-12);
    double ds = U(1);
    double H = ds / th;

    Vector4d H_U = Vector4d::Zero();
    H_U(0) = -ds / (th * th);  // dH/dth
    H_U(1) = 1.0 / th;         // dH/dds

    return {H, H_U};
}

ScalarLin4 GetHk(const Vector4d& U, const Param& param) {
    auto [H, H_U] = GetH(U);

    if (param.Minf > 0) {
        auto [M2, M2_U] = GetMach2(U, param);
        double den = 1.0 + 0.113 * M2;
        double Hk = (H - 0.29 * M2) / den;

        Vector4d Hk_U = (H_U - 0.29 * M2_U) / den -
                        Hk / den * 0.113 * M2_U;
        return {Hk, Hk_U};
    }

    return {H, H_U};
}

ScalarLin4 GetHw(const Vector4d& U, double wgap) {
    auto [H, H_U] = GetH(U);
    double th = std::max(U(0), 1e-12);

    // Wake gap correction
    double Hw = H - wgap / th;

    Vector4d Hw_U = H_U;
    Hw_U(0) += wgap / (th * th);  // Additional term from gap

    return {Hw, Hw_U};
}

ScalarLin4 GetHs(const Vector4d& U, const Param& param) {
    auto [Hk, Hk_U] = GetHk(U, param);

    // Kinetic energy shape parameter correlation
    double Hs;
    Vector4d Hs_U;

    if (param.turb) {
        // Turbulent correlation
        double Hk0 = std::max(Hk, 1.05);
        if (Hk < 5.5) {
            Hs = 0.028 * (Hk0 - 1.0) + 0.014 / std::pow(Hk0 - 1.0, 2);
        } else {
            Hs = 0.028 * (Hk0 - 1.0) + 0.014 / 20.25;
        }
        double dHs_dHk = 0.028 - 0.028 / std::pow(std::max(Hk0 - 1.0, 0.05), 3);
        Hs_U = dHs_dHk * Hk_U;
    } else {
        // Laminar correlation
        Hs = 1.515 + 0.076 * (Hk - 4.0) * (Hk - 4.0) / Hk;
        double dHs_dHk = 0.076 * (2.0 * (Hk - 4.0) / Hk - (Hk - 4.0) * (Hk - 4.0) / (Hk * Hk));
        Hs_U = dHs_dHk * Hk_U;
    }

    return {Hs, Hs_U};
}

ScalarLin4 GetHss(const Vector4d& U, const Param& param) {
    // Density shape parameter (simplified)
    auto [Hk, Hk_U] = GetHk(U, param);
    double Hss = 0.064 / (Hk - 0.8) + 0.251;
    double dHss_dHk = -0.064 / ((Hk - 0.8) * (Hk - 0.8));

    return {Hss, dHss_dHk * Hk_U};
}

ScalarLin4 GetRet(const Vector4d& U, const Param& param) {
    double th = std::max(U(0), 1e-12);
    double ue = U(3);

    double Ret = param.rho0 * std::abs(ue) * th / param.mu0;

    Vector4d Ret_U = Vector4d::Zero();
    Ret_U(0) = param.rho0 * std::abs(ue) / param.mu0;
    Ret_U(3) = param.rho0 * th * (ue > 0 ? 1.0 : -1.0) / param.mu0;

    return {Ret, Ret_U};
}

ScalarLin4 GetCf(const Vector4d& U, const Param& param) {
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double cf;
    Vector4d cf_U;

    if (param.turb) {
        // Turbulent: Ludwieg-Tillmann correlation (simplified)
        double Fc = 1.0;  // Compressibility factor
        cf = Fc * 0.246 / std::pow(std::max(Ret, 1.0), 0.268) *
             std::pow(10.0, -0.678 * Hk);

        double dcf_dRet = -0.268 * cf / std::max(Ret, 1.0);
        double dcf_dHk = -0.678 * std::log(10.0) * cf;

        cf_U = dcf_dRet * Ret_U + dcf_dHk * Hk_U;
    } else {
        // Laminar: Thwaites correlation (simplified)
        double lam = (Hk - 1.0) / (Hk - 1.0 + 1.0);  // Simplified lambda
        cf = 2.0 * (0.22 + 1.402 * lam - 0.018 * lam * lam / (0.107 + lam)) /
             std::max(Ret, 1.0);

        // Simplified linearization
        double dcf_dRet = -cf / std::max(Ret, 1.0);
        cf_U = dcf_dRet * Ret_U;
    }

    if (param.wake) {
        cf = 0.0;
        cf_U.setZero();
    }

    return {cf, cf_U};
}

ScalarLin4 GetCDi(const Vector4d& U, const Param& param) {
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [cf, cf_U] = GetCf(U, param);

    // Dissipation coefficient
    double cDi = 0.5 * cf * Hs;
    Vector4d cDi_U = 0.5 * (cf_U * Hs + cf * Hs_U);

    return {cDi, cDi_U};
}

ScalarLin4 GetCpFromUe(double ue, const Param& param) {
    double Vinf = std::max(param.Vinf, 1e-10);
    double cp = 1.0 - (ue / Vinf) * (ue / Vinf);

    // Karman-Tsien correction
    if (param.KTl > 0) {
        double denom = param.KTb + param.KTl * (1.0 + param.KTb) * 0.5 * cp;
        if (std::abs(denom) > 1e-10) {
            cp = cp / denom;
        }
    }

    double dcp_due = -2.0 * ue / (Vinf * Vinf);

    Vector4d cp_U = Vector4d::Zero();
    cp_U(3) = dcp_due;

    return {cp, cp_U};
}

ScalarLin4 GetUk(double ue, const Param& param) {
    // Compressible velocity correction (simplified)
    double uk = ue;
    Vector4d uk_U = Vector4d::Zero();
    uk_U(3) = 1.0;

    return {uk, uk_U};
}

ScalarLin4 GetMach2(const Vector4d& U, const Param& param) {
    double ue = U(3);
    double a_inf = param.Vinf / std::max(param.Minf, 1e-10);

    double M2 = (ue / a_inf) * (ue / a_inf);

    Vector4d M2_U = Vector4d::Zero();
    M2_U(3) = 2.0 * ue / (a_inf * a_inf);

    return {M2, M2_U};
}

ScalarLin4 GetRho(const Vector4d& U, const Param& param) {
    // Density from isentropic relations (simplified)
    double rho = param.rho0;
    Vector4d rho_U = Vector4d::Zero();

    if (param.Minf > 0) {
        auto [M2, M2_U] = GetMach2(U, param);
        double gam = param.gam;
        double fac = 1.0 + 0.5 * (gam - 1.0) * param.Minf * param.Minf;
        double fac_local = 1.0 + 0.5 * (gam - 1.0) * M2;
        rho = param.rho0 * std::pow(fac / fac_local, 1.0 / (gam - 1.0));

        double drho_dM2 = -0.5 * rho / fac_local;
        rho_U = drho_dM2 * M2_U;
    }

    return {rho, rho_U};
}

ScalarLin4 GetDamp(const Vector4d& U, const Param& param) {
    // Amplification factor rate (eN method)
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    // Gaster criterion (simplified)
    double Hk_crit = 2.8;
    double damp = 0.0;
    Vector4d damp_U = Vector4d::Zero();

    if (Hk < Hk_crit && Ret > 0) {
        double f = (Hk_crit - Hk) / (Hk_crit - 1.0);
        damp = 0.01 * f * f * std::sqrt(Ret);

        double ddamp_dHk = -0.02 * f / (Hk_crit - 1.0) * std::sqrt(Ret);
        double ddamp_dRet = 0.005 * f * f / std::sqrt(std::max(Ret, 1e-10));

        damp_U = ddamp_dHk * Hk_U + ddamp_dRet * Ret_U;
    }

    return {damp, damp_U};
}

ScalarLin4 GetCttr(const Vector4d& U, const Param& param) {
    // Initial ctau at transition
    auto [Hk, Hk_U] = GetHk(U, param);

    double cttr = param.CtauC * std::pow(Hk, -param.CtauE);
    double dcttr_dHk = -param.CtauE * cttr / Hk;

    return {cttr, dcttr_dHk * Hk_U};
}

ScalarLin4 GetCteq(const Vector4d& U, const Param& param) {
    // Equilibrium ctau
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [cf, cf_U] = GetCf(U, param);

    double Hk0 = std::max(Hk, 1.05);
    double cteq = 0.03 * std::sqrt(0.5 * cf) * (Hk0 - 1.0);

    Vector4d cteq_U = Vector4d::Zero();
    if (cf > 1e-12) {
        double dcteq_dcf = 0.015 / std::sqrt(2.0 * cf) * (Hk0 - 1.0);
        double dcteq_dHk = 0.03 * std::sqrt(0.5 * cf);
        cteq_U = dcteq_dcf * cf_U + dcteq_dHk * Hk_U;
    }

    return {cteq, cteq_U};
}

ScalarLin8 GetUpw(const Vector4d& U1, const Vector4d& U2, const Param& param) {
    // Upwind factor for stability
    auto [Hk1, Hk1_U] = GetHk(U1, param);
    auto [Hk2, Hk2_U] = GetHk(U2, param);

    // Simple upwind based on shape parameter
    double upw = 0.5;
    Vector8d upw_U = Vector8d::Zero();

    if (Hk1 > 1.05 && Hk2 > 1.05) {
        double dH = Hk2 - Hk1;
        if (std::abs(dH) > 0.1) {
            upw = (dH > 0) ? 0.7 : 0.3;
            // Simplified linearization
            upw_U.head<4>() = -0.4 * Hk1_U / std::abs(dH);
            upw_U.tail<4>() = 0.4 * Hk2_U / std::abs(dH);
        }
    }

    return {upw, upw_U};
}

ScalarLin8 Upwind(double upw, const Vector8d& upw_U,
                   double f1, const Vector4d& f1_U1,
                   double f2, const Vector4d& f2_U2) {
    double f = (1.0 - upw) * f1 + upw * f2;

    Vector8d f_U;
    f_U.head<4>() = (1.0 - upw) * f1_U1;
    f_U.tail<4>() = upw * f2_U2;
    f_U += upw_U * (f2 - f1);

    return {f, f_U};
}

}  // namespace cppfoil
