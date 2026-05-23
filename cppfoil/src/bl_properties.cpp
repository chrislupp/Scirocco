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

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-XXXX.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "cppfoil/properties/bl_properties.h"

#include <algorithm>
#include <cmath>

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

        Vector4d Hk_U = (H_U - 0.29 * M2_U) / den - Hk / den * 0.113 * M2_U;
        return {Hk, Hk_U};
    }

    return {H, H_U};
}

ScalarLin4 GetHw(const Vector4d& U, double wgap) {
    // Wake gap shape parameter Hw = wgap/theta
    // This is added to H = (ds - wgap)/th to recover the full ds/th
    double th = std::max(U(0), 1e-12);

    double Hw = wgap / th;

    Vector4d Hw_U = Vector4d::Zero();
    Hw_U(0) = -Hw / th;  // dHw/dth = -wgap/th^2

    return {Hw, Hw_U};
}

ScalarLin4 GetHs(const Vector4d& U, const Param& param) {
    // Kinetic energy shape parameter (matches mfoil get_Hs)
    auto [Hk, Hk_U] = GetHk(U, param);

    // Limit Hk
    if (param.wake && Hk < 1.00005) {
        Hk = 1.00005;
        Hk_U.setZero();
    }
    if (!param.wake && Hk < 1.05) {
        Hk = 1.05;
        Hk_U.setZero();
    }

    double Hs;
    Vector4d Hs_U;

    if (param.turb) {
        // Turbulent correlation
        constexpr double Hsmin = 1.5;
        constexpr double dHsinf = 0.015;

        auto [Ret, Ret_U] = GetRet(U, param);

        // Limit Re_theta and dependence
        double Ho = 4.0;
        Vector4d Ho_U = Vector4d::Zero();
        if (Ret > 400.0) {
            Ho = 3.0 + 400.0 / Ret;
            Ho_U = -400.0 / (Ret * Ret) * Ret_U;
        }

        double Reb = Ret;
        Vector4d Reb_U = Ret_U;
        if (Ret < 200.0) {
            Reb = 200.0;
            Reb_U.setZero();
        }

        if (Hk < Ho) {
            // Attached branch
            double Hr = (Ho - Hk) / (Ho - 1.0);
            Vector4d Hr_U = (Ho_U - Hk_U) / (Ho - 1.0) -
                            (Ho - Hk) / ((Ho - 1.0) * (Ho - 1.0)) * Ho_U;

            // Match mfoil's exact FP ordering: Hr**2 is computed as a single
            // scalar, then the vector is multiplied by it ONCE — not
            // (vector*Hr)*Hr which has different rounding.
            double Hr2 = Hr * Hr;
            double aa = (2.0 - Hsmin - 4.0 / Reb) * Hr2;
            Vector4d aa_U = (4.0 / (Reb * Reb) * Reb_U) * Hr2 +
                            (2.0 - Hsmin - 4.0 / Reb) * 2.0 * Hr * Hr_U;

            Hs = Hsmin + 4.0 / Reb + aa * 1.5 / (Hk + 0.5);
            Hs_U = -4.0 / (Reb * Reb) * Reb_U + aa_U * 1.5 / (Hk + 0.5) -
                   aa * 1.5 / ((Hk + 0.5) * (Hk + 0.5)) * Hk_U;
        } else {
            // Separated branch
            double lrb = std::log(Reb);
            Vector4d lrb_U = Reb_U / Reb;

            double aa = Hk - Ho + 4.0 / lrb;
            Vector4d aa_U = Hk_U - Ho_U - 4.0 / (lrb * lrb) * lrb_U;

            double bb = 0.007 * lrb / (aa * aa) + dHsinf / Hk;
            Vector4d bb_U = 0.007 * (lrb_U / (aa * aa) -
                                     2.0 * lrb / (aa * aa * aa) * aa_U) -
                            dHsinf / (Hk * Hk) * Hk_U;

            Hs = Hsmin + 4.0 / Reb + (Hk - Ho) * (Hk - Ho) * bb;
            Hs_U = -4.0 / (Reb * Reb) * Reb_U +
                   2.0 * (Hk - Ho) * (Hk_U - Ho_U) * bb +
                   (Hk - Ho) * (Hk - Ho) * bb_U;
        }

        // Mach number correction
        auto [M2, M2_U] = GetMach2(U, param);
        double den = 1.0 + 0.014 * M2;
        Hs = (Hs + 0.028 * M2) / den;
        Hs_U = (Hs_U + 0.028 * M2_U) / den - Hs / den * 0.014 * M2_U;
    } else {
        // Laminar correlation
        double a = Hk - 4.35;
        double Hs_Hk;

        if (Hk < 4.35) {
            double num = 0.0111 * a * a - 0.0278 * a * a * a;
            Hs = num / (Hk + 1.0) + 1.528 - 0.0002 * (a * Hk) * (a * Hk);
            Hs_Hk = (0.0111 * 2.0 * a - 0.0278 * 3.0 * a * a) / (Hk + 1.0) -
                    num / ((Hk + 1.0) * (Hk + 1.0)) -
                    0.0002 * 2.0 * a * Hk * (Hk + a);
        } else {
            Hs = 0.015 * a * a / Hk + 1.528;
            Hs_Hk = 0.015 * 2.0 * a / Hk - 0.015 * a * a / (Hk * Hk);
        }
        Hs_U = Hs_Hk * Hk_U;
    }

    return {Hs, Hs_U};
}

ScalarLin4 GetHss(const Vector4d& U, const Param& param) {
    // Density shape parameter Hss = M2 * (0.064/(Hk-0.8) + 0.251)
    // For incompressible flow (M2=0), Hss = 0
    auto [M2, M2_U] = GetMach2(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);

    double num = 0.064 / (Hk - 0.8) + 0.251;
    Vector4d num_U = -0.064 / ((Hk - 0.8) * (Hk - 0.8)) * Hk_U;

    double Hss = M2 * num;
    Vector4d Hss_U = M2_U * num + M2 * num_U;

    return {Hss, Hss_U};
}

ScalarLin4 GetRet(const Vector4d& U, const Param& param) {
    double th = std::max(U(0), 1e-12);
    double ue = U(3);

    if (param.Minf > 0) {
        auto [M2, M2_U] = GetMach2(U, param);
        auto [uk, uk_U] = GetUk(ue, param);
        double H0 = param.H0;
        double gmi = param.gam - 1.0;
        double Ts = param.Tsrat;
        // Temperature ratio
        double Tr = 1.0 - 0.5 * uk * uk / H0;
        double Tr_uk = -uk / H0;
        // Sutherland's law
        double f = std::pow(Tr, 1.5) * (1.0 + Ts) / (Tr + Ts);
        double f_Tr = 1.5 * f / Tr - f / (Tr + Ts);
        double mu = param.mu0 * f;
        double mu_uk = param.mu0 * f_Tr * Tr_uk;
        // Isentropic density
        double den = 1.0 + 0.5 * gmi * M2;
        double den_M2 = 0.5 * gmi;
        double rho = param.rho0 * std::pow(den, -1.0 / gmi);
        Vector4d rho_U = (-1.0 / gmi) * rho / den * den_M2 * M2_U;
        // Ret
        double Ret = rho * uk * th / mu;
        Vector4d Ret_U = rho_U * uk * th / mu;
        double uk_u = uk_U(3);
        Vector4d due = Vector4d::Zero();
        due(3) = uk_u;
        Ret_U += (rho * th / mu - Ret / mu * mu_uk) * due;
        Ret_U(0) += rho * uk / mu;
        return {Ret, Ret_U};
    }

    double Ret = param.rho0 * ue * th / param.mu0;
    Vector4d Ret_U = Vector4d::Zero();
    Ret_U(0) = param.rho0 * ue / param.mu0;
    Ret_U(3) = param.rho0 * th / param.mu0;
    return {Ret, Ret_U};
}

ScalarLin4 GetCf(const Vector4d& U, const Param& param) {
    // Skin friction coefficient (matches mfoil get_cf)
    if (param.wake) {
        return {0.0, Vector4d::Zero()};
    }

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double cf;
    Vector4d cf_U;

    if (param.turb) {
        // Turbulent cf (mfoil correlation)
        auto [M2, M2_U] = GetMach2(U, param);
        double Fc = std::sqrt(1.0 + 0.5 * (param.gam - 1.0) * M2);
        Vector4d Fc_U = 0.5 / Fc * 0.5 * (param.gam - 1.0) * M2_U;

        double aa = -1.33 * Hk;
        Vector4d aa_U = -1.33 * Hk_U;
        // Smooth limiting of aa
        if (aa < -17.0) {
            double old_aa = aa;
            aa = -20.0 + 3.0 * std::exp((old_aa + 17.0) / 3.0);
            aa_U = (aa + 20.0) / 3.0 * aa_U;
        }

        double bb = std::log(Ret / Fc);
        Vector4d bb_U = Ret_U / Ret - Fc_U / Fc;
        if (bb < 3.0) {
            bb = 3.0;
            bb_U.setZero();
        }
        bb /= std::log(10.0);
        bb_U /= std::log(10.0);

        double cc = -1.74 - 0.31 * Hk;
        Vector4d cc_U = -0.31 * Hk_U;

        double dd = std::tanh(4.0 - Hk / 0.875);
        Vector4d dd_U = (1.0 - dd * dd) * (-Hk_U / 0.875);

        double cf0 = 0.3 * std::exp(aa) * std::pow(bb, cc);
        Vector4d cf0_U =
            cf0 * aa_U +
            0.3 * std::exp(aa) * cc * std::pow(bb, cc - 1.0) * bb_U +
            cf0 * std::log(bb) * cc_U;

        cf = (cf0 + 1.1e-4 * (dd - 1.0)) / Fc;
        cf_U = (cf0_U + 1.1e-4 * dd_U) / Fc - cf / Fc * Fc_U;
    } else {
        // Laminar cf (mfoil correlation)
        double num, num_Hk;
        if (Hk < 5.5) {
            double term = std::pow(5.5 - Hk, 3) / (Hk + 1.0);
            num = 0.0727 * term - 0.07;
            num_Hk =
                0.0727 * (3.0 * std::pow(5.5 - Hk, 2) / (Hk + 1.0) * (-1.0) -
                          std::pow(5.5 - Hk, 3) / std::pow(Hk + 1.0, 2));
        } else {
            double invterm = 1.0 - 1.0 / (Hk - 4.5);
            num = 0.015 * invterm * invterm - 0.07;
            num_Hk = 0.015 * 2.0 * invterm / std::pow(Hk - 4.5, 2);
        }
        cf = num / Ret;
        cf_U = num_Hk / Ret * Hk_U - num / (Ret * Ret) * Ret_U;
    }

    return {cf, cf_U};
}

// Helper: laminar dissipation (used by laminar flow and as floor for turbulent)
ScalarLin4 GetCDi_Lam(const Vector4d& U, const Param& param) {
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double num, num_Hk;
    if (Hk < 4.0) {
        double term = std::pow(4.0 - Hk, 5.5);
        num = 0.00205 * term + 0.207;
        num_Hk = 0.00205 * 5.5 * std::pow(4.0 - Hk, 4.5) * (-1.0);
    } else {
        double Hk1 = Hk - 4.0;
        double denom = 1.0 + 0.02 * Hk1 * Hk1;
        num = -0.0016 * Hk1 * Hk1 / denom + 0.207;
        num_Hk = -0.0016 * (2.0 * Hk1 / denom -
                            Hk1 * Hk1 / (denom * denom) * 0.02 * 2.0 * Hk1);
    }

    double cDi = num / Ret;
    Vector4d cDi_U = num_Hk / Ret * Hk_U - num / (Ret * Ret) * Ret_U;

    return {cDi, cDi_U};
}

// Helper: turbulent wall contribution
ScalarLin4 GetCDi_TurbWall(const Vector4d& U, const Param& param) {
    if (param.wake) return {0.0, Vector4d::Zero()};

    auto [cf, cf_U] = GetCf(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double lr = std::log(Ret);
    Vector4d lr_U = Ret_U / Ret;

    double Hmin = 1.0 + 2.1 / lr;
    Vector4d Hmin_U = -2.1 / (lr * lr) * lr_U;

    double aa = std::tanh((Hk - 1.0) / (Hmin - 1.0));
    double fac = 0.5 + 0.5 * aa;
    Vector4d fac_U = 0.5 * (1.0 - aa * aa) *
                     (Hk_U / (Hmin - 1.0) -
                      (Hk - 1.0) / ((Hmin - 1.0) * (Hmin - 1.0)) * Hmin_U);

    double cDi = 0.5 * cf * Us * (2.0 / Hs) * fac;
    Vector4d cDi_U = cf_U * Us / Hs * fac + cf * Us_U / Hs * fac -
                     cDi / Hs * Hs_U + cf * Us / Hs * fac_U;

    return {cDi, cDi_U};
}

// Helper: outer layer contribution (turbulent)
ScalarLin4 GetCDi_Outer(const Vector4d& U, const Param& param) {
    if (!param.turb) return {0.0, Vector4d::Zero()};

    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);

    // Shear stress: state stores ct^0.5, so ct = sa^2
    double ct = U(2) * U(2);
    Vector4d ct_U = Vector4d::Zero();
    ct_U(2) = 2.0 * U(2);

    double cDi = ct * (0.995 - Us) * 2.0 / Hs;
    Vector4d cDi_U = ct_U * (0.995 - Us) * 2.0 / Hs + ct * (-Us_U) * 2.0 / Hs -
                     ct * (0.995 - Us) * 2.0 / (Hs * Hs) * Hs_U;

    return {cDi, cDi_U};
}

// Helper: laminar stress contribution (turbulent)
ScalarLin4 GetCDi_LamStress(const Vector4d& U, const Param& param) {
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Us, Us_U] = GetUs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double HsRet = Hs * Ret;
    Vector4d HsRet_U = Hs_U * Ret + Hs * Ret_U;

    double num = 0.15 * (0.995 - Us) * (0.995 - Us) * 2.0;
    double num_Us = 0.15 * 2.0 * (0.995 - Us) * (-1.0) * 2.0;

    double cDi = num / HsRet;
    Vector4d cDi_U = num_Us * Us_U / HsRet - num / (HsRet * HsRet) * HsRet_U;

    return {cDi, cDi_U};
}

// Helper: laminar wake dissipation (matches mfoil get_cDi_lamwake)
ScalarLin4 GetCDi_LamWake(const Vector4d& U, const Param& paramin) {
    // Force laminar mode for sub-function calls
    Param param = paramin;
    param.turb = false;

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);

    double HsRet = Hs * Ret;
    Vector4d HsRet_U = Hs_U * Ret + Hs * Ret_U;

    // num = 2*1.1*(1-1/Hk)^2*(1/Hk)
    double Hk_inv = 1.0 / Hk;
    double one_minus_Hk_inv = 1.0 - Hk_inv;
    double num = 2.0 * 1.1 * one_minus_Hk_inv * one_minus_Hk_inv * Hk_inv;
    // d(num)/dHk = 2*1.1 * (2*(1-1/Hk)*(1/Hk^2)*(1/Hk) + (1-1/Hk)^2*(-1/Hk^2))
    double num_Hk = 2.0 * 1.1 *
                    (2.0 * one_minus_Hk_inv * (Hk_inv * Hk_inv) * Hk_inv +
                     one_minus_Hk_inv * one_minus_Hk_inv * (-Hk_inv * Hk_inv));

    double cDi = num / HsRet;
    Vector4d cDi_U = num_Hk * Hk_U / HsRet - num / (HsRet * HsRet) * HsRet_U;

    return {cDi, cDi_U};
}

ScalarLin4 GetCDi(const Vector4d& U, const Param& param) {
    // Dissipation function cDi = 2*cD/H* (matches mfoil get_cDi)
    if (param.turb) {
        double cDi = 0.0;
        Vector4d cDi_U = Vector4d::Zero();

        if (!param.wake) {
            // Turbulent wall contribution
            auto [cDi0, cDi0_U] = GetCDi_TurbWall(U, param);
            cDi += cDi0;
            cDi_U += cDi0_U;
        }

        // Outer layer contribution
        auto [cDi_outer, cDi_outer_U] = GetCDi_Outer(U, param);
        cDi += cDi_outer;
        cDi_U += cDi_outer_U;

        // Laminar stress contribution
        auto [cDi_lamstress, cDi_lamstress_U] = GetCDi_LamStress(U, param);
        cDi += cDi_lamstress;
        cDi_U += cDi_lamstress_U;

        // Maximum check against laminar value
        // (use laminar wake dissipation in wake, normal laminar otherwise)
        auto [cDil, cDil_U] =
            param.wake ? GetCDi_LamWake(U, param) : GetCDi_Lam(U, param);
        if (cDil > cDi) {
            cDi = cDil;
            cDi_U = cDil_U;
        }

        // Double dissipation in wake
        if (param.wake) {
            cDi *= 2.0;
            cDi_U *= 2.0;
        }

        return {cDi, cDi_U};
    } else {
        // Laminar dissipation
        return GetCDi_Lam(U, param);
    }
}

ScalarLin4 GetCpFromUe(double ue, const Param& param) {
    double Vinf = std::max(param.Vinf, 1e-10);
    double cp = 1.0 - (ue / Vinf) * (ue / Vinf);
    double dcp_due = -2.0 * ue / (Vinf * Vinf);

    // Karman-Tsien correction. Match mfoil's get_cp: compute den from
    // cp_inc, then update cp to KT-corrected, then chain dcp/due using
    // (1 - cp_KT * den_cp) / den_inc.
    if (param.KTl > 0) {
        double b = param.KTb;
        double l = param.KTl;
        double den = b + 0.5 * l * (1.0 + b) * cp;
        double den_cp = 0.5 * l * (1.0 + b);
        if (std::abs(den) > 1e-10) {
            cp = cp / den;
            dcp_due *= (1.0 - cp * den_cp) / den;
        }
    }

    Vector4d cp_U = Vector4d::Zero();
    cp_U(3) = dcp_due;

    return {cp, cp_U};
}

ScalarLin4 GetUk(double ue, const Param& param) {
    if (param.Minf > 0) {
        double l = param.KTl;
        double Vinf = param.Vinf;
        double den = 1.0 - l * (ue / Vinf) * (ue / Vinf);
        double den_u = -2.0 * l * ue / (Vinf * Vinf);
        double uk = ue * (1.0 - l) / den;
        double uk_u = (1.0 - l) / den - (uk / den) * den_u;
        Vector4d uk_U = Vector4d::Zero();
        uk_U(3) = uk_u;
        return {uk, uk_U};
    }
    Vector4d uk_U = Vector4d::Zero();
    uk_U(3) = 1.0;
    return {ue, uk_U};
}

ScalarLin4 GetMach2(const Vector4d& U, const Param& param) {
    if (param.Minf > 0) {
        double H0 = param.H0;
        double g = param.gam;
        auto [uk, uk_U] = GetUk(U(3), param);
        double c2 = (g - 1.0) * (H0 - 0.5 * uk * uk);
        double c2_uk = (g - 1.0) * (-uk);
        double M2 = uk * uk / c2;
        double M2_uk = 2.0 * uk / c2 - M2 / c2 * c2_uk;
        Vector4d M2_U = Vector4d::Zero();
        M2_U(3) = M2_uk * uk_U(3);
        return {M2, M2_U};
    }
    return {0.0, Vector4d::Zero()};
}

ScalarLin4 GetRho(const Vector4d& U, const Param& param) {
    // Local density from stagnation density + isentropic relation:
    //   rho = rho0 / (1 + 0.5·(γ-1)·M²)^(1/(γ-1))
    // Matches mfoil's get_rho (param.rho0 IS the stagnation density set by
    // init_thermo as rhoinf·(1+0.5·(γ-1)·Minf²)^(1/(γ-1))).
    double rho = param.rho0;
    Vector4d rho_U = Vector4d::Zero();

    if (param.Minf > 0) {
        auto [M2, M2_U] = GetMach2(U, param);
        double gmi = param.gam - 1.0;
        double den = 1.0 + 0.5 * gmi * M2;
        double den_M2 = 0.5 * gmi;
        rho = param.rho0 * std::pow(den, -1.0 / gmi);
        double drho_dM2 = (-1.0 / gmi) * rho / den * den_M2;
        rho_U = drho_dM2 * M2_U;
    }

    return {rho, rho_U};
}

ScalarLin4 GetDamp(const Vector4d& U, const Param& param) {
    // Amplification rate dn/dx (matches mfoil get_damp)
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Ret, Ret_U] = GetRet(U, param);
    double th = U(0);

    // Limit Hk
    if (Hk < 1.05) {
        Hk = 1.05;
        Hk_U.setZero();
    }

    double Hmi = 1.0 / (Hk - 1.0);
    Vector4d Hmi_U = -Hmi * Hmi * Hk_U;

    double aa = 2.492 * std::pow(Hmi, 0.43);
    Vector4d aa_U = 0.43 * aa / Hmi * Hmi_U;

    double bb = std::tanh(14.0 * Hmi - 9.24);
    Vector4d bb_U = (1.0 - bb * bb) * 14.0 * Hmi_U;

    double lrc = aa + 0.7 * (bb + 1.0);
    Vector4d lrc_U = aa_U + 0.7 * bb_U;

    constexpr double lten = 2.302585092994046;  // log(10)
    double lr = std::log(Ret) / lten;
    Vector4d lr_U = (1.0 / Ret) * Ret_U / lten;

    double damp = 0.0;
    Vector4d damp_U = Vector4d::Zero();

    constexpr double dl = 0.1;
    if (lr >= lrc - dl) {
        double rn = (lr - (lrc - dl)) / (2.0 * dl);
        Vector4d rn_U = (lr_U - lrc_U) / (2.0 * dl);

        double rf, rf_U_coeff;
        if (rn >= 1.0) {
            rf = 1.0;
            rf_U_coeff = 0.0;
        } else {
            rf = 3.0 * rn * rn - 2.0 * rn * rn * rn;
            rf_U_coeff = 6.0 * rn - 6.0 * rn * rn;
        }
        Vector4d rf_U = rf_U_coeff * rn_U;

        double ar = 3.87 * Hmi - 2.52;
        Vector4d ar_U = 3.87 * Hmi_U;

        double ex = std::exp(-ar * ar);
        Vector4d ex_U = ex * (-2.0 * ar * ar_U);

        double da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
        Vector4d da_U = 0.028 * Hk_U - 0.0345 * ex_U;

        double af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi +
                    3.0 * Hmi * Hmi * Hmi + 0.1 * std::exp(-20.0 * Hmi);
        Vector4d af_U =
            (2.7 - 11.0 * Hmi + 9.0 * Hmi * Hmi - 2.0 * std::exp(-20.0 * Hmi)) *
            Hmi_U;

        damp = rf * af * da / th;
        damp_U = (rf_U * af * da + rf * af_U * da + rf * af * da_U) / th;
        damp_U(0) -= damp / th;
    }

    // Extra amplification to ensure dn/dx > 0 near ncrit
    double ncrit = param.ncrit;
    constexpr double Cea = 5.0;
    double nx = Cea * (U(2) - ncrit);
    Vector4d nx_U = Vector4d::Zero();
    nx_U(2) = Cea;

    double eex = 1.0 + std::tanh(nx);
    Vector4d eex_U = (1.0 - std::tanh(nx) * std::tanh(nx)) * nx_U;

    double ed = eex * 0.001 / th;
    Vector4d ed_U = eex_U * 0.001 / th;
    ed_U(0) -= ed / th;

    damp += ed;
    damp_U += ed_U;

    return {damp, damp_U};
}

ScalarLin4 GetCteq(const Vector4d& U, const Param& param) {
    // Equilibrium sqrt(ctau) (matches mfoil get_cteq)
    double CC = 0.5 / (param.GA * param.GA * param.GB);
    double C = param.GC;

    auto [Hk, Hk_U] = GetHk(U, param);
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [H, H_U] = GetH(U);
    auto [Ret, Ret_U] = GetRet(U, param);
    auto [Us, Us_U] = GetUs(U, param);

    // Limit Hk and compute Hkc
    double Hkc;
    Vector4d Hkc_U;
    if (param.wake) {
        if (Hk < 1.00005) {
            Hk = 1.00005;
            Hk_U.setZero();
        }
        Hkc = Hk - 1.0;
        Hkc_U = Hk_U;
    } else {
        if (Hk < 1.05) {
            Hk = 1.05;
            Hk_U.setZero();
        }
        Hkc = Hk - 1.0 - C / Ret;
        Hkc_U = Hk_U + C / (Ret * Ret) * Ret_U;
        if (Hkc < 0.01) {
            Hkc = 0.01;
            Hkc_U.setZero();
        }
    }

    // Match mfoil's exact FP ordering: precompute Hkc**2 and Hk**2 as
    // single scalars so vector * scalar^2 is one multiplication, not two
    // sequential vector*scalar steps with different rounding.
    double Hkc2 = Hkc * Hkc;
    double Hk2 = Hk * Hk;
    double num = CC * Hs * (Hk - 1.0) * Hkc2;
    Vector4d num_U = CC * (Hs_U * (Hk - 1.0) * Hkc2 + Hs * Hk_U * Hkc2 +
                           Hs * (Hk - 1.0) * 2.0 * Hkc * Hkc_U);

    double den = (1.0 - Us) * H * Hk2;
    Vector4d den_U = (-Us_U) * H * Hk2 + (1.0 - Us) * H_U * Hk2 +
                     (1.0 - Us) * H * 2.0 * Hk * Hk_U;

    double cteq = std::sqrt(num / den);
    Vector4d cteq_U = 0.5 / cteq * (num_U / den - num / (den * den) * den_U);

    return {cteq, cteq_U};
}

ScalarLin4 GetCttr(const Vector4d& U, const Param& param) {
    // Initial sqrt(ctau) at transition (matches mfoil get_cttr)
    // Note: transition happens just before the wake starts, so force wake=false
    Param param_nowake = param;
    param_nowake.wake = false;

    auto [cteq, cteq_U] = GetCteq(U, param_nowake);
    auto [Hk, Hk_U] = GetHk(U, param_nowake);

    // Limit Hk
    if (Hk < 1.05) {
        Hk = 1.05;
        Hk_U.setZero();
    }

    double C = param.CtauC;
    double E = param.CtauE;

    double c = C * std::exp(-E / (Hk - 1.0));
    Vector4d c_U = c * E / ((Hk - 1.0) * (Hk - 1.0)) * Hk_U;

    double cttr = c * cteq;
    Vector4d cttr_U = c_U * cteq + c * cteq_U;

    return {cttr, cttr_U};
}

ScalarLin8 GetUpw(const Vector4d& U1, const Vector4d& U2, const Param& param) {
    // Upwind factor for stability (matches mfoil get_upw)
    // Decision to upwind is made based on the shape factor change
    auto [Hk1, Hk1_U] = GetHk(U1, param);
    auto [Hk2, Hk2_U] = GetHk(U2, param);

    Vector4d Z = Vector4d::Zero();

    // Triggering constant and wake factor
    double Hut = 1.0;
    double C = param.wake ? 1.0 : 5.0;
    double Huc = C * Hut / (Hk2 * Hk2);  // only depends on U2

    Vector8d Huc_U;
    Huc_U.head<4>() = Z;
    Huc_U.tail<4>() = -2.0 * Huc / Hk2 * Hk2_U;

    // aa = (Hk2-1)/(Hk1-1), log-based upwinding
    double aa = (Hk2 - 1.0) / (Hk1 - 1.0);
    double sga = (aa >= 0) ? 1.0 : -1.0;
    double la = std::log(sga * aa);

    Vector8d la_U;
    la_U.head<4>() = -1.0 / (Hk1 - 1.0) * Hk1_U;
    la_U.tail<4>() = 1.0 / (Hk2 - 1.0) * Hk2_U;

    double Hls = la * la;
    Vector8d Hls_U = 2.0 * la * la_U;

    // Cap at 15
    if (Hls > 15.0) {
        Hls = 15.0;
        Hls_U.setZero();
    }

    double upw = 1.0 - 0.5 * std::exp(-Hls * Huc);
    Vector8d upw_U = -0.5 * std::exp(-Hls * Huc) * (-Hls_U * Huc - Hls * Huc_U);

    return {upw, upw_U};
}

ScalarLin8 Upwind(double upw, const Vector8d& upw_U, double f1,
                  const Vector4d& f1_U1, double f2, const Vector4d& f2_U2) {
    double f = (1.0 - upw) * f1 + upw * f2;

    // Match mfoil's exact FP ordering:
    //   f_U = (-upw_U)*f1 + upw_U*f2 + concat((1-upw)*f1_U1, upw*f2_U2)
    // Refactoring to `upw_U*(f2 - f1)` is mathematically equivalent but
    // produces different floating-point rounding when f1 ≈ f2 (typical
    // for upwinded BL closures), which the ForwardRx reference test
    // detects at the 1e-9 rel level.
    Vector8d f_U = (-upw_U) * f1 + upw_U * f2;
    f_U.head<4>() += (1.0 - upw) * f1_U1;
    f_U.tail<4>() += upw * f2_U2;

    return {f, f_U};
}

ScalarLin4 GetDe(const Vector4d& U, const Param& param) {
    // BL thickness measure: delta* + weighted theta
    auto [Hk, Hk_U] = GetHk(U, param);

    double th = std::max(U(0), 1e-12);
    double ds = U(1);

    // Weight depends on Hk
    double aa = 3.15 + 1.72 / (Hk - 1.0);
    Vector4d aa_U = -1.72 / ((Hk - 1.0) * (Hk - 1.0)) * Hk_U;

    double de = th * aa + ds;
    Vector4d de_U;
    de_U << aa, 1.0, 0.0, 0.0;
    de_U += th * aa_U;

    // Cap at 12*theta
    double dmx = 12.0;
    if (de > dmx * th) {
        de = dmx * th;
        de_U << dmx, 0.0, 0.0, 0.0;
    }

    return {de, de_U};
}

ScalarLin4 GetUs(const Vector4d& U, const Param& param) {
    // Normalized wall slip velocity
    auto [Hs, Hs_U] = GetHs(U, param);
    auto [Hk, Hk_U] = GetHk(U, param);
    auto [H, H_U] = GetH(U);

    // Limit Hk
    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) {
        Hk = Hk_lim;
        Hk_U.setZero();
    }

    double beta = param.GB;
    double bi = 1.0 / beta;

    double Us = 0.5 * Hs * (1.0 - bi * (Hk - 1.0) / H);
    Vector4d Us_U =
        0.5 * Hs_U * (1.0 - bi * (Hk - 1.0) / H) +
        0.5 * Hs * (-bi * Hk_U / H + bi * (Hk - 1.0) / (H * H) * H_U);

    // Limits
    if (!param.wake) {
        if (Us > 0.95) {
            Us = 0.98;
            Us_U.setZero();
        }
    }

    return {Us, Us_U};
}

ScalarLinX GetCfxt(const Vector4d& U, double x, const Param& param) {
    // cf * x / theta
    auto [cf, cf_U] = GetCf(U, param);
    double th = std::max(U(0), 1e-12);

    double cfxt = cf * x / th;
    Vector4d cfxt_U = cf_U * x / th;
    cfxt_U(0) -= cfxt / th;
    double cfxt_x = cf / th;

    return {cfxt, cfxt_U, cfxt_x};
}

ScalarLinX GetCDixt(const Vector4d& U, double x, const Param& param) {
    // cDi * x / theta (dissipation function)
    auto [cDi, cDi_U] = GetCDi(U, param);
    double th = std::max(U(0), 1e-12);

    double cDixt = cDi * x / th;
    Vector4d cDixt_U = cDi_U * x / th;
    cDixt_U(0) -= cDixt / th;
    double cDixt_x = cDi / th;

    return {cDixt, cDixt_U, cDixt_x};
}

ScalarLin8 GetUq(double ds, const Vector8d& ds_U, double cf,
                 const Vector8d& cf_U, double Hk, const Vector8d& Hk_U,
                 double Ret, const Vector8d& Ret_U, const Param& param) {
    // Equilibrium pressure gradient 1/ue * due/dx
    double beta = param.GB;
    double A = param.GA;
    double C = param.GC;

    if (param.wake) {
        A *= param.Dlr;
        C = 0.0;
    }

    // Limit Hk
    double Hk_lim = param.wake ? 1.00005 : 1.05;
    Vector8d Hk_U_use = Hk_U;
    if (Hk < Hk_lim) {
        Hk = Hk_lim;
        Hk_U_use.setZero();
    }

    double Hkc = Hk - 1.0 - C / std::max(Ret, 1e-10);
    Vector8d Hkc_U = Hk_U_use + C / (Ret * Ret) * Ret_U;

    if (Hkc < 0.01) {
        Hkc = 0.01;
        Hkc_U.setZero();
    }

    // Match mfoil ordering: compute scalars once (avoid `A*Hk*Hk` left-to-right
    // which has different FP rounding than `A * (Hk**2)`).
    double AHk = A * Hk;
    double Hk2 = Hk * Hk;
    double AHk2 = A * Hk2;
    double ut = 0.5 * cf - (Hkc / AHk) * (Hkc / AHk);
    Vector8d ut_U =
        0.5 * cf_U - 2.0 * (Hkc / AHk) * (Hkc_U / AHk - Hkc / AHk2 * Hk_U_use);

    double ds_safe = std::max(ds, 1e-12);
    double uq = ut / (beta * ds_safe);
    Vector8d uq_U = ut_U / (beta * ds_safe) - uq / ds_safe * ds_U;

    return {uq, uq_U};
}

}  // namespace cppfoil
