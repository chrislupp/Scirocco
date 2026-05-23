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
#include "cppfoil/core/boundary_layer.h"

#include <cmath>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/properties/bl_properties.h"

namespace cppfoil {

std::pair<double, double> ThwaitesInit(double K, double nu) {
    // Thwaites correlation for stagnation point
    // K = due/dx at stagnation (ue = K*x)
    // nu = kinematic viscosity
    //
    // th^2 = ue^(-6) * 0.45 * nu * int_0^x ue^5 dx = 0.45*nu/(6*K)
    // ds = Hstag*th = 2.2*th

    double th = std::sqrt(0.45 * nu / (6.0 * std::max(K, 1e-10)));
    double ds = 2.2 * th;

    return {th, ds};
}

Param BuildParam(const FoilState& state, int si) {
    Param param = state.param;

    // Set surface-specific flags
    param.wake = (si == 2);
    param.turb = param.wake;  // Wake is fully turbulent (matches mfoil)
    param.simi = false;

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

    // Similarity flag: true at stagnation point (matches mfoil)
    param.simi = (i == state.isol.Istag[0] || i == state.isol.Istag[1]);
}

void InitBoundaryLayer(FoilState& state) {
    // Initialize boundary layer by marching from stagnation point
    // Matches mfoil's init_boundary_layer approach

    constexpr double Hmaxl = 3.8;  // Max H before laminar separation
    constexpr double Hmaxt = 2.5;  // Max H before turbulent separation
    constexpr int nNewton = 30;
    constexpr int iNswitch =
        12;  // Switch to inverse after this many iterations

    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    // Get inviscid edge velocity
    VectorXd ueinv = GetUeInv(state);

    state.glob.Nsys = Nsys;
    state.glob.U.resize(4, Nsys);
    state.glob.U.setZero();
    state.vsol.turb.resize(Nsys, false);

    // Kinematic viscosity
    double nu = state.param.mu0 / state.param.rho0;

    // Loop over surfaces: 0=lower, 1=upper, 2=wake
    for (int si = 0; si < 3; ++si) {
        const auto& Is = state.vsol.Is[si];
        int Nsurf = static_cast<int>(Is.size());
        if (Nsurf == 0) continue;

        // Local arrays for this surface
        VectorXd xi(Nsurf);
        VectorXd ue(Nsurf);
        VectorXd Aux(Nsurf);
        MatrixXd U(4, Nsurf);
        U.setZero();
        Aux.setZero();

        // Fill local arrays
        for (int i = 0; i < Nsurf; ++i) {
            int idx = Is[i];
            xi(i) = state.isol.xi(idx);
            ue(i) = ueinv(idx);
            if (si == 2 && i < static_cast<int>(state.vsol.wgap.size())) {
                Aux(i) = state.vsol.wgap(i);
            }
        }

        // Ensure edge velocities are not tiny
        double uemax = ue.cwiseAbs().maxCoeff();
        for (int i = 0; i < Nsurf; ++i) {
            ue(i) = std::max(ue(i), 1e-8 * uemax);
        }

        // Build parameter structure for this surface
        Param param = BuildParam(state, si);

        // Initialize first point
        int i0 = 0;

        if (si < 2) {
            // Airfoil surface: stagnation point initialization

            // Stagnation point flow constant K = ue/x
            double K;
            bool hitstag;
            if (xi(0) < 1e-8 * xi(Nsurf - 1)) {
                K = ue(1) / std::max(xi(1), 1e-12);
                hitstag = true;
            } else {
                K = ue(0) / std::max(xi(0), 1e-12);
                hitstag = false;
            }
            // Thwaites initialization
            auto [th, ds] = ThwaitesInit(K, nu);
            double xst = 1e-6;  // Small but nonzero
            Vector4d Ust(th, ds, 0.0, K * xst);

            // Newton refinement at stagnation (matching mfoil)
            param.turb = false;
            param.simi = true;  // Similarity station flag

            for (int iNewton = 0; iNewton < 20; ++iNewton) {
                // Call residual at stagnation with simi=true
                Vector2d xpair(xst, xst);
                Matrix42d Upair;
                Upair.col(0) = Ust;
                Upair.col(1) = Ust;
                Vector2d Auxpair(0.0, 0.0);

                auto res = ResidualStation(param, xpair, Upair, Auxpair);

                if (res.value.norm() < 1e-10) break;

                // Solve: A = R_U[:,4:7] + R_U[:,0:3] (both sides same state)
                Matrix3d A =
                    res.lin_U.block<3, 3>(0, 4) + res.lin_U.block<3, 3>(0, 0);
                Vector3d b = -res.value;
                Vector3d dU3 = A.partialPivLu().solve(b);

                // Under-relaxation
                double dm = std::max(std::abs(dU3(0) / Ust(0)),
                                     std::abs(dU3(1) / Ust(1)));
                double omega = (dm < 0.2) ? 1.0 : 0.2 / dm;
                dU3 *= omega;

                Ust(0) += dU3(0);
                Ust(1) += dU3(1);
                Ust(2) += dU3(2);
            }

            param.simi = false;  // Reset similarity flag

            // Store stagnation state (mfoil approach)
            // When hitstag, both point 0 and point 1 get the stagnation state,
            // and marching starts from point 2
            if (hitstag) {
                U.col(0) = Ust;
                U(3, 0) = ue(0);
                U.col(1) = Ust;
                U(3, 1) = ue(1);
                i0 = 1;  // March will start from i0+1 = 2
            } else {
                U.col(0) = Ust;
                U(3, 0) = ue(0);
                i0 = 0;  // March will start from i0+1 = 1
            }

        } else {
            // Wake: initialize from TE states
            i0 = 0;  // Reset for wake
            auto wres = WakeSys(state, param);
            // WakeSys modifies param.wake to false for cttr calculation;
            // reset it to true for wake marching (matches mfoil behavior)
            param.wake = true;
            Vector4d Uw = state.glob.U.col(Is[0]);
            Uw.head<3>() -= wres.R;
            Uw(3) = ue(0);
            U.col(0) = Uw;

            param.turb = true;
            state.vsol.turb[Is[0]] = true;
        }

        // March over rest of points
        bool tran = false;  // Transition flag
        int i = i0 + 1;

        while (i < Nsurf) {
            // Initialize with previous state and new ue
            U.col(i) = U.col(i - 1);
            U(3, i) = ue(i);

            // Set ctau at transition (but don't set param.turb yet - matches
            // mfoil)
            if (tran) {
                auto [ct, ct_U] = GetCttr(U.col(i), param);
                U(2, i) = ct;
            }

            state.vsol.turb[Is[i]] = (tran || param.turb);
            bool direct = true;
            double Hktgt = 0.0;

            // Newton iteration at this station
            bool converged = false;
            for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
                Vector2d xpair(xi(i - 1), xi(i));
                Matrix42d Upair;
                Upair.col(0) = U.col(i - 1);
                Upair.col(1) = U.col(i);
                Vector2d Auxpair(Aux(i - 1), Aux(i));

                // Use transition residual at transition stations
                ResidualLin res;
                if (tran) {
                    res =
                        ResidualTransition(state, param, xpair, Upair, Auxpair);
                } else {
                    res = ResidualStation(param, xpair, Upair, Auxpair);
                }

                if (res.value.norm() < 1e-10) {
                    converged = true;
                    break;
                }

                Vector4d dU;
                if (direct) {
                    // Direct mode: solve for th, ds, sa (ue prescribed)
                    Matrix3d A = res.lin_U.block<3, 3>(0, 4);
                    dU.head<3>() = A.partialPivLu().solve(-res.value);
                    dU(3) = 0.0;
                } else {
                    // Inverse mode: prescribe Hk, solve for th, ds, sa, ue
                    auto [Hk, Hk_U] = GetHk(U.col(i), param);

                    Eigen::Matrix<double, 4, 4> A;
                    A.topRows<3>() = res.lin_U.block<3, 4>(0, 4);
                    A.row(3) = Hk_U.transpose();

                    Vector4d b;
                    b.head<3>() = -res.value;
                    b(3) = Hktgt - Hk;

                    dU = A.partialPivLu().solve(b);
                }

                // Under-relaxation
                double dm =
                    std::max(std::abs(dU(0) / std::max(U(0, i - 1), 1e-12)),
                             std::abs(dU(1) / std::max(U(1, i - 1), 1e-12)));
                if (!direct)
                    dm = std::max(
                        dm, std::abs(dU(3) / std::max(U(3, i - 1), 1e-10)));
                if (param.turb)
                    dm = std::max(
                        dm, std::abs(dU(2) / std::max(U(2, i - 1), 1e-10)));
                else if (direct)
                    dm = std::max(dm, std::abs(dU(2) / 10.0));

                double omega = (dm > 0.3) ? 0.3 / dm : 1.0;

                dU *= omega;

                // Trial update
                Vector4d Ui = U.col(i) + dU;

                // Clip extreme values
                Ui(0) = std::max(Ui(0), 1e-12);  // theta > 0
                Ui(1) = std::max(Ui(1), 1e-12);  // delta* > 0
                if (param.turb) {
                    Ui(2) = std::max(std::min(Ui(2), 0.3), 1e-7);
                }

                // Check if about to separate
                double Hmax = param.turb ? Hmaxt : Hmaxl;
                auto [Hk, Hk_U] = GetHk(Ui, param);

                if (direct && ((Hk > Hmax) || (iNewton > iNswitch))) {
                    // Switch to inverse mode
                    direct = false;
                    auto [Hkprev, Hkprev_U] = GetHk(U.col(i - 1), param);
                    double Hkr =
                        (xi(i) - xi(i - 1)) / std::max(U(0, i - 1), 1e-12);

                    if (param.wake) {
                        double H2 = Hkprev;
                        for (int k = 0; k < 6; ++k) {
                            H2 -= (H2 + 0.03 * Hkr * std::pow(H2 - 1, 3) -
                                   Hkprev) /
                                  (1.0 + 0.09 * Hkr * std::pow(H2 - 1, 2));
                        }
                        Hktgt = std::max(H2, 1.01);
                    } else if (param.turb) {
                        Hktgt = Hkprev - 0.15 * Hkr;
                    } else {
                        Hktgt = Hkprev + 0.03 * Hkr;
                    }
                    if (!param.wake) {
                        Hktgt = std::max(Hktgt, Hmax);
                    }
                    if (iNewton > iNswitch) {
                        // Reinitialize
                        U.col(i) = U.col(i - 1);
                        U(3, i) = ue(i);
                    }
                } else {
                    U.col(i) = Ui;
                }
            }

            // If not converged, extrapolate or use fallback
            if (!converged) {
                U.col(i) = U.col(i - 1);
                U(3, i) = ue(i);

                // At transition station, preserve cttr instead of laminar amp
                if (tran) {
                    auto [ct, ct_U] = GetCttr(U.col(i), param);
                    U(2, i) = ct;
                }
                // Extrapolate theta and delta* (matches mfoil: si<3 is always
                // true)
                double xratio = std::sqrt(xi(i) / xi(i - 1));
                U(0, i) = U(0, i - 1) * xratio;
                U(1, i) = U(1, i - 1) * xratio;
            }

            // Check for transition (laminar only)
            if (!param.turb && !tran && U(2, i) > param.ncrit) {
                tran = true;
                continue;  // Redo station with transition
            }

            if (tran) {
                // Store transition location and set turbulent flag (matches
                // mfoil)
                state.vsol.Xt(si, 0) = 0.5 * (xi(i - 1) + xi(i));
                param.turb = true;
                tran = false;
            }

            ++i;
        }

        // Store states in global array
        for (int j = 0; j < Nsurf; ++j) {
            state.glob.U.col(Is[j]) = U.col(j);
        }
    }
}

// ============================================================================
// ∂closure/∂Ret and ∂closure/∂Ma helpers at fixed U. Used by
// ComputeStationLinRet / ComputeStationLinMa to build the per-station
// lin_Ret / lin_Ma for the analytic Re / Ma chains. Wrapped in detail
// namespace (not anonymous) so they can be probed from tests.
//
// Each helper inspects the closure's existing formula in bl_properties.cpp
// and extracts the explicit Ret- or Ma-derivative path.
// ============================================================================
namespace detail {

double GetCf_dRet(const Vector4d& U, const Param& param) {
    if (param.wake) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    if (param.turb) {
        auto [M2, _M2_U] = GetMach2(U, param);
        (void)_M2_U;
        double Fc = std::sqrt(1.0 + 0.5 * (param.gam - 1.0) * M2);
        double aa = -1.33 * Hk;
        if (aa < -17.0) aa = -20.0 + 3.0 * std::exp((aa + 17.0) / 3.0);
        double bb_nat = std::log(Ret / Fc);
        double bb_nat_dRet = (bb_nat < 3.0) ? 0.0 : 1.0 / Ret;
        double bb = (bb_nat < 3.0) ? 3.0 : bb_nat;
        bb /= std::log(10.0);
        double bb_dRet = bb_nat_dRet / std::log(10.0);
        double cc = -1.74 - 0.31 * Hk;
        double dd = std::tanh(4.0 - Hk / 0.875);
        double cf0 = 0.3 * std::exp(aa) * std::pow(bb, cc);
        double cf0_dRet = (bb > 0) ? 0.3 * std::exp(aa) * cc *
                                         std::pow(bb, cc - 1.0) * bb_dRet
                                   : 0.0;
        (void)dd;
        (void)cf0;
        return cf0_dRet / Fc;
    } else {
        double num;
        if (Hk < 5.5) {
            num = 0.0727 * std::pow(5.5 - Hk, 3) / (Hk + 1.0) - 0.07;
        } else {
            double invterm = 1.0 - 1.0 / (Hk - 4.5);
            num = 0.015 * invterm * invterm - 0.07;
        }
        return -num / (Ret * Ret);
    }
}

double GetHs_dRet(const Vector4d& U, const Param& param) {
    if (!param.turb) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    if (param.wake && Hk < 1.00005) Hk = 1.00005;
    if (!param.wake && Hk < 1.05) Hk = 1.05;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    constexpr double Hsmin = 1.5;
    constexpr double dHsinf = 0.015;
    double Ho = 4.0;
    double Ho_dRet = 0.0;
    if (Ret > 400.0) {
        Ho = 3.0 + 400.0 / Ret;
        Ho_dRet = -400.0 / (Ret * Ret);
    }
    double Reb = Ret;
    double Reb_dRet = 1.0;
    if (Ret < 200.0) {
        Reb = 200.0;
        Reb_dRet = 0.0;
    }
    double Hs_dRet;
    if (Hk < Ho) {
        double Hr = (Ho - Hk) / (Ho - 1.0);
        double Hr_dRet = Ho_dRet / (Ho - 1.0) -
                         (Ho - Hk) / ((Ho - 1.0) * (Ho - 1.0)) * Ho_dRet;
        double Hr2 = Hr * Hr;
        double aa_dRet = (4.0 / (Reb * Reb) * Reb_dRet) * Hr2 +
                         (2.0 - Hsmin - 4.0 / Reb) * 2.0 * Hr * Hr_dRet;
        Hs_dRet = -4.0 / (Reb * Reb) * Reb_dRet + aa_dRet * 1.5 / (Hk + 0.5);
    } else {
        double lrb = std::log(Reb);
        double lrb_dRet = Reb_dRet / Reb;
        double aa = Hk - Ho + 4.0 / lrb;
        double aa_dRet = -Ho_dRet - 4.0 / (lrb * lrb) * lrb_dRet;
        double bb_dRet = 0.007 * (lrb_dRet / (aa * aa) -
                                  2.0 * lrb / (aa * aa * aa) * aa_dRet);
        Hs_dRet = -4.0 / (Reb * Reb) * Reb_dRet +
                  2.0 * (Hk - Ho) * (-Ho_dRet) *
                      (0.007 * lrb / (aa * aa) + dHsinf / Hk) +
                  (Hk - Ho) * (Hk - Ho) * bb_dRet;
    }
    auto [M2, _M2_U] = GetMach2(U, param);
    (void)_M2_U;
    double den = 1.0 + 0.014 * M2;
    return Hs_dRet / den;
}

double GetUs_dRet(const Vector4d& U, const Param& param) {
    // Us = 0.5·Hs·(1 - bi·(Hk-1)/H) with bi = 1/param.GB. Limits: if
    // !wake and Us > 0.95 then Us = 0.98 (zero derivative).
    // ∂Us/∂Ret at fixed U: Hk and H have no Ret dep, so only Hs chains.
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [H, _H_U] = GetH(U);
    (void)_H_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) Hk = Hk_lim;
    double bi = 1.0 / param.GB;
    double factor = 1.0 - bi * (Hk - 1.0) / H;
    double Us = 0.5 * Hs * factor;
    if (!param.wake && Us > 0.95) return 0.0;
    return 0.5 * GetHs_dRet(U, param) * factor;
}

double GetCDi_Lam_dRet(const Vector4d& U, const Param& param) {
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double num;
    if (Hk < 4.0) {
        num = 0.00205 * std::pow(4.0 - Hk, 5.5) + 0.207;
    } else {
        double Hk1 = Hk - 4.0;
        double denom = 1.0 + 0.02 * Hk1 * Hk1;
        num = -0.0016 * Hk1 * Hk1 / denom + 0.207;
    }
    return -num / (Ret * Ret);
}

double GetCDi_TurbWall_dRet(const Vector4d& U, const Param& param) {
    if (param.wake) return 0.0;
    auto [cf, _cf_U] = GetCf(U, param);
    (void)_cf_U;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double cf_dRet = GetCf_dRet(U, param);
    double Hs_dRet = GetHs_dRet(U, param);
    double Us_dRet = GetUs_dRet(U, param);
    double lr = std::log(Ret);
    double lr_dRet = 1.0 / Ret;
    double Hmin = 1.0 + 2.1 / lr;
    double Hmin_dRet = -2.1 / (lr * lr) * lr_dRet;
    double aa = std::tanh((Hk - 1.0) / (Hmin - 1.0));
    double fac_dRet = 0.5 * (1.0 - aa * aa) *
                      (-(Hk - 1.0) / ((Hmin - 1.0) * (Hmin - 1.0)) * Hmin_dRet);
    double fac = 0.5 + 0.5 * aa;
    double cDi = 0.5 * cf * Us * (2.0 / Hs) * fac;
    return cf_dRet * Us / Hs * fac + cf * Us_dRet / Hs * fac -
           cDi / Hs * Hs_dRet + cf * Us / Hs * fac_dRet;
}

double GetCDi_Outer_dRet(const Vector4d& U, const Param& param) {
    if (!param.turb) return 0.0;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    double Hs_dRet = GetHs_dRet(U, param);
    double Us_dRet = GetUs_dRet(U, param);
    double ct = U(2) * U(2);
    return ct * (-Us_dRet) * 2.0 / Hs -
           ct * (0.995 - Us) * 2.0 / (Hs * Hs) * Hs_dRet;
}

double GetCDi_LamStress_dRet(const Vector4d& U, const Param& param) {
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double Hs_dRet = GetHs_dRet(U, param);
    double Us_dRet = GetUs_dRet(U, param);
    double HsRet = Hs * Ret;
    double HsRet_dRet = Hs_dRet * Ret + Hs;
    double num = 0.15 * (0.995 - Us) * (0.995 - Us) * 2.0;
    double num_dRet = 0.15 * 2.0 * (0.995 - Us) * (-Us_dRet) * 2.0;
    return num_dRet / HsRet - num / (HsRet * HsRet) * HsRet_dRet;
}

double GetCDi_LamWake_dRet(const Vector4d& U, const Param& paramin) {
    Param param = paramin;
    param.turb = false;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    // Hs_dRet in laminar = 0 (Hs has no Ret dep in laminar branch)
    double HsRet = Hs * Ret;
    double HsRet_dRet = Hs;  // Hs_dRet = 0 for laminar
    double Hk_inv = 1.0 / Hk;
    double one_m = 1.0 - Hk_inv;
    double num = 2.0 * 1.1 * one_m * one_m * Hk_inv;
    return -num / (HsRet * HsRet) * HsRet_dRet;
}

double GetCDi_dRet(const Vector4d& U, const Param& param) {
    if (param.turb) {
        double cDi_dRet = 0.0;
        double cDi_val = 0.0;
        if (!param.wake) {
            auto [v, _vU] = GetCDi_TurbWall(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dRet += GetCDi_TurbWall_dRet(U, param);
        }
        {
            auto [v, _vU] = GetCDi_Outer(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dRet += GetCDi_Outer_dRet(U, param);
        }
        {
            auto [v, _vU] = GetCDi_LamStress(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dRet += GetCDi_LamStress_dRet(U, param);
        }
        auto [cDil, _cDil_U] =
            param.wake ? GetCDi_LamWake(U, param) : GetCDi_Lam(U, param);
        (void)_cDil_U;
        double cDil_dRet = param.wake ? GetCDi_LamWake_dRet(U, param)
                                      : GetCDi_Lam_dRet(U, param);
        if (cDil > cDi_val) cDi_dRet = cDil_dRet;
        if (param.wake) cDi_dRet *= 2.0;
        return cDi_dRet;
    }
    return GetCDi_Lam_dRet(U, param);
}

double GetCfxt_dRet(const Vector4d& U, double x, const Param& param) {
    double th = std::max(U(0), 1e-12);
    return GetCf_dRet(U, param) * x / th;
}

double GetCDixt_dRet(const Vector4d& U, double x, const Param& param) {
    double th = std::max(U(0), 1e-12);
    return GetCDi_dRet(U, param) * x / th;
}

// GetCteq, GetUq, GetDamp _dRet — for the turbulent and transition cases.
// cteq depends on Hs, Us, Hk, H, Ret directly. Uq takes Ret as a direct
// arg (not via U). Damp uses Ret in the laminar amplification formula.
// For Mach=0 NACA 0012 alpha=3 Re=1e6 ViscousRe test, these matter only
// at turbulent foil stations and wake. For now keep the chain there.
//
// GetCteq_dRet inspection: from bl_properties.cpp, cteq² has a complex
// formula. Inspect later if needed.
double GetCteq_dRet(const Vector4d& U, const Param& param) {
    if (!param.turb) return 0.0;
    // cteq = sqrt(num/den) where
    //   num = CC * Hs * (Hk-1) * Hkc²
    //   den = (1-Us) * H * Hk²
    //   Hkc = (param.wake) ? Hk-1 : Hk - 1 - C/Ret  (clamped)
    // ∂cteq/∂Ret at fixed U: Hk and H are Ret-indep. Hs, Us, Hkc carry Ret-dep.
    double CC = 0.5 / (param.GA * param.GA * param.GB);
    double C = param.GC;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [H, _H_U] = GetH(U);
    (void)_H_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    double Hs_dRet = GetHs_dRet(U, param);
    double Us_dRet = GetUs_dRet(U, param);

    double Hkc;
    double Hkc_dRet;
    if (param.wake) {
        if (Hk < 1.00005) Hk = 1.00005;
        Hkc = Hk - 1.0;
        Hkc_dRet = 0.0;  // no Ret dep in wake branch
    } else {
        if (Hk < 1.05) Hk = 1.05;
        Hkc = Hk - 1.0 - C / Ret;
        Hkc_dRet = C / (Ret * Ret);
        if (Hkc < 0.01) {
            Hkc = 0.01;
            Hkc_dRet = 0.0;
        }
    }
    double Hkc2 = Hkc * Hkc;
    double Hk2 = Hk * Hk;
    double num = CC * Hs * (Hk - 1.0) * Hkc2;
    double num_dRet = CC * (Hs_dRet * (Hk - 1.0) * Hkc2 +
                            Hs * (Hk - 1.0) * 2.0 * Hkc * Hkc_dRet);
    double den = (1.0 - Us) * H * Hk2;
    double den_dRet = (-Us_dRet) * H * Hk2;
    double cteq = std::sqrt(num / den);
    return 0.5 / cteq * (num_dRet / den - num / (den * den) * den_dRet);
}

// ∂uq/∂Ret_arg at fixed (ds, cf, Hk). GetUq's direct Ret argument enters
// only via Hkc = Hk - 1 - C/Ret (foil non-wake) and the C/Ret²·Ret_U term.
double GetUq_dRet_arg(double cf, double Hk, double Ret, double ds,
                      const Param& param) {
    double beta = param.GB;
    double A = param.GA;
    double C = param.GC;
    if (param.wake) {
        A *= param.Dlr;
        C = 0.0;
    }
    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) Hk = Hk_lim;
    double Hkc = Hk - 1.0 - C / std::max(Ret, 1e-10);
    double Hkc_dRet = C / (Ret * Ret);
    if (Hkc < 0.01) {
        Hkc = 0.01;
        Hkc_dRet = 0.0;
    }
    double AHk = A * Hk;
    double ut_dRet = -2.0 * (Hkc / AHk) * (Hkc_dRet / AHk);
    double ds_safe = std::max(ds, 1e-12);
    return ut_dRet / (beta * ds_safe);
}

// ∂uq/∂cf at fixed (ds, Hk, Ret). ut = 0.5*cf - (Hkc/AHk)². ∂ut/∂cf = 0.5.
double GetUq_dcf(double ds, const Param& param) {
    double ds_safe = std::max(ds, 1e-12);
    return 0.5 / (param.GB * ds_safe);
}

double GetDamp_dRet(const Vector4d& U, const Param& param) {
    // GetDamp uses Ret only via lr = log(Ret)/log(10). Chain:
    //   lr_dRet = 1/(Ret · log(10))
    //   rn = (lr - (lrc - dl)) / (2*dl), lrc has no Ret dep → rn_dRet =
    //   lr_dRet/(2*dl) rf (cubic blend) has Ret dep through rn:
    //     rn < 1: rf = 3rn² - 2rn³ → ∂rf/∂rn = 6rn - 6rn²
    //     rn ≥ 1: rf = 1 (flat) → ∂rf/∂rn = 0
    //   af, da, ar, ex have no Ret dep.
    //   damp = rf · af · da / th → damp_dRet = rf_dRet · af · da / th
    //   ed term: tanh of sa-only, no Ret dep.
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double th = U(0);
    if (Hk < 1.05) Hk = 1.05;

    double Hmi = 1.0 / (Hk - 1.0);
    double aa = 2.492 * std::pow(Hmi, 0.43);
    double bb = std::tanh(14.0 * Hmi - 9.24);
    double lrc = aa + 0.7 * (bb + 1.0);

    constexpr double lten = 2.302585092994046;
    double lr = std::log(Ret) / lten;
    constexpr double dl = 0.1;
    if (lr < lrc - dl) return 0.0;

    double rn = (lr - (lrc - dl)) / (2.0 * dl);
    double lr_dRet = 1.0 / (Ret * lten);
    double rn_dRet = lr_dRet / (2.0 * dl);

    double rf_dRet;
    if (rn >= 1.0) {
        rf_dRet = 0.0;
    } else {
        rf_dRet = (6.0 * rn - 6.0 * rn * rn) * rn_dRet;
    }

    double ar = 3.87 * Hmi - 2.52;
    double ex = std::exp(-ar * ar);
    double da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
    double af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi + 3.0 * Hmi * Hmi * Hmi +
                0.1 * std::exp(-20.0 * Hmi);
    return rf_dRet * af * da / th;
}

// ============================================================================
// ∂(thermo params)/∂Ma — closed-form derivatives of param.{KTb,KTl,H0,mu0,rho0}
// w.r.t. Ma, used by the analytic Ma chain. Mirrors InitThermo's formulas
// (cppfoil/src/inviscid.cpp:42-95).
//
// Returns zero struct for Ma <= 0 (incompressible case).
// ============================================================================
struct ThermoMaDeriv {
    double KTb_dMa = 0;
    double KTl_dMa = 0;
    double H0_dMa = 0;
    double mu0_dMa = 0;
    double rho0_dMa = 0;
};

ThermoMaDeriv ComputeThermoMaDeriv(const Param& param, double Vinf,
                                   double rhoinf) {
    ThermoMaDeriv d;
    double Ma = param.Minf;
    if (Ma <= 0) return d;
    double gam = param.gam;
    double gmi = gam - 1.0;
    double Tsrat = param.Tsrat;

    double KTb = param.KTb;
    // KTb = sqrt(1 - Ma²)  → KTb_dMa = -Ma/KTb
    d.KTb_dMa = -Ma / KTb;

    // KTl = Ma²/(1+KTb)²  → KTl_dMa = 2*Ma/(1+KTb)² - 2*Ma²/(1+KTb)³ · KTb_dMa
    double opb = 1.0 + KTb;
    d.KTl_dMa =
        2.0 * Ma / (opb * opb) - 2.0 * Ma * Ma / (opb * opb * opb) * d.KTb_dMa;

    // H0 = (1 + 0.5*gmi*Ma²)*Vinf²/(gmi*Ma²) = Vinf²/(gmi*Ma²) + 0.5*Vinf²
    //   → H0_dMa = -2·Vinf²/(gmi·Ma³)
    d.H0_dMa = -2.0 * Vinf * Vinf / (gmi * Ma * Ma * Ma);

    // mu0 = muinf/finf, finf = Tr^1.5 · (1+Ts)/(Tr+Ts), Tr = 1 - 0.5·Vinf²/H0
    double H0 = param.H0;
    double muinf = param.muinf;
    double Tr = 1.0 - 0.5 * Vinf * Vinf / H0;
    double Tr_dMa = 0.5 * Vinf * Vinf / (H0 * H0) * d.H0_dMa;
    double finf = std::pow(Tr, 1.5) * (1.0 + Tsrat) / (Tr + Tsrat);
    double finf_dTr = 1.5 * finf / Tr - finf / (Tr + Tsrat);
    double finf_dMa = finf_dTr * Tr_dMa;
    // mu0 = muinf / finf (muinf is Ma-indep)
    d.mu0_dMa = -muinf / (finf * finf) * finf_dMa;

    // rho0 = rhoinf · (1 + 0.5·gmi·Ma²)^(1/gmi)
    //   → rho0_dMa = Ma·rho0/(1 + 0.5·gmi·Ma²)
    double fac = 1.0 + 0.5 * gmi * Ma * Ma;
    d.rho0_dMa = Ma * param.rho0 / fac;
    (void)rhoinf;
    return d;
}

// ∂uk/∂Ma at fixed ue. From GetUk formula with KTl(Ma).
double GetUk_dMa(double ue, const Param& param, const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    double l = param.KTl;
    double Vinf = param.Vinf;
    double x = ue / Vinf;
    double den = 1.0 - l * x * x;
    // ∂uk/∂l with l = KTl:
    //   uk = ue·(1-l)/den, ∂uk/∂l = ue·(x²-1)/den²
    double duk_dKTl = ue * (x * x - 1.0) / (den * den);
    return duk_dKTl * d.KTl_dMa;
}

// ∂M2/∂Ma at fixed U.
double GetMach2_dMa(const Vector4d& U, const Param& param,
                    const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    double gam = param.gam;
    double H0 = param.H0;
    auto [uk, _uk_U] = GetUk(U(3), param);
    (void)_uk_U;
    double uk_dMa = GetUk_dMa(U(3), param, d);
    double c2 = (gam - 1.0) * (H0 - 0.5 * uk * uk);
    double c2_dMa = (gam - 1.0) * (d.H0_dMa - uk * uk_dMa);
    double M2 = uk * uk / c2;
    return (2.0 * uk * uk_dMa - M2 * c2_dMa) / c2;
}

// ∂Hk/∂Ma at fixed U. Hk = (H - 0.29·M2)/(1 + 0.113·M2). H is Ma-indep.
double GetHk_dMa(const Vector4d& U, const Param& param,
                 const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [M2, _M2_U] = GetMach2(U, param);
    (void)_M2_U;
    auto [H, _H_U] = GetH(U);
    (void)_H_U;
    double M2_dMa = GetMach2_dMa(U, param, d);
    double den = 1.0 + 0.113 * M2;
    double Hk = (H - 0.29 * M2) / den;
    return (-0.29 * M2_dMa) / den - Hk / den * 0.113 * M2_dMa;
}

// ∂Hss/∂Ma at fixed U. Hss = M² · (0.064/(Hk-0.8) + 0.251).
double GetHss_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [M2, _M2_U] = GetMach2(U, param);
    (void)_M2_U;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    double M2_dMa = GetMach2_dMa(U, param, d);
    double Hk_dMa = GetHk_dMa(U, param, d);
    double num = 0.064 / (Hk - 0.8) + 0.251;
    double num_dMa = -0.064 / ((Hk - 0.8) * (Hk - 0.8)) * Hk_dMa;
    return M2_dMa * num + M2 * num_dMa;
}

// ∂Ret/∂Ma at fixed U. Uses param.{Minf, H0, mu0, rho0, KTl, gam}.
double GetRet_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    double th = std::max(U(0), 1e-12);
    double H0 = param.H0;
    double gmi = param.gam - 1.0;
    double Tsrat = param.Tsrat;
    auto [M2, _M2_U] = GetMach2(U, param);
    (void)_M2_U;
    auto [uk, _uk_U] = GetUk(U(3), param);
    (void)_uk_U;
    double M2_dMa = GetMach2_dMa(U, param, d);
    double uk_dMa = GetUk_dMa(U(3), param, d);

    double Tr = 1.0 - 0.5 * uk * uk / H0;
    double Tr_dMa_atU =
        -uk * uk_dMa / H0 + 0.5 * uk * uk / (H0 * H0) * d.H0_dMa;
    double f = std::pow(Tr, 1.5) * (1.0 + Tsrat) / (Tr + Tsrat);
    double f_dTr = 1.5 * f / Tr - f / (Tr + Tsrat);
    double f_dMa_atU = f_dTr * Tr_dMa_atU;
    double mu = param.mu0 * f;
    // ∂mu/∂Ma = mu0_dMa · f + mu0 · f_dMa
    double mu_dMa = d.mu0_dMa * f + param.mu0 * f_dMa_atU;

    double den_iso = 1.0 + 0.5 * gmi * M2;
    double rho = param.rho0 * std::pow(den_iso, -1.0 / gmi);
    // rho = rho0 · den_iso^(-1/gmi).
    // ∂rho/∂Ma = rho0_dMa · den_iso^(-1/gmi) + rho0 · (-1/gmi) ·
    // den_iso^(-1/gmi-1) · 0.5·gmi·M2_dMa
    double rho_dMa = d.rho0_dMa * std::pow(den_iso, -1.0 / gmi) +
                     rho * (-0.5 * M2_dMa / den_iso);

    double Ret = rho * uk * th / mu;
    // ∂Ret/∂Ma at fixed U: rho_dMa·uk·th/mu + rho·uk_dMa·th/mu - Ret·mu_dMa/mu
    return rho_dMa * uk * th / mu + rho * uk_dMa * th / mu - Ret * mu_dMa / mu;
}

// ∂Hs/∂Ma at fixed U. Chains through Hk(Ma), M2(Ma), and (turb only) Ret(Ma).
double GetHs_dMa(const Vector4d& U, const Param& param,
                 const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    double Hk_dMa = GetHk_dMa(U, param, d);
    if (param.wake && Hk < 1.00005) {
        Hk = 1.00005;
        Hk_dMa = 0;
    }
    if (!param.wake && Hk < 1.05) {
        Hk = 1.05;
        Hk_dMa = 0;
    }

    double Hs_base_dMa = 0;  // ∂Hs(no Mach correction)/∂Ma via Hk and Ret
    if (param.turb) {
        constexpr double Hsmin = 1.5;
        constexpr double dHsinf = 0.015;
        auto [Ret, _Ret_U] = GetRet(U, param);
        (void)_Ret_U;
        double Ret_dMa = GetRet_dMa(U, param, d);
        double Ho = 4.0, Ho_dRet = 0.0;
        if (Ret > 400.0) {
            Ho = 3.0 + 400.0 / Ret;
            Ho_dRet = -400.0 / (Ret * Ret);
        }
        double Ho_dMa_part = Ho_dRet * Ret_dMa;
        double Reb = Ret, Reb_dRet = 1.0;
        if (Ret < 200.0) {
            Reb = 200.0;
            Reb_dRet = 0.0;
        }
        double Reb_dMa = Reb_dRet * Ret_dMa;
        if (Hk < Ho) {
            double Hr = (Ho - Hk) / (Ho - 1.0);
            double Hr_dMa = (Ho_dMa_part - Hk_dMa) / (Ho - 1.0) -
                            (Ho - Hk) / ((Ho - 1.0) * (Ho - 1.0)) * Ho_dMa_part;
            double Hr2 = Hr * Hr;
            double aa = (2.0 - Hsmin - 4.0 / Reb) * Hr2;
            double aa_dMa = (4.0 / (Reb * Reb) * Reb_dMa) * Hr2 +
                            (2.0 - Hsmin - 4.0 / Reb) * 2.0 * Hr * Hr_dMa;
            Hs_base_dMa = -4.0 / (Reb * Reb) * Reb_dMa +
                          aa_dMa * 1.5 / (Hk + 0.5) -
                          aa * 1.5 / ((Hk + 0.5) * (Hk + 0.5)) * Hk_dMa;
        } else {
            double lrb = std::log(Reb);
            double lrb_dMa = Reb_dMa / Reb;
            double aa = Hk - Ho + 4.0 / lrb;
            double aa_dMa = Hk_dMa - Ho_dMa_part - 4.0 / (lrb * lrb) * lrb_dMa;
            double bb = 0.007 * lrb / (aa * aa) + dHsinf / Hk;
            double bb_dMa = 0.007 * (lrb_dMa / (aa * aa) -
                                     2.0 * lrb / (aa * aa * aa) * aa_dMa) -
                            dHsinf / (Hk * Hk) * Hk_dMa;
            Hs_base_dMa = -4.0 / (Reb * Reb) * Reb_dMa +
                          2.0 * (Hk - Ho) * (Hk_dMa - Ho_dMa_part) * bb +
                          (Hk - Ho) * (Hk - Ho) * bb_dMa;
        }
    } else {
        // Laminar: Hs depends on Hk only (no Mach correction in GetHs's
        // laminar branch — verified by reading bl_properties.cpp:159-175).
        // ∂Hs/∂Ma = Hs_Hk · Hk_dMa, no further correction.
        double a = Hk - 4.35;
        double Hs_Hk;
        if (Hk < 4.35) {
            double num = 0.0111 * a * a - 0.0278 * a * a * a;
            Hs_Hk = (0.0111 * 2.0 * a - 0.0278 * 3.0 * a * a) / (Hk + 1.0) -
                    num / ((Hk + 1.0) * (Hk + 1.0)) -
                    0.0002 * 2.0 * a * Hk * (Hk + a);
        } else {
            Hs_Hk = 0.015 * 2.0 * a / Hk - 0.015 * a * a / (Hk * Hk);
        }
        return Hs_Hk * Hk_dMa;
    }
    // Turbulent branch: GetHs applies the Mach correction
    //   Hs = (Hs_base + 0.028·M²) / (1 + 0.014·M²)
    // ∂Hs/∂Ma = (Hs_base_dMa + 0.028·M2_dMa)/den - Hs_corrected/den ·
    // 0.014·M2_dMa
    auto [M2, _M2_U] = GetMach2(U, param);
    (void)_M2_U;
    double M2_dMa = GetMach2_dMa(U, param, d);
    auto [Hs_corrected, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    double den = 1.0 + 0.014 * M2;
    return (Hs_base_dMa + 0.028 * M2_dMa) / den -
           Hs_corrected / den * 0.014 * M2_dMa;
}

double GetUs_dMa(const Vector4d& U, const Param& param,
                 const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    // Us = 0.5·Hs·(1 - bi·(Hk-1)/H), bi = 1/param.GB.
    // ∂Us/∂Ma at fixed U: Hs(Ma) AND Hk(Ma); H has no Ma dep.
    //   ∂Us/∂Ma = 0.5·Hs_dMa·factor + 0.5·Hs·factor_dMa
    //   factor = 1 - bi·(Hk-1)/H  → factor_dMa = -bi·Hk_dMa/H
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [H, _H_U] = GetH(U);
    (void)_H_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    double Hk_dMa = GetHk_dMa(U, param, d);
    double Hs_dMa = GetHs_dMa(U, param, d);
    double Hk_lim = param.wake ? 1.00005 : 1.05;
    if (Hk < Hk_lim) {
        Hk = Hk_lim;
        Hk_dMa = 0.0;
    }
    double bi = 1.0 / param.GB;
    double factor = 1.0 - bi * (Hk - 1.0) / H;
    double factor_dMa = -bi * Hk_dMa / H;
    double Us = 0.5 * Hs * factor;
    if (!param.wake && Us > 0.95) return 0.0;
    return 0.5 * Hs_dMa * factor + 0.5 * Hs * factor_dMa;
}

double GetCf_dMa(const Vector4d& U, const Param& param,
                 const ThermoMaDeriv& d) {
    if (param.wake || param.Minf <= 0) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double Hk_dMa = GetHk_dMa(U, param, d);
    double Ret_dMa = GetRet_dMa(U, param, d);
    if (param.turb) {
        auto [M2, _M2_U] = GetMach2(U, param);
        (void)_M2_U;
        double M2_dMa = GetMach2_dMa(U, param, d);
        double Fc = std::sqrt(1.0 + 0.5 * (param.gam - 1.0) * M2);
        double Fc_dMa = 0.5 / Fc * 0.5 * (param.gam - 1.0) * M2_dMa;
        double aa = -1.33 * Hk;
        double aa_dMa = -1.33 * Hk_dMa;
        if (aa < -17.0) {
            double old_aa = aa;
            aa = -20.0 + 3.0 * std::exp((old_aa + 17.0) / 3.0);
            aa_dMa = (aa + 20.0) / 3.0 * aa_dMa;
        }
        double bb_nat = std::log(Ret / Fc);
        double bb_nat_dMa = Ret_dMa / Ret - Fc_dMa / Fc;
        double bb = bb_nat;
        double bb_dMa = bb_nat_dMa;
        if (bb < 3.0) {
            bb = 3.0;
            bb_dMa = 0;
        }
        bb /= std::log(10.0);
        bb_dMa /= std::log(10.0);
        double cc = -1.74 - 0.31 * Hk;
        double cc_dMa = -0.31 * Hk_dMa;
        double dd = std::tanh(4.0 - Hk / 0.875);
        double dd_dMa = (1.0 - dd * dd) * (-Hk_dMa / 0.875);
        double cf0 = 0.3 * std::exp(aa) * std::pow(bb, cc);
        double cf0_dMa =
            cf0 * aa_dMa +
            0.3 * std::exp(aa) * cc * std::pow(bb, cc - 1.0) * bb_dMa +
            cf0 * std::log(std::max(bb, 1e-15)) * cc_dMa;
        double cf = (cf0 + 1.1e-4 * (dd - 1.0)) / Fc;
        return (cf0_dMa + 1.1e-4 * dd_dMa) / Fc - cf / Fc * Fc_dMa;
    } else {
        // Laminar: cf = num/Ret, num depends only on Hk → Ma chain through
        // Hk_dMa and Ret_dMa
        double num;
        double num_Hk;
        if (Hk < 5.5) {
            num = 0.0727 * std::pow(5.5 - Hk, 3) / (Hk + 1.0) - 0.07;
            num_Hk =
                0.0727 * (3.0 * std::pow(5.5 - Hk, 2) / (Hk + 1.0) * (-1.0) -
                          std::pow(5.5 - Hk, 3) / ((Hk + 1.0) * (Hk + 1.0)));
        } else {
            double invterm = 1.0 - 1.0 / (Hk - 4.5);
            num = 0.015 * invterm * invterm - 0.07;
            num_Hk = 0.015 * 2.0 * invterm / ((Hk - 4.5) * (Hk - 4.5));
        }
        return num_Hk * Hk_dMa / Ret - num / (Ret * Ret) * Ret_dMa;
    }
}

// CDi sub-functions — chain rule each. For laminar (cDi=num/Ret) only Ret_dMa
// matters; for turb the chain goes through cf, Hk, Hs, Us.
double GetCDi_Lam_dMa(const Vector4d& U, const Param& param,
                      const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double Hk_dMa = GetHk_dMa(U, param, d);
    double Ret_dMa = GetRet_dMa(U, param, d);
    double num, num_Hk;
    if (Hk < 4.0) {
        num = 0.00205 * std::pow(4.0 - Hk, 5.5) + 0.207;
        num_Hk = 0.00205 * 5.5 * std::pow(4.0 - Hk, 4.5) * (-1.0);
    } else {
        double Hk1 = Hk - 4.0;
        double denom = 1.0 + 0.02 * Hk1 * Hk1;
        num = -0.0016 * Hk1 * Hk1 / denom + 0.207;
        num_Hk = -0.0016 * (2.0 * Hk1 / denom -
                            Hk1 * Hk1 / (denom * denom) * 0.02 * 2.0 * Hk1);
    }
    return num_Hk * Hk_dMa / Ret - num / (Ret * Ret) * Ret_dMa;
}

double GetCDi_TurbWall_dMa(const Vector4d& U, const Param& param,
                           const ThermoMaDeriv& d) {
    if (param.wake || param.Minf <= 0) return 0.0;
    auto [cf, _cf_U] = GetCf(U, param);
    (void)_cf_U;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double cf_dMa = GetCf_dMa(U, param, d);
    double Hs_dMa = GetHs_dMa(U, param, d);
    double Us_dMa = GetUs_dMa(U, param, d);
    double Hk_dMa = GetHk_dMa(U, param, d);
    double Ret_dMa = GetRet_dMa(U, param, d);

    double lr = std::log(Ret);
    double lr_dMa = Ret_dMa / Ret;
    double Hmin = 1.0 + 2.1 / lr;
    double Hmin_dMa = -2.1 / (lr * lr) * lr_dMa;
    double aa = std::tanh((Hk - 1.0) / (Hmin - 1.0));
    double fac = 0.5 + 0.5 * aa;
    double fac_dMa = 0.5 * (1.0 - aa * aa) *
                     (Hk_dMa / (Hmin - 1.0) -
                      (Hk - 1.0) / ((Hmin - 1.0) * (Hmin - 1.0)) * Hmin_dMa);
    double cDi = 0.5 * cf * Us * (2.0 / Hs) * fac;
    return cf_dMa * Us / Hs * fac + cf * Us_dMa / Hs * fac - cDi / Hs * Hs_dMa +
           cf * Us / Hs * fac_dMa;
}

double GetCDi_Outer_dMa(const Vector4d& U, const Param& param,
                        const ThermoMaDeriv& d) {
    if (!param.turb || param.Minf <= 0) return 0.0;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    double Hs_dMa = GetHs_dMa(U, param, d);
    double Us_dMa = GetUs_dMa(U, param, d);
    double ct = U(2) * U(2);
    return ct * (-Us_dMa) * 2.0 / Hs -
           ct * (0.995 - Us) * 2.0 / (Hs * Hs) * Hs_dMa;
}

double GetCDi_LamStress_dMa(const Vector4d& U, const Param& param,
                            const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Us, _Us_U] = GetUs(U, param);
    (void)_Us_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double Hs_dMa = GetHs_dMa(U, param, d);
    double Us_dMa = GetUs_dMa(U, param, d);
    double Ret_dMa = GetRet_dMa(U, param, d);
    double HsRet = Hs * Ret;
    double HsRet_dMa = Hs_dMa * Ret + Hs * Ret_dMa;
    double num = 0.15 * (0.995 - Us) * (0.995 - Us) * 2.0;
    double num_dMa = 0.15 * 2.0 * (0.995 - Us) * (-Us_dMa) * 2.0;
    return num_dMa / HsRet - num / (HsRet * HsRet) * HsRet_dMa;
}

double GetCDi_LamWake_dMa(const Vector4d& U, const Param& paramin,
                          const ThermoMaDeriv& d) {
    if (paramin.Minf <= 0) return 0.0;
    Param param = paramin;
    param.turb = false;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    auto [Hs, _Hs_U] = GetHs(U, param);
    (void)_Hs_U;
    auto [Ret, _Ret_U] = GetRet(U, param);
    (void)_Ret_U;
    double Hk_dMa = GetHk_dMa(U, param, d);
    double Hs_dMa = GetHs_dMa(U, param, d);  // laminar Hs
    double Ret_dMa = GetRet_dMa(U, param, d);
    double HsRet = Hs * Ret;
    double HsRet_dMa = Hs_dMa * Ret + Hs * Ret_dMa;
    double Hk_inv = 1.0 / Hk;
    double one_m = 1.0 - Hk_inv;
    double num = 2.0 * 1.1 * one_m * one_m * Hk_inv;
    double num_dHk = 2.0 * 1.1 *
                     (2.0 * one_m * (Hk_inv * Hk_inv) * Hk_inv +
                      one_m * one_m * (-Hk_inv * Hk_inv));
    return num_dHk * Hk_dMa / HsRet - num / (HsRet * HsRet) * HsRet_dMa;
}

double GetCDi_dMa(const Vector4d& U, const Param& param,
                  const ThermoMaDeriv& d) {
    if (param.turb) {
        double cDi_dMa = 0.0;
        double cDi_val = 0.0;
        if (!param.wake) {
            auto [v, _vU] = GetCDi_TurbWall(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dMa += GetCDi_TurbWall_dMa(U, param, d);
        }
        {
            auto [v, _vU] = GetCDi_Outer(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dMa += GetCDi_Outer_dMa(U, param, d);
        }
        {
            auto [v, _vU] = GetCDi_LamStress(U, param);
            (void)_vU;
            cDi_val += v;
            cDi_dMa += GetCDi_LamStress_dMa(U, param, d);
        }
        auto [cDil, _cDil_U] =
            param.wake ? GetCDi_LamWake(U, param) : GetCDi_Lam(U, param);
        (void)_cDil_U;
        double cDil_dMa = param.wake ? GetCDi_LamWake_dMa(U, param, d)
                                     : GetCDi_Lam_dMa(U, param, d);
        if (cDil > cDi_val) cDi_dMa = cDil_dMa;
        if (param.wake) cDi_dMa *= 2.0;
        return cDi_dMa;
    }
    return GetCDi_Lam_dMa(U, param, d);
}

double GetCfxt_dMa(const Vector4d& U, double x, const Param& param,
                   const ThermoMaDeriv& d) {
    double th = std::max(U(0), 1e-12);
    return GetCf_dMa(U, param, d) * x / th;
}

double GetCDixt_dMa(const Vector4d& U, double x, const Param& param,
                    const ThermoMaDeriv& d) {
    double th = std::max(U(0), 1e-12);
    return GetCDi_dMa(U, param, d) * x / th;
}

// ∂GetDamp/∂Ma at fixed U. Chain: GetDamp uses Hk(Ma) via Hmi, ar, da, af;
// and Ret(Ma) via lr = log(Ret)/log(10). Laminar amplification only.
double GetDamp_dMa(const Vector4d& U, const Param& param,
                   const ThermoMaDeriv& d) {
    if (param.Minf <= 0) return 0.0;
    auto [Hk, _Hk_U] = GetHk(U, param);
    (void)_Hk_U;
    double Hk_dM = GetHk_dMa(U, param, d);
    if (Hk < 1.05) {
        Hk = 1.05;
        Hk_dM = 0.0;
    }
    double Hmi = 1.0 / (Hk - 1.0);
    double Hmi_dM = -Hmi * Hmi * Hk_dM;
    double aa = 2.492 * std::pow(Hmi, 0.43);
    double aa_dM = 0.43 * aa / Hmi * Hmi_dM;
    double bb = std::tanh(14.0 * Hmi - 9.24);
    double bb_dM = (1.0 - bb * bb) * 14.0 * Hmi_dM;
    double lrc = aa + 0.7 * (bb + 1.0);
    double lrc_dM = aa_dM + 0.7 * bb_dM;
    auto [Ret_v, _r] = GetRet(U, param);
    (void)_r;
    double Ret_dM = GetRet_dMa(U, param, d);
    constexpr double lten = 2.302585092994046;
    double lr = std::log(Ret_v) / lten;
    double lr_dM = Ret_dM / (Ret_v * lten);
    constexpr double dl = 0.1;
    if (lr < lrc - dl) return 0.0;
    double rn = (lr - (lrc - dl)) / (2.0 * dl);
    double rn_dM = (lr_dM - lrc_dM) / (2.0 * dl);
    double rf, rf_dM;
    if (rn >= 1.0) {
        rf = 1.0;
        rf_dM = 0.0;
    } else {
        rf = 3.0 * rn * rn - 2.0 * rn * rn * rn;
        rf_dM = (6.0 * rn - 6.0 * rn * rn) * rn_dM;
    }
    double ar = 3.87 * Hmi - 2.52;
    double ar_dM = 3.87 * Hmi_dM;
    double ex = std::exp(-ar * ar);
    double ex_dM = ex * (-2.0 * ar * ar_dM);
    double da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
    double da_dM = 0.028 * Hk_dM - 0.0345 * ex_dM;
    double af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi + 3.0 * Hmi * Hmi * Hmi +
                0.1 * std::exp(-20.0 * Hmi);
    double af_dM =
        (2.7 - 11.0 * Hmi + 9.0 * Hmi * Hmi - 2.0 * std::exp(-20.0 * Hmi)) *
        Hmi_dM;
    double th_loc = U(0);
    return (rf_dM * af * da + rf * af_dM * da + rf * af * da_dM) / th_loc;
}

// ∂GetCteq/∂Ma — lifted from inline lambda in ComputeStationLinMa's turb
// branch. Chains through Hs(Ma), Us(Ma), Hk(Ma), Ret(Ma), Hkc(Ma).
double GetCteq_dMa(const Vector4d& U_in, const Param& p_in,
                   const ThermoMaDeriv& d) {
    if (!p_in.turb || p_in.Minf <= 0) return 0.0;
    double CC = 0.5 / (p_in.GA * p_in.GA * p_in.GB);
    double C = p_in.GC;
    auto [Hk, _h] = GetHk(U_in, p_in);
    (void)_h;
    auto [Hs, _hs] = GetHs(U_in, p_in);
    (void)_hs;
    auto [H, _H] = GetH(U_in);
    (void)_H;
    auto [Ret, _r] = GetRet(U_in, p_in);
    (void)_r;
    auto [Us, _us] = GetUs(U_in, p_in);
    (void)_us;
    double Hk_dM = GetHk_dMa(U_in, p_in, d);
    double Hs_dM = GetHs_dMa(U_in, p_in, d);
    double Us_dM = GetUs_dMa(U_in, p_in, d);
    double Ret_dM = GetRet_dMa(U_in, p_in, d);
    double Hkc, Hkc_dM;
    if (p_in.wake) {
        if (Hk < 1.00005) {
            Hk = 1.00005;
            Hk_dM = 0;
        }
        Hkc = Hk - 1.0;
        Hkc_dM = Hk_dM;
    } else {
        if (Hk < 1.05) {
            Hk = 1.05;
            Hk_dM = 0;
        }
        Hkc = Hk - 1.0 - C / Ret;
        Hkc_dM = Hk_dM + C / (Ret * Ret) * Ret_dM;
        if (Hkc < 0.01) {
            Hkc = 0.01;
            Hkc_dM = 0.0;
        }
    }
    double Hkc2 = Hkc * Hkc;
    double Hk2 = Hk * Hk;
    double num = CC * Hs * (Hk - 1.0) * Hkc2;
    double num_dM = CC * (Hs_dM * (Hk - 1.0) * Hkc2 + Hs * Hk_dM * Hkc2 +
                          Hs * (Hk - 1.0) * 2.0 * Hkc * Hkc_dM);
    double den_v = (1.0 - Us) * H * Hk2;
    double den_dM = (-Us_dM) * H * Hk2 + (1.0 - Us) * H * 2.0 * Hk * Hk_dM;
    double cteq_v = std::sqrt(num / den_v);
    return 0.5 / cteq_v * (num_dM / den_v - num / (den_v * den_v) * den_dM);
}

// ∂GetCttr/∂Ma. GetCttr = c·cteq where c = CtauC·exp(-CtauE/(Hk-1)).
// Both c and cteq depend on Ma (via Hk and the cteq chain). p must have
// wake=false (matching GetCttr's internal `param_nowake`).
double GetCttr_dMa(const Vector4d& U, const Param& p_nowake,
                   const ThermoMaDeriv& d) {
    if (p_nowake.Minf <= 0) return 0.0;
    auto [Hk, _h] = GetHk(U, p_nowake);
    (void)_h;
    double Hk_dM = GetHk_dMa(U, p_nowake, d);
    if (Hk < 1.05) {
        Hk = 1.05;
        Hk_dM = 0.0;
    }
    auto [cteq, _ct] = GetCteq(U, p_nowake);
    (void)_ct;
    double cteq_dM = GetCteq_dMa(U, p_nowake, d);
    double E = p_nowake.CtauE;
    double Cc = p_nowake.CtauC;
    double c = Cc * std::exp(-E / (Hk - 1.0));
    double c_dM = c * E / ((Hk - 1.0) * (Hk - 1.0)) * Hk_dM;
    return c_dM * cteq + c * cteq_dM;
}

// ============================================================================
// ComputeStationLinRet
//   Returns Vector3d lin_Ret = sum_k ∂R/∂Ret_k · Ret_k for the three residual
//   equations (Rmom, Rshape, Rlag) at the station.
//
//   Mirrors ResidualStation's formula structure but tracks ONLY the Ret-
//   derivative path (closures' explicit ∂/∂Ret · Ret contribution). U and
//   xi are held fixed (this is the partial for the adjoint).
// ============================================================================
Vector3d ComputeStationLinRet(const Param& param, const Vector2d& x,
                              const Matrix42d& Uin, const Vector2d& Aux) {
    // Copy U and remove wake gap from delta* (matching ResidualStation).
    Matrix42d U = Uin;
    U.row(1) -= Aux.transpose();
    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);

    // xlog (matches ResidualStation including simi override).
    double xlog;
    if (param.simi) {
        xlog = 1.0;
    } else {
        xlog = std::log(std::max(x(1), 1e-12) / std::max(x(0), 1e-12));
    }
    double dx;
    if (param.simi) {
        dx = 0.5 * (x(0) + x(1));
    } else {
        dx = x(1) - x(0);
    }

    // uelog (used in Rshape via "(2*Hss/Hs + 1 - H - Hw)*uelog" — Hs is the
    // only Ret-dep factor here).
    double uelog;
    if (param.simi) {
        uelog = 1.0;
    } else {
        auto [uk1, _uk1] = GetUk(U1(3), param);
        auto [uk2, _uk2] = GetUk(U2(3), param);
        (void)_uk1;
        (void)_uk2;
        uelog = std::log(std::max(uk2, 1e-12) / std::max(uk1, 1e-12));
    }

    // Upwind factor (Hk-only → no Ret dep, but we need the value for
    // weighting Ret-deps of cfxtu/cDixt).
    auto [upw, _upw_U] = GetUpw(U1, U2, param);
    (void)_upw_U;

    // Hs (turb-dep on Ret), needed for Rshape coefficients
    auto [Hs1, _Hs1_U] = GetHs(U1, param);
    auto [Hs2, _Hs2_U] = GetHs(U2, param);
    (void)_Hs1_U;
    (void)_Hs2_U;
    double Hs1_dRet = GetHs_dRet(U1, param);
    double Hs2_dRet = GetHs_dRet(U2, param);
    auto [Ret1_at_1, _R1_U] = GetRet(U1, param);
    auto [Ret2_at_2, _R2_U] = GetRet(U2, param);
    (void)_R1_U;
    (void)_R2_U;
    double Hs1_dRet_Ret = Hs1_dRet * Ret1_at_1;
    double Hs2_dRet_Ret = Hs2_dRet * Ret2_at_2;
    double Hs = 0.5 * (Hs1 + Hs2);
    double Hs_dRet_Ret = 0.5 * (Hs1_dRet_Ret + Hs2_dRet_Ret);

    // Hslog = log(Hs2/Hs1). ∂Hslog/∂Ret_1·Ret_1 = (-1/Hs1)*Hs1_dRet_Ret
    double Hslog_dRet_Ret = 0.0;
    if (!param.simi) {
        Hslog_dRet_Ret = -Hs1_dRet_Ret / std::max(Hs1, 1e-12) +
                         Hs2_dRet_Ret / std::max(Hs2, 1e-12);
    }

    // Hss is Mach-only, no Ret dep.
    auto [Hss1, _Hss1] = GetHss(U1, param);
    auto [Hss2, _Hss2] = GetHss(U2, param);
    (void)_Hss1;
    (void)_Hss2;
    double Hss = 0.5 * (Hss1 + Hss2);

    // cfxt symmetric-average derivative
    double cfxt1_dRet_at_U1 = GetCfxt_dRet(U1, x(0), param);
    double cfxt2_dRet_at_U2 = GetCfxt_dRet(U2, x(1), param);
    Vector4d Um = 0.5 * (U1 + U2);
    double xm = 0.5 * (x(0) + x(1));
    double cfxtm_dRet_at_Um = GetCfxt_dRet(Um, xm, param);
    auto [Retm_at_Um, _Retm_U] = GetRet(Um, param);
    (void)_Retm_U;
    double cfxt_dRet_Ret = 0.25 * cfxt1_dRet_at_U1 * Ret1_at_1 +
                           0.5 * cfxtm_dRet_at_Um * Retm_at_Um +
                           0.25 * cfxt2_dRet_at_U2 * Ret2_at_2;

    // cfxtu = upwind(cfxt1, cfxt2) — no middle.
    double cfxt1_dRet_Ret = cfxt1_dRet_at_U1 * Ret1_at_1;
    double cfxt2_dRet_Ret = cfxt2_dRet_at_U2 * Ret2_at_2;
    double cfxtu_dRet_Ret = (1.0 - upw) * cfxt1_dRet_Ret + upw * cfxt2_dRet_Ret;

    // cDixt = upwind(cDixt1, cDixt2)
    double cDixt1_dRet_at_U1 = GetCDixt_dRet(U1, x(0), param);
    double cDixt2_dRet_at_U2 = GetCDixt_dRet(U2, x(1), param);
    double cDixt_dRet_Ret = (1.0 - upw) * cDixt1_dRet_at_U1 * Ret1_at_1 +
                            upw * cDixt2_dRet_at_U2 * Ret2_at_2;

    // ---- Rmom ----
    // Rmom = thlog + (2+H+Hw-Ms)*uelog - 0.5*xlog*cfxt
    // Only cfxt has Ret dep.
    double Rmom_dRet_Ret = -0.5 * xlog * cfxt_dRet_Ret;

    // ---- Rshape ----
    // Rshape = Hslog + (2*Hss/Hs + 1 - H - Hw)*uelog + xlog*(0.5*cfxtu - cDixt)
    // ∂Rshape/∂Ret_k · Ret_k contributions:
    //   Hslog: Hslog_dRet_Ret
    //   2*Hss/Hs → ∂(2*Hss/Hs)/∂Ret · Ret = 2*Hss * (-1/Hs^2) * Hs_dRet_Ret
    //   xlog * (0.5*cfxtu - cDixt): xlog*(0.5*cfxtu_dRet_Ret - cDixt_dRet_Ret)
    double Rshape_dRet_Ret = Hslog_dRet_Ret +
                             (-2.0 * Hss / (Hs * Hs)) * uelog * Hs_dRet_Ret +
                             xlog * (0.5 * cfxtu_dRet_Ret - cDixt_dRet_Ret);

    // ---- Rlag ----
    double Rlag_dRet_Ret = 0.0;
    if (param.turb) {
        // Rlag = Clag·(cteq - ald·saa)·dx - 2·de·salog + 2·de·(uq·dx -
        // uelog)·Cuq Ret-dep components: Clag(Us), cteq (Hs,Hkc,Us), uq (cf
        // upwind + Ret_avg direct arg). de, saa, salog, uelog are Ret-indep.
        auto [Us1, _Us1_U] = GetUs(U1, param);
        auto [Us2, _Us2_U] = GetUs(U2, param);
        (void)_Us1_U;
        (void)_Us2_U;
        double Us1_dRet = GetUs_dRet(U1, param);
        double Us2_dRet = GetUs_dRet(U2, param);
        double Us_avg = 0.5 * (Us1 + Us2);
        double Us_dRet_Ret =
            0.5 * (Us1_dRet * Ret1_at_1 + Us2_dRet * Ret2_at_2);

        auto [de1, _de1_U] = GetDe(U1, param);
        auto [de2, _de2_U] = GetDe(U2, param);
        (void)_de1_U;
        (void)_de2_U;
        double de = 0.5 * (de1 + de2);

        auto [cteq1, _cteq1_U] = GetCteq(U1, param);
        auto [cteq2, _cteq2_U] = GetCteq(U2, param);
        (void)_cteq1_U;
        (void)_cteq2_U;
        double cteq1_dRet = GetCteq_dRet(U1, param);
        double cteq2_dRet = GetCteq_dRet(U2, param);
        double cteq_v = (1.0 - upw) * cteq1 + upw * cteq2;
        double cteq_dRet_Ret =
            (1.0 - upw) * cteq1_dRet * Ret1_at_1 + upw * cteq2_dRet * Ret2_at_2;

        // uq inputs: dsa, cf (upwind), Hk (upwind), Ret_avg = 0.5(Ret1+Ret2)
        auto [Hk1, _Hk1_U] = GetHk(U1, param);
        auto [Hk2, _Hk2_U] = GetHk(U2, param);
        (void)_Hk1_U;
        (void)_Hk2_U;
        double Hk_up = (1.0 - upw) * Hk1 + upw * Hk2;
        auto [cf1, _cf1_U] = GetCf(U1, param);
        auto [cf2, _cf2_U] = GetCf(U2, param);
        (void)_cf1_U;
        (void)_cf2_U;
        double cf_up = (1.0 - upw) * cf1 + upw * cf2;
        double cf1_dRet = GetCf_dRet(U1, param);
        double cf2_dRet = GetCf_dRet(U2, param);
        double cf_dRet_Ret =
            (1.0 - upw) * cf1_dRet * Ret1_at_1 + upw * cf2_dRet * Ret2_at_2;
        double Ret_avg = 0.5 * (Ret1_at_1 + Ret2_at_2);
        double Ret_avg_dRet_Ret = 0.5 * (Ret1_at_1 + Ret2_at_2);
        double dsa = 0.5 * (U1(1) + U2(1));
        double uq_dRet_arg = GetUq_dRet_arg(cf_up, Hk_up, Ret_avg, dsa, param);
        double uq_dcf = GetUq_dcf(dsa, param);
        double uq_dRet_Ret =
            uq_dcf * cf_dRet_Ret + uq_dRet_arg * Ret_avg_dRet_Ret;

        // saa value (no Ret dep)
        double sa1_v = U1(2), sa2_v = U2(2);
        double saa = (1.0 - upw) * sa1_v + upw * sa2_v;

        double Klag = param.SlagK;
        double beta = param.GB;
        double Clag = Klag / beta / (1.0 + Us_avg);
        double Clag_dRet_Ret = -Clag / (1.0 + Us_avg) * Us_dRet_Ret;
        double ald = param.wake ? param.Dlr : 1.0;
        double Cuq = param.Cuq;

        // T1 = Clag·(cteq - ald·saa)·dx; T3 = 2·de·(uq·dx - uelog)·Cuq
        double T1_dRet_Ret = Clag_dRet_Ret * (cteq_v - ald * saa) * dx +
                             Clag * cteq_dRet_Ret * dx;
        double T3_dRet_Ret = 2.0 * de * uq_dRet_Ret * dx * Cuq;
        Rlag_dRet_Ret = T1_dRet_Ret + T3_dRet_Ret;
    } else if (!param.simi) {
        // Laminar non-simi: Rlag = sa2 - sa1 - damp*dx
        double damp1_dRet = GetDamp_dRet(U1, param);
        double damp2_dRet = GetDamp_dRet(U2, param);
        double damp_dRet_Ret =
            0.5 * (damp1_dRet * Ret1_at_1 + damp2_dRet * Ret2_at_2);
        Rlag_dRet_Ret = -damp_dRet_Ret * dx;
    }
    // Laminar simi: Rlag = sa1 + sa2 — no Ret dep, stays at 0.

    return Vector3d(Rmom_dRet_Ret, Rshape_dRet_Ret, Rlag_dRet_Ret);
}

// ============================================================================
// ComputeStationLinMa
//   Returns Vector3d lin_Ma = ∂R/∂Ma · Ma for (Rmom, Rshape, Rlag) at a
//   station. Mirrors ComputeStationLinRet but tracks the Ma chain through
//   the 6 thermo params (KTb, KTl, H0, cps, mu0, rho0). For Ma <= 0 returns
//   zero (incompressible case).
// ============================================================================
Vector3d ComputeStationLinMa(const Param& param, const Vector2d& x,
                             const Matrix42d& Uin, const Vector2d& Aux) {
    if (param.Minf <= 0) return Vector3d::Zero();

    // Recover constants from param to compute thermo derivatives.
    double Vinf = param.Vinf;
    double Ma = param.Minf;
    // muinf = ρ_inf · Vinf · chord / Re, but param doesn't carry rhoinf or
    // chord directly. ComputeThermoMaDeriv only uses param.muinf, param.mu0,
    // param.rho0, param.KTb, param.H0, gam, Tsrat — all in Param.
    ThermoMaDeriv d = ComputeThermoMaDeriv(param, Vinf, 0.0);

    // Set up U pairs (mirror ResidualStation prelude).
    Matrix42d U = Uin;
    U.row(1) -= Aux.transpose();
    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);

    double xlog, dx;
    if (param.simi) {
        xlog = 1.0;
        dx = 0.5 * (x(0) + x(1));
    } else {
        xlog = std::log(std::max(x(1), 1e-12) / std::max(x(0), 1e-12));
        dx = x(1) - x(0);
    }

    // uelog and ∂uelog/∂Ma: uelog = log(uk2/uk1), uk depends on Ma via KTl.
    double uelog;
    double uelog_dMa = 0;
    if (param.simi) {
        uelog = 1.0;
    } else {
        auto [uk1, _u1] = GetUk(U1(3), param);
        (void)_u1;
        auto [uk2, _u2] = GetUk(U2(3), param);
        (void)_u2;
        uelog = std::log(std::max(uk2, 1e-12) / std::max(uk1, 1e-12));
        double uk1_dMa = GetUk_dMa(U1(3), param, d);
        double uk2_dMa = GetUk_dMa(U2(3), param, d);
        uelog_dMa =
            uk2_dMa / std::max(uk2, 1e-12) - uk1_dMa / std::max(uk1, 1e-12);
    }

    // Upwind factor — depends on Hk → Ma. Need upw_dMa.
    auto [upw, _upw_U] = GetUpw(U1, U2, param);
    (void)_upw_U;
    // upw formula: upw = 1 - 0.5·exp(-Hls·Huc). Hls depends on Hk1, Hk2
    // (Ma-dep). Huc depends on Hk2 (Ma-dep). Compute upw_dMa from GetUpw's
    // chain:
    //   Hk1, Hk2 both have Ma dep. Other constants don't.
    double upw_dMa;
    {
        auto [Hk1, _Hk1_U] = GetHk(U1, param);
        (void)_Hk1_U;
        auto [Hk2, _Hk2_U] = GetHk(U2, param);
        (void)_Hk2_U;
        double Hk1_dMa = GetHk_dMa(U1, param, d);
        double Hk2_dMa = GetHk_dMa(U2, param, d);
        double Hut = 1.0;
        double C = param.wake ? 1.0 : 5.0;
        double Huc = C * Hut / (Hk2 * Hk2);
        double Huc_dMa = -2.0 * C * Hut / (Hk2 * Hk2 * Hk2) * Hk2_dMa;
        double aa_ratio = (Hk2 - 1.0) / (Hk1 - 1.0);
        double sga = aa_ratio >= 0 ? 1.0 : -1.0;
        double la = std::log(sga * aa_ratio);
        double la_dMa = Hk2_dMa / (Hk2 - 1.0) - Hk1_dMa / (Hk1 - 1.0);
        double Hls = la * la;
        double Hls_dMa = 2.0 * la * la_dMa;
        if (Hls > 15.0) {
            Hls = 15.0;
            Hls_dMa = 0.0;
        }
        double ef = std::exp(-Hls * Huc);
        upw_dMa = -0.5 * ef * (-Hls_dMa * Huc - Hls * Huc_dMa);
    }

    // Hs averaged
    auto [Hs1, _Hs1_U] = GetHs(U1, param);
    (void)_Hs1_U;
    auto [Hs2, _Hs2_U] = GetHs(U2, param);
    (void)_Hs2_U;
    double Hs1_dMa = GetHs_dMa(U1, param, d);
    double Hs2_dMa = GetHs_dMa(U2, param, d);
    double Hs = 0.5 * (Hs1 + Hs2);
    double Hs_dMa = 0.5 * (Hs1_dMa + Hs2_dMa);

    // Hslog
    double Hslog_dMa = 0.0;
    if (!param.simi) {
        Hslog_dMa =
            -Hs1_dMa / std::max(Hs1, 1e-12) + Hs2_dMa / std::max(Hs2, 1e-12);
    }

    // H, Hw — no Ma dep (Hw depends on Aux, no Ma)
    auto [H1, _H1_U] = GetH(U1);
    (void)_H1_U;
    auto [H2, _H2_U] = GetH(U2);
    (void)_H2_U;
    double H = 0.5 * (H1 + H2);
    auto [Hw1, _Hw1_U] = GetHw(U1, Aux(0));
    (void)_Hw1_U;
    auto [Hw2, _Hw2_U] = GetHw(U2, Aux(1));
    (void)_Hw2_U;
    double Hw = 0.5 * (Hw1 + Hw2);

    // Ms = upwinded M2 (0.5 symmetric average per the formula in GetMach2 use)
    auto [Ms1, _Ms1_U] = GetMach2(U1, param);
    (void)_Ms1_U;
    auto [Ms2, _Ms2_U] = GetMach2(U2, param);
    (void)_Ms2_U;
    double Ms1_dMa = GetMach2_dMa(U1, param, d);
    double Ms2_dMa = GetMach2_dMa(U2, param, d);
    double Ms = 0.5 * (Ms1 + Ms2);
    double Ms_dMa = 0.5 * (Ms1_dMa + Ms2_dMa);

    // cfxt — symmetric average + middle
    double cfxt1_dMa = GetCfxt_dMa(U1, x(0), param, d);
    double cfxt2_dMa = GetCfxt_dMa(U2, x(1), param, d);
    Vector4d Um = 0.5 * (U1 + U2);
    double xm = 0.5 * (x(0) + x(1));
    double cfxtm_dMa = GetCfxt_dMa(Um, xm, param, d);
    double cfxt_dMa = 0.25 * cfxt1_dMa + 0.5 * cfxtm_dMa + 0.25 * cfxt2_dMa;

    // cfxtu, cDixt — upwinded (upw has Ma dep)
    auto [cfxt1_v, _cfxt1_U, _cfxt1_x] = GetCfxt(U1, x(0), param);
    auto [cfxt2_v, _cfxt2_U, _cfxt2_x] = GetCfxt(U2, x(1), param);
    (void)_cfxt1_U;
    (void)_cfxt1_x;
    (void)_cfxt2_U;
    (void)_cfxt2_x;
    double cfxtu_dMa = -upw_dMa * cfxt1_v + upw_dMa * cfxt2_v +
                       (1.0 - upw) * cfxt1_dMa + upw * cfxt2_dMa;
    auto [cDixt1_v, _cD1_U, _cD1_x] = GetCDixt(U1, x(0), param);
    auto [cDixt2_v, _cD2_U, _cD2_x] = GetCDixt(U2, x(1), param);
    (void)_cD1_U;
    (void)_cD1_x;
    (void)_cD2_U;
    (void)_cD2_x;
    double cDixt1_dMa = GetCDixt_dMa(U1, x(0), param, d);
    double cDixt2_dMa = GetCDixt_dMa(U2, x(1), param, d);
    double cDixt_dMa = -upw_dMa * cDixt1_v + upw_dMa * cDixt2_v +
                       (1.0 - upw) * cDixt1_dMa + upw * cDixt2_dMa;

    // Hss averaged
    auto [Hss1, _Hss1_U] = GetHss(U1, param);
    (void)_Hss1_U;
    auto [Hss2, _Hss2_U] = GetHss(U2, param);
    (void)_Hss2_U;
    double Hss = 0.5 * (Hss1 + Hss2);
    double Hss_dMa =
        0.5 * (GetHss_dMa(U1, param, d) + GetHss_dMa(U2, param, d));

    // ---- Rmom = thlog + (2+H+Hw-Ms)*uelog - 0.5·xlog·cfxt ----
    // thlog has no Ma dep. H, Hw have no Ma dep. Ms, uelog, cfxt do.
    // ∂Rmom/∂Ma = (-Ms_dMa)·uelog + (2+H+Hw-Ms)·uelog_dMa - 0.5·xlog·cfxt_dMa
    double Rmom_dMa = (-Ms_dMa) * uelog + (2.0 + H + Hw - Ms) * uelog_dMa -
                      0.5 * xlog * cfxt_dMa;

    // ---- Rshape = Hslog + (2·Hss/Hs + 1 - H - Hw)·uelog + xlog·(0.5·cfxtu -
    // cDixt) ∂Rshape/∂Ma = Hslog_dMa
    //               + d/dMa[2·Hss/Hs]·uelog + (2·Hss/Hs+1-H-Hw)·uelog_dMa
    //               + xlog·(0.5·cfxtu_dMa - cDixt_dMa)
    double Hss_over_Hs_dMa =
        2.0 * Hss_dMa / Hs - 2.0 * Hss / (Hs * Hs) * Hs_dMa;
    double Rshape_dMa = Hslog_dMa + Hss_over_Hs_dMa * uelog +
                        (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog_dMa +
                        xlog * (0.5 * cfxtu_dMa - cDixt_dMa);

    // ---- Rlag ----
    double Rlag_dMa = 0.0;
    if (param.turb) {
        // Rlag = Clag·(cteq - ald·saa)·dx - 2·de·salog + 2·de·(uq·dx -
        // uelog)·Cuq de uses Hk (Ma-dep). cteq, Us, cf all Ma-dep. uq via cf,
        // Hk, Ret.
        auto [Us1, _u1U] = GetUs(U1, param);
        (void)_u1U;
        auto [Us2, _u2U] = GetUs(U2, param);
        (void)_u2U;
        double Us1_dMa = GetUs_dMa(U1, param, d);
        double Us2_dMa = GetUs_dMa(U2, param, d);
        double Us_avg = 0.5 * (Us1 + Us2);
        double Us_dMa = 0.5 * (Us1_dMa + Us2_dMa);

        auto [de1, _de1U] = GetDe(U1, param);
        (void)_de1U;
        auto [de2, _de2U] = GetDe(U2, param);
        (void)_de2U;
        double de = 0.5 * (de1 + de2);
        // Correct ∂de/∂Ma chain. From GetDe (bl_properties.cpp:778-802):
        //   aa = 3.15 + 1.72/(Hk-1)  →  aa_dMa = -1.72/(Hk-1)² · Hk_dMa
        //   de = th·aa + ds          [capped at 12·th]
        // At fixed U: ∂de/∂Ma = th · aa_dMa  (if not capped).
        auto [Hk1, _Hk1U] = GetHk(U1, param);
        (void)_Hk1U;
        auto [Hk2, _Hk2U] = GetHk(U2, param);
        (void)_Hk2U;
        double Hk1_dMa = GetHk_dMa(U1, param, d);
        double Hk2_dMa = GetHk_dMa(U2, param, d);
        double th1 = std::max(U1(0), 1e-12);
        double th2 = std::max(U2(0), 1e-12);
        auto de_dMa_at = [&](double Hk_v, double Hk_dMa_v, double th_v,
                             double ds_v) {
            if (Hk_v < 1.05) return 0.0;
            double aa = 3.15 + 1.72 / (Hk_v - 1.0);
            // Cap condition: de = th·aa + ds compared against 12·th.
            if (th_v * aa + ds_v >= 12.0 * th_v) return 0.0;
            return th_v * (-1.72 / ((Hk_v - 1.0) * (Hk_v - 1.0))) * Hk_dMa_v;
        };
        double de1_dMa = de_dMa_at(Hk1, Hk1_dMa, th1, U1(1));
        double de2_dMa = de_dMa_at(Hk2, Hk2_dMa, th2, U2(1));
        double de_dMa = 0.5 * (de1_dMa + de2_dMa);

        // cteq upwinded with upw
        auto [cteq1, _ct1U] = GetCteq(U1, param);
        (void)_ct1U;
        auto [cteq2, _ct2U] = GetCteq(U2, param);
        (void)_ct2U;
        // GetCteq_dMa: need to derive. For brevity, chain through Hs, Us, Ret,
        // Hkc — all already available via *_dMa helpers. Inline:
        auto compute_cteq_dMa = [&](const Vector4d& Uin, const Param& p_in) {
            double CC = 0.5 / (p_in.GA * p_in.GA * p_in.GB);
            double C = p_in.GC;
            auto [Hk, _h] = GetHk(Uin, p_in);
            (void)_h;
            auto [Hs, _hs] = GetHs(Uin, p_in);
            (void)_hs;
            auto [H, _H] = GetH(Uin);
            (void)_H;
            auto [Ret, _r] = GetRet(Uin, p_in);
            (void)_r;
            auto [Us, _us] = GetUs(Uin, p_in);
            (void)_us;
            double Hk_dM = GetHk_dMa(Uin, p_in, d);
            double Hs_dM = GetHs_dMa(Uin, p_in, d);
            double Us_dM = GetUs_dMa(Uin, p_in, d);
            double Ret_dM = GetRet_dMa(Uin, p_in, d);
            double Hkc, Hkc_dM;
            if (p_in.wake) {
                if (Hk < 1.00005) {
                    Hk = 1.00005;
                    Hk_dM = 0;
                }
                Hkc = Hk - 1.0;
                Hkc_dM = Hk_dM;
            } else {
                if (Hk < 1.05) {
                    Hk = 1.05;
                    Hk_dM = 0;
                }
                Hkc = Hk - 1.0 - C / Ret;
                Hkc_dM = Hk_dM + C / (Ret * Ret) * Ret_dM;
                if (Hkc < 0.01) {
                    Hkc = 0.01;
                    Hkc_dM = 0.0;
                }
            }
            double Hkc2 = Hkc * Hkc;
            double Hk2 = Hk * Hk;
            double num = CC * Hs * (Hk - 1.0) * Hkc2;
            double num_dM =
                CC * (Hs_dM * (Hk - 1.0) * Hkc2 + Hs * Hk_dM * Hkc2 +
                      Hs * (Hk - 1.0) * 2.0 * Hkc * Hkc_dM);
            double den_v = (1.0 - Us) * H * Hk2;
            double den_dM =
                (-Us_dM) * H * Hk2 + (1.0 - Us) * H * 2.0 * Hk * Hk_dM;
            double cteq_v = std::sqrt(num / den_v);
            return 0.5 / cteq_v *
                   (num_dM / den_v - num / (den_v * den_v) * den_dM);
        };
        double cteq1_dMa = compute_cteq_dMa(U1, param);
        double cteq2_dMa = compute_cteq_dMa(U2, param);
        double cteq_v = (1.0 - upw) * cteq1 + upw * cteq2;
        double cteq_dMa = -upw_dMa * cteq1 + upw_dMa * cteq2 +
                          (1.0 - upw) * cteq1_dMa + upw * cteq2_dMa;

        // cf upwinded, Hk upwinded
        auto [cf1, _cf1U] = GetCf(U1, param);
        (void)_cf1U;
        auto [cf2, _cf2U] = GetCf(U2, param);
        (void)_cf2U;
        double cf1_dMa = GetCf_dMa(U1, param, d);
        double cf2_dMa = GetCf_dMa(U2, param, d);
        double cf_up = (1.0 - upw) * cf1 + upw * cf2;
        double cf_up_dMa = -upw_dMa * cf1 + upw_dMa * cf2 +
                           (1.0 - upw) * cf1_dMa + upw * cf2_dMa;
        double Hk_up = (1.0 - upw) * Hk1 + upw * Hk2;
        double Hk_up_dMa = -upw_dMa * Hk1 + upw_dMa * Hk2 +
                           (1.0 - upw) * Hk1_dMa + upw * Hk2_dMa;
        auto [Ret1_v, _R1] = GetRet(U1, param);
        (void)_R1;
        auto [Ret2_v, _R2] = GetRet(U2, param);
        (void)_R2;
        double Ret_avg = 0.5 * (Ret1_v + Ret2_v);
        double Ret_avg_dMa =
            0.5 * (GetRet_dMa(U1, param, d) + GetRet_dMa(U2, param, d));
        double dsa = 0.5 * (U1(1) + U2(1));
        // uq chain: ∂uq/∂cf · cf_dMa + ∂uq/∂Hk · Hk_dMa + ∂uq/∂Ret · Ret_dMa
        // From GetUq formula: uq = (0.5·cf - (Hkc/AHk)²) / (β·ds),
        // Hkc=Hk-1-C/Ret, AHk = A·Hk. ∂uq/∂Hk has multiple terms. For brevity,
        // use GetUq_dRet_arg path for Ret, and finite-difference approach is
        // too expensive — fallback: build uq from cf, Hk, Ret, and chain
        // analytically.
        double beta = param.GB;
        double A = param.GA;
        double Cc = param.GC;
        if (param.wake) {
            A *= param.Dlr;
            Cc = 0.0;
        }
        double Hk_lim2 = param.wake ? 1.00005 : 1.05;
        double Hk_for_uq = Hk_up;
        double Hk_for_uq_dMa = Hk_up_dMa;
        if (Hk_for_uq < Hk_lim2) {
            Hk_for_uq = Hk_lim2;
            Hk_for_uq_dMa = 0.0;
        }
        double Hkc_uq = Hk_for_uq - 1.0 - Cc / std::max(Ret_avg, 1e-10);
        double Hkc_uq_dMa =
            Hk_for_uq_dMa + Cc / (Ret_avg * Ret_avg) * Ret_avg_dMa;
        if (Hkc_uq < 0.01) {
            Hkc_uq = 0.01;
            Hkc_uq_dMa = 0.0;
        }
        double AHk = A * Hk_for_uq;
        double AHk_dMa = A * Hk_for_uq_dMa;
        double Hkc_over_AHk = Hkc_uq / AHk;
        double Hkc_over_AHk_dMa =
            Hkc_uq_dMa / AHk - Hkc_uq / (AHk * AHk) * AHk_dMa;
        double ut_dMa = 0.5 * cf_up_dMa - 2.0 * Hkc_over_AHk * Hkc_over_AHk_dMa;
        double ds_safe = std::max(dsa, 1e-12);
        double uq_dMa = ut_dMa / (beta * ds_safe);

        // saa, salog have no Ma dep.
        double saa = (1.0 - upw) * U1(2) + upw * U2(2);
        double saa_dMa =
            -upw_dMa * U1(2) + upw_dMa * U2(2);  // upwind has Ma-dep

        double sa1_safe = std::max(U1(2), 1e-10);
        double sa2_safe = std::max(U2(2), 1e-10);
        double salog = std::log(sa2_safe / sa1_safe);
        (void)salog;

        double Klag = param.SlagK;
        double Clag = Klag / beta / (1.0 + Us_avg);
        double Clag_dMa = -Clag / (1.0 + Us_avg) * Us_dMa;
        double ald = param.wake ? param.Dlr : 1.0;
        double Cuq = param.Cuq;

        // T1 = Clag·(cteq - ald·saa)·dx
        // T2 = -2·de·salog
        // T3 = 2·de·(uq·dx - uelog)·Cuq
        double T1_dMa = Clag_dMa * (cteq_v - ald * saa) * dx +
                        Clag * (cteq_dMa - ald * saa_dMa) * dx;
        double T2_dMa = -2.0 * de_dMa * salog;
        // need salog value...we have it from `std::log(sa2_safe/sa1_safe)`
        double salog_v = std::log(sa2_safe / sa1_safe);
        T2_dMa = -2.0 * de_dMa * salog_v;
        // uq value (needed for T3)
        double ut = 0.5 * cf_up - (Hkc_uq / AHk) * (Hkc_uq / AHk);
        double uq_v = ut / (beta * ds_safe);
        double T3_dMa = 2.0 * de_dMa * (uq_v * dx - uelog) * Cuq +
                        2.0 * de * (uq_dMa * dx - uelog_dMa) * Cuq;
        Rlag_dMa = T1_dMa + T2_dMa + T3_dMa;
    } else if (!param.simi) {
        // Laminar non-simi: Rlag = sa2 - sa1 - damp·dx
        // GetDamp's Ma-dep is through Ret (via lr = log(Ret)/log(10)) and Hk
        // (via Hmi). Both have _dMa. Use full chain.
        // For brevity, compute damp_dMa via GetDamp using ε-style helpers:
        // This is small effort; the chain mirrors GetDamp_dRet but with
        // Use the lifted GetDamp_dMa helper.
        double damp1_dMa = GetDamp_dMa(U1, param, d);
        double damp2_dMa = GetDamp_dMa(U2, param, d);
        double damp_dMa = 0.5 * (damp1_dMa + damp2_dMa);
        Rlag_dMa = -damp_dMa * dx;
    }
    // Laminar simi: Rlag = sa1 + sa2, no Ma dep.

    // Pre-multiply by Ma (matches lin_Ret · Ret convention).
    return Vector3d(Rmom_dMa, Rshape_dMa, Rlag_dMa) * Ma;
}

}  // namespace detail
using namespace detail;

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
    // sa: NO outer clamp. mfoil uses raw U[2]. Clamping sa at 1e-12 prevents
    // Newton from converging the laminar-stagnation residual (Rlag = sa1+sa2)
    // below 2e-12, pinning |R| at the plateau. Turbulent salog has its own
    // inner clamp at 1e-10 below.
    double sa1 = U1(2);
    double sa2 = U2(2);

    // Compressibility-corrected velocities
    auto [uk1, uk1_U] = GetUk(U1(3), param);
    auto [uk2, uk2_U] = GetUk(U2(3), param);

    // Log changes (with protection against small/negative values)
    double th1_safe = std::max(th1, 1e-12);
    double th2_safe = std::max(th2, 1e-12);
    double thlog = std::log(th2_safe / th1_safe);
    Vector8d thlog_U = Vector8d::Zero();
    thlog_U(0) = (th1 > 1e-12) ? -1.0 / th1 : 0.0;
    thlog_U(4) = (th2 > 1e-12) ? 1.0 / th2 : 0.0;

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
    cfxt_U << 0.25 * cfxt1_U1 + 0.25 * cfxtm_Um,
        0.25 * cfxtm_Um + 0.25 * cfxt2_U2;
    Vector2d cfxt_x(0.25 * cfxt1_x1 + 0.25 * cfxtm_xm,
                    0.25 * cfxtm_xm + 0.25 * cfxt2_x2);

    // Momentum equation
    double Rmom = thlog + (2.0 + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt;
    Vector8d Rmom_U = thlog_U + (H_U + Hw_U - Ms_U) * uelog +
                      (2.0 + H + Hw - Ms) * uelog_U - 0.5 * xlog * cfxt_U;
    Vector2d Rmom_x = -0.5 * xlog_x * cfxt - 0.5 * xlog * cfxt_x;

    // cDi * x / theta (dissipation)
    auto [cDixt1, cDixt1_U1, cDixt1_x1] = GetCDixt(U1, x(0), param);
    auto [cDixt2, cDixt2_U2, cDixt2_x2] = GetCDixt(U2, x(1), param);
    auto [cDixt, cDixt_U] =
        Upwind(upw, upw_U, cDixt1, cDixt1_U1, cDixt2, cDixt2_U2);
    Vector2d cDixt_x((1.0 - upw) * cDixt1_x1, upw * cDixt2_x2);

    // Upwinded cf*x/theta
    auto [cfxtu, cfxtu_U] =
        Upwind(upw, upw_U, cfxt1, cfxt1_U1, cfxt2, cfxt2_U2);
    Vector2d cfxtu_x((1.0 - upw) * cfxt1_x1, upw * cfxt2_x2);

    // Density shape parameter Hss
    auto [Hss1, Hss1_U1] = GetHss(U1, param);
    auto [Hss2, Hss2_U2] = GetHss(U2, param);
    auto [Hss, Hss_U] =
        Upwind(0.5, Vector8d::Zero(), Hss1, Hss1_U1, Hss2, Hss2_U2);

    // Shape parameter equation
    double Rshape = Hslog + (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog +
                    xlog * (0.5 * cfxtu - cDixt);
    Vector8d Rshape_U =
        Hslog_U +
        (2.0 * Hss_U / Hs - 2.0 * Hss / (Hs * Hs) * Hs_U - H_U - Hw_U) * uelog +
        (2.0 * Hss / Hs + 1.0 - H - Hw) * uelog_U +
        xlog * (0.5 * cfxtu_U - cDixt_U);
    Vector2d Rshape_x =
        xlog_x * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x - cDixt_x);

    // Amplification or shear lag equation
    double Rlag;
    Vector8d Rlag_U;
    Vector2d Rlag_x;

    if (param.turb) {
        // Turbulent: shear lag equation
        // Protect against small/negative ctau values
        double sa1_safe = std::max(sa1, 1e-10);
        double sa2_safe = std::max(sa2, 1e-10);
        double salog = std::log(sa2_safe / sa1_safe);
        Vector8d salog_U = Vector8d::Zero();
        salog_U(2) = (sa1 > 1e-10) ? -1.0 / sa1 : 0.0;
        salog_U(6) = (sa2 > 1e-10) ? 1.0 / sa2 : 0.0;

        auto [de1, de1_U1] = GetDe(U1, param);
        auto [de2, de2_U2] = GetDe(U2, param);
        auto [de, de_U] =
            Upwind(0.5, Vector8d::Zero(), de1, de1_U1, de2, de2_U2);

        auto [Us1, Us1_U1] = GetUs(U1, param);
        auto [Us2, Us2_U2] = GetUs(U2, param);
        auto [Us, Us_U] =
            Upwind(0.5, Vector8d::Zero(), Us1, Us1_U1, Us2, Us2_U2);

        // Hk, upwinded (needed for uq)
        auto [Hk1, Hk1_U1] = GetHk(U1, param);
        auto [Hk2, Hk2_U2] = GetHk(U2, param);
        auto [Hk, Hk_U] = Upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2);

        // Ret, averaged
        auto [Ret1, Ret1_U1] = GetRet(U1, param);
        auto [Ret2, Ret2_U2] = GetRet(U2, param);
        auto [Ret, Ret_U] =
            Upwind(0.5, Vector8d::Zero(), Ret1, Ret1_U1, Ret2, Ret2_U2);

        // cf, upwinded
        auto [cf1, cf1_U1] = GetCf(U1, param);
        auto [cf2, cf2_U2] = GetCf(U2, param);
        auto [cf, cf_U] = Upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2);

        // ds averaged
        double dsa = 0.5 * (U1(1) + U2(1));
        Vector8d dsa_U = Vector8d::Zero();
        dsa_U(1) = 0.5;
        dsa_U(5) = 0.5;

        // uq = equilibrium velocity gradient
        auto [uq, uq_U] =
            GetUq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param);

        auto [cteq1, cteq1_U1] = GetCteq(U1, param);
        auto [cteq2, cteq2_U2] = GetCteq(U2, param);
        auto [cteq, cteq_U] =
            Upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2);

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
        double Cuq = param.Cuq;

        // Full shear lag equation with uq correction
        Rlag = Clag * (cteq - ald * saa) * dx - 2.0 * de * salog +
               2.0 * de * (uq * dx - uelog) * Cuq;
        Rlag_U = Clag_U * (cteq - ald * saa) * dx +
                 Clag * (cteq_U - ald * saa_U) * dx - 2.0 * de_U * salog -
                 2.0 * de * salog_U + 2.0 * de_U * (uq * dx - uelog) * Cuq +
                 2.0 * de * (uq_U * dx - uelog_U) * Cuq;
        Rlag_x = Clag * (cteq - ald * saa) * dx_x + 2.0 * de * uq * dx_x;
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
            auto [damp, damp_U] =
                Upwind(0.5, Vector8d::Zero(), damp1, damp1_U1, damp2, damp2_U2);

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

    // ∂R/∂Re·Re via the analytic Ret chain. Re only enters R through
    // Ret = ρ·ue·θ/μ_0 in BL closures, with ∂Ret/∂Re = Ret/Re, so
    // ∂R/∂Re = lin_Ret/Re. ComputeStationLinRet replicates the residual
    // structure but tracks only the Ret-derivative path.
    result.lin_Ret = ComputeStationLinRet(param, x, Uin, Aux);

    // ∂R/∂Ma·Ma via the analytic Ma chain through KTb, KTl, H0, mu0, rho0.
    result.lin_Ma = ComputeStationLinMa(param, x, Uin, Aux);

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

    // Get TE gap (recompute from foil.x, matching mfoil's wake_sys)
    TEInfo te = GetTEInfo(state.foil.x);
    double hTE = te.hTE;

    // Get shear stress at wake start
    // Use transition ctau if not already turbulent
    param.turb = true;
    param.wake = false;

    double ctl, ctu;
    Vector4d ctl_Ul = Vector4d::Zero();
    Vector4d ctu_Uu = Vector4d::Zero();
    double ctl_dMa = 0.0;
    double ctu_dMa = 0.0;
    bool any_lam_TE = !state.vsol.turb[il] || !state.vsol.turb[iu];
    // GetCttr internally forces wake=false; param has turb=true already.
    Param p_cttr = param;
    p_cttr.wake = false;
    ThermoMaDeriv d_Ma;
    if (any_lam_TE) {
        d_Ma = ComputeThermoMaDeriv(p_cttr, p_cttr.Vinf, 0.0);
    }

    if (state.vsol.turb[il]) {
        ctl = Ul(2);
        ctl_Ul(2) = 1.0;
    } else {
        auto [ct, ct_U] = GetCttr(Ul, param);
        ctl = ct;
        ctl_Ul = ct_U;
        ctl_dMa = GetCttr_dMa(Ul, p_cttr, d_Ma);
    }

    if (state.vsol.turb[iu]) {
        ctu = Uu(2);
        ctu_Uu(2) = 1.0;
    } else {
        auto [ct, ct_U] = GetCttr(Uu, param);
        ctu = ct;
        ctu_Uu = ct_U;
        ctu_dMa = GetCttr_dMa(Uu, p_cttr, d_Ma);
    }

    // Theta-weighted average of ctau
    double thsum = Ul(0) + Uu(0);
    double ctw = (ctl * Ul(0) + ctu * Uu(0)) / thsum;

    Vector4d ctw_Ul =
        (ctl_Ul * Ul(0) + (ctl - ctw) * Vector4d(1, 0, 0, 0)) / thsum;
    Vector4d ctw_Uu =
        (ctu_Uu * Uu(0) + (ctu - ctw) * Vector4d(1, 0, 0, 0)) / thsum;

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

    // ∂R/∂Ma · Ma. Only R(2) = Uw(2) - ctw has Ma dependence (through ctl/ctu
    // when the corresponding TE side is laminar). R(0), R(1) are theta/ds
    // balances with no closure Ma dep.
    if (any_lam_TE) {
        double ctw_dMa = (ctl_dMa * Ul(0) + ctu_dMa * Uu(0)) / thsum;
        double Ma_v = param.Minf;
        result.lin_Ma_times_Ma(2) = -ctw_dMa * Ma_v;
    }

    return result;
}

ResidualLin ResidualTransition(FoilState& state, const Param& param_base,
                               const Vector2d& x, const Matrix42d& U,
                               const Vector2d& Aux) {
    // Calculates the combined lam + turb residual for a transition station
    // U1 should be laminar; U2 should be turbulent
    // Calculates and linearizes the transition location in the process

    ResidualLin result;
    result.value = Vector3d::Zero();
    result.lin_U = Matrix38d::Zero();
    result.lin_x = Matrix32d::Zero();

    // States
    Vector4d U1 = U.col(0);
    Vector4d U2 = U.col(1);
    Vector4d Z4 = Vector4d::Zero();

    // Interval
    double x1 = x(0);
    double x2 = x(1);
    double dx = x2 - x1;

    // Determine transition location (xt) using amplification equation
    double xt = x1 + 0.5 * dx;  // Initial guess
    double ncrit = param_base.ncrit;
    constexpr int nNewton = 30;

    double sa1 = U1(2);  // Amplification at point 1

    // Newton iteration to find xt
    for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
        double w2 = (xt - x1) / dx;
        double w1 = 1.0 - w2;

        // State at xt (linear interpolation)
        Vector4d Ut = w1 * U1 + w2 * U2;
        Vector4d Ut_xt = (U2 - U1) / dx;  // d(Ut)/d(xt)
        Ut(2) = ncrit;                    // Amplification at transition = ncrit
        Ut_xt(2) = 0.0;                   // sa at xt is always ncrit

        // Get damping at U1 and Ut
        Param param = param_base;
        param.turb = false;
        auto [damp1, damp1_U1] = GetDamp(U1, param);
        auto [dampt, dampt_Ut] = GetDamp(Ut, param);
        dampt_Ut(2) = 0.0;  // No dependence on sa at transition

        // Residual for transition location: ncrit - sa1 - 0.5*(xt-x1)*(damp1 +
        // dampt)
        double Rxt = ncrit - sa1 - 0.5 * (xt - x1) * (damp1 + dampt);
        double Rxt_xt =
            -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);

        double dxt = -Rxt / Rxt_xt;

        // Limit step size
        double dmax = 0.2 * dx * (1.1 - static_cast<double>(iNewton) / nNewton);
        if (std::abs(dxt) > dmax) {
            dxt = dxt * dmax / std::abs(dxt);
        }

        if (std::abs(Rxt) < 1e-14) break;
        if (iNewton < nNewton - 1) xt += dxt;
    }

    // Store transition location
    state.vsol.xt = xt;

    // Prepare for xt linearizations
    double w2 = (xt - x1) / dx;
    double w1 = 1.0 - w2;
    Vector4d Ut = w1 * U1 + w2 * U2;
    Vector4d Ut_xt = (U2 - U1) / dx;
    Ut(2) = ncrit;
    Ut_xt(2) = 0.0;

    Param param = param_base;
    param.turb = false;
    auto [damp1, damp1_U1] = GetDamp(U1, param);
    auto [dampt, dampt_Ut] = GetDamp(Ut, param);
    dampt_Ut(2) = 0.0;

    double Rxt_xt =
        -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);

    // Rxt_U = d(Rxt)/d(U1, U2) at fixed xt
    Vector8d Rxt_U = Vector8d::Zero();
    Rxt_U.head<4>() = -0.5 * (xt - x1) * (damp1_U1 + dampt_Ut * w1);
    Rxt_U.tail<4>() = -0.5 * (xt - x1) * dampt_Ut * w2;
    Rxt_U(2) -= 1.0;  // d(Rxt)/d(sa1) = -1

    // Ut_x at fixed xt
    Vector4d Ut_x1 = (U2 - U1) * (w2 - 1.0) / dx;
    Vector4d Ut_x2 = (U2 - U1) * (-w2) / dx;
    Ut_x1(2) = 0.0;
    Ut_x2(2) = 0.0;

    double Rxt_x1 =
        0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_x1);
    double Rxt_x2 = -0.5 * (xt - x1) * dampt_Ut.dot(Ut_x2);

    // Sensitivity of xt w.r.t. U,x from Rxt(xt,U,x) = 0 constraint
    Vector8d xt_U = -Rxt_U / Rxt_xt;
    Vector4d xt_U1 = xt_U.head<4>();
    Vector4d xt_U2 = xt_U.tail<4>();
    double xt_x1 = -Rxt_x1 / Rxt_xt;
    double xt_x2 = -Rxt_x2 / Rxt_xt;

    // Include derivatives w.r.t. xt in Ut_x
    Ut_x1 += Ut_xt * xt_x1;
    Ut_x2 += Ut_xt * xt_x2;

    // Sensitivity of Ut w.r.t. U1 and U2
    Matrix4d Ut_U1 =
        w1 * Matrix4d::Identity() + (U2 - U1) * xt_U1.transpose() / dx;
    Matrix4d Ut_U2 =
        w2 * Matrix4d::Identity() + (U2 - U1) * xt_U2.transpose() / dx;

    // Laminar and turbulent states at transition
    Vector4d Utl = Ut;
    Matrix4d Utl_U1 = Ut_U1;
    Matrix4d Utl_U2 = Ut_U2;
    Vector4d Utl_x1 = Ut_x1;
    Vector4d Utl_x2 = Ut_x2;
    Utl(2) = ncrit;
    Utl_U1.row(2).setZero();
    Utl_U2.row(2).setZero();
    Utl_x1(2) = 0.0;
    Utl_x2(2) = 0.0;

    Vector4d Utt = Ut;
    Matrix4d Utt_U1 = Ut_U1;
    Matrix4d Utt_U2 = Ut_U2;
    Vector4d Utt_x1 = Ut_x1;
    Vector4d Utt_x2 = Ut_x2;

    // Set turbulent shear coefficient in Utt
    param.turb = true;
    auto [cttr, cttr_Ut] = GetCttr(Ut, param);
    Utt(2) = cttr;
    Utt_U1.row(2) = (cttr_Ut.transpose() * Ut_U1).transpose();
    Utt_U2.row(2) = (cttr_Ut.transpose() * Ut_U2).transpose();
    Utt_x1(2) = cttr_Ut.dot(Ut_x1);
    Utt_x2(2) = cttr_Ut.dot(Ut_x2);

    // Laminar residual: x1 to xt
    param.turb = false;
    Matrix42d U_lam;
    U_lam.col(0) = U1;
    U_lam.col(1) = Utl;
    Vector2d x_lam(x1, xt);
    auto Rl = ResidualStation(param, x_lam, U_lam, Aux);

    // Turbulent residual: xt to x2
    param.turb = true;
    Matrix42d U_turb;
    U_turb.col(0) = Utt;
    U_turb.col(1) = U2;
    Vector2d x_turb(xt, x2);
    auto Rt = ResidualStation(param, x_turb, U_turb, Aux);

    // Combined residual
    result.value = Rl.value + Rt.value;

    // Combined Jacobian using chain rule
    // R_U1 = Rl_U1 + Rl_Utl * Utl_U1 + Rl_x[:,1] * xt_U1 + Rt_Utt * Utt_U1 +
    // Rt_x[:,0] * xt_U1
    Matrix34d Rl_U1 = Rl.lin_U.block<3, 4>(0, 0);
    Matrix34d Rl_Utl = Rl.lin_U.block<3, 4>(0, 4);
    Matrix34d Rt_Utt = Rt.lin_U.block<3, 4>(0, 0);
    Matrix34d Rt_U2 = Rt.lin_U.block<3, 4>(0, 4);

    Matrix34d R_U1 = Rl_U1 + Rl_Utl * Utl_U1 +
                     Rl.lin_x.col(1) * xt_U1.transpose() + Rt_Utt * Utt_U1 +
                     Rt.lin_x.col(0) * xt_U1.transpose();

    Matrix34d R_U2 = Rl_Utl * Utl_U2 + Rl.lin_x.col(1) * xt_U2.transpose() +
                     Rt_Utt * Utt_U2 + Rt_U2 +
                     Rt.lin_x.col(0) * xt_U2.transpose();

    result.lin_U.block<3, 4>(0, 0) = R_U1;
    result.lin_U.block<3, 4>(0, 4) = R_U2;

    // Combined x linearization
    result.lin_x.col(0) = Rl.lin_x.col(0) + Rl.lin_x.col(1) * xt_x1 +
                          Rt.lin_x.col(0) * xt_x1 + Rl_Utl * Utl_x1 +
                          Rt_Utt * Utt_x1;
    result.lin_x.col(1) = Rt.lin_x.col(1) + Rl.lin_x.col(1) * xt_x2 +
                          Rt.lin_x.col(0) * xt_x2 + Rl_Utl * Utl_x2 +
                          Rt_Utt * Utt_x2;

    // ∂R/∂Ma·Ma for transition station — mirrors the xt(Re) chain below
    // exactly, just with GetDamp_dMa / GetCttr_dMa / GetRet/Ret value
    // replaced by their Ma-specific analogues. At fixed U, xt depends on
    // Ma via damp1, dampt (both have Ma dep through GetRet, GetHk in
    // GetDamp); Utt(2) = cttr also depends on Ma via GetCttr (c·cteq).
    //
    // Steps mirror the Re chain (a)–(d):
    //   (a) Rl.lin_Ma + Rt.lin_Ma (direct closure Ma-deps; xt held fixed)
    //   (b) ∂xt/∂Ma · (Rl.lin_x.col(1) + Rt.lin_x.col(0))
    //   (c) ∂xt/∂Ma · (Rl_Utl · Utl_xt + Rt_Utt · Utt_xt)
    //   (d) (Rt_Utt · ê_2) · cttr_dMa · Ma
    {
        Param p_lam = param_base;
        p_lam.turb = false;
        ThermoMaDeriv d_Ma =
            ComputeThermoMaDeriv(param_base, param_base.Vinf, 0.0);
        double damp1_dMa = GetDamp_dMa(U1, p_lam, d_Ma);
        double dampt_dMa = GetDamp_dMa(Ut, p_lam, d_Ma);
        double Ma_v = param_base.Minf;
        double Rxt_dMa_Ma = -0.5 * (xt - x1) * (damp1_dMa + dampt_dMa) * Ma_v;
        double xt_dMa_Ma = -Rxt_dMa_Ma / Rxt_xt;

        Vector3d dR_dxt_direct = Rl.lin_x.col(1) + Rt.lin_x.col(0);
        Vector4d Utl_xt_ma = Ut_xt;
        Utl_xt_ma(2) = 0;
        Vector4d Utt_xt_ma = Ut_xt;
        Utt_xt_ma(2) = cttr_Ut.dot(Ut_xt);
        Vector3d dR_dxt_via_Ut = Rl_Utl * Utl_xt_ma + Rt_Utt * Utt_xt_ma;

        // cttr direct Ma chain (GetCttr uses non-wake turbulent param).
        Param p_nowake_turb = param_base;
        p_nowake_turb.wake = false;
        p_nowake_turb.turb = true;
        double cttr_dMa_Ma = GetCttr_dMa(Ut, p_nowake_turb, d_Ma) * Ma_v;
        Vector3d dR_via_cttr = Rt_Utt.col(2) * cttr_dMa_Ma;

        result.lin_Ma = Rl.lin_Ma + Rt.lin_Ma +
                        xt_dMa_Ma * (dR_dxt_direct + dR_dxt_via_Ut) +
                        dR_via_cttr;
    }

    // ∂R/∂Re·Re: at fixed U, xt is NOT independent of Re because Rxt(xt, U, Re)
    // = 0 uses GetDamp which has Ret-dep. By implicit FN: ∂xt/∂Re =
    // -Rxt_Re/Rxt_xt.
    //
    //   Rxt = ncrit - sa1 - 0.5·(xt-x1)·(damp1 + dampt)
    //   ∂Rxt/∂Re·Re at fixed U, xt = -0.5·(xt-x1)·(damp1_dRet·Ret1 +
    //   dampt_dRet·Ret_t)
    //
    // The combined transition residual depends on Re via:
    //   (a) Rl.lin_Ret + Rt.lin_Ret (direct closure Ret-deps with x_lam, U_lam
    //   fixed) (b) ∂xt/∂Re · (Rl.lin_x.col(1) + Rt.lin_x.col(0))      [direct x
    //   dep] (c) ∂xt/∂Re · (Rl_Utl · Utl_xt + Rt_Utt · Utt_xt)      [Utl, Utt
    //   depend on xt] (d) (Rt_Utt · ê_2) · cttr_dRet·Ret_t [cttr(Re) directly]
    {
        // damp_dRet at the two endpoints used in Rxt
        Param p_lam = param_base;
        p_lam.turb = false;
        double damp1_dRet = GetDamp_dRet(U1, p_lam);
        double dampt_dRet = GetDamp_dRet(Ut, p_lam);
        auto [Ret1_v, _Ret1_U] = GetRet(U1, p_lam);
        (void)_Ret1_U;
        auto [Rett_v, _Rett_U] = GetRet(Ut, p_lam);
        (void)_Rett_U;
        double Rxt_Re_Re =
            -0.5 * (xt - x1) * (damp1_dRet * Ret1_v + dampt_dRet * Rett_v);
        double xt_dRe_Re = -Rxt_Re_Re / Rxt_xt;

        // ∂R/∂xt direct (Rl uses x[1]=xt, Rt uses x[0]=xt):
        Vector3d dR_dxt_direct = Rl.lin_x.col(1) + Rt.lin_x.col(0);

        // ∂R/∂xt via Utl/Utt = Ut(xt):
        Vector4d Utl_xt = Ut_xt;
        Utl_xt(2) = 0.0;  // ncrit constant
        Vector4d Utt_xt = Ut_xt;
        Utt_xt(2) = cttr_Ut.dot(Ut_xt);  // cttr(Ut(xt)) chain through Ut_xt

        Vector3d dR_dxt_via_Ut = Rl_Utl * Utl_xt + Rt_Utt * Utt_xt;

        // ∂R/∂Re via cttr direct Ret-dep (cttr uses cteq + GetHk, no other Ret
        // path) cttr_dRet from chain: GetCttr returns C·exp(-E/(Hk-1))·cteq
        // with Hk and cteq from non-wake param. ∂cttr/∂Ret at fixed U =
        // C·exp(...)·cteq_dRet.
        Param p_nowake = param_base;
        p_nowake.wake = false;
        p_nowake.turb = true;
        auto [Hk_at_Ut, _Hk_U] = GetHk(Ut, p_nowake);
        (void)_Hk_U;
        double Hk_lim = Hk_at_Ut;
        if (Hk_lim < 1.05) Hk_lim = 1.05;
        double c =
            param_base.CtauC * std::exp(-param_base.CtauE / (Hk_lim - 1.0));
        double cteq_dRet_at_Ut = GetCteq_dRet(Ut, p_nowake);
        double cttr_dRet_Ret_t = c * cteq_dRet_at_Ut * Rett_v;
        // ∂R/∂Utt(2) is Rt_Utt.col(2) — apply cttr Ret-chain to it.
        Vector3d dR_via_cttr = Rt_Utt.col(2) * cttr_dRet_Ret_t;

        result.lin_Ret = Rl.lin_Ret + Rt.lin_Ret +
                         xt_dRe_Re * (dR_dxt_direct + dR_dxt_via_Ut) +
                         dR_via_cttr;
    }

    return result;
}

}  // namespace cppfoil
