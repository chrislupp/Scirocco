// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/boundary_layer.h"

#include <cmath>
#include <iostream>
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
    // Initialize boundary layer by marching from stagnation point
    // Matches mfoil's init_boundary_layer approach

    constexpr double Hmaxl = 3.8;  // Max H before laminar separation
    constexpr double Hmaxt = 2.5;  // Max H before turbulent separation
    constexpr int nNewton = 30;
    constexpr int iNswitch = 12;   // Switch to inverse after this many iterations

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
            std::cout << "  BL init si=" << si << " hitstag=" << hitstag
                      << " xi(0)=" << xi(0) << " xi(Nsurf-1)=" << xi(Nsurf-1)
                      << " K=" << K << std::endl;

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
                Matrix3d A = res.lin_U.block<3, 3>(0, 4) + res.lin_U.block<3, 3>(0, 0);
                Vector3d b = -res.value;
                Vector3d dU3 = A.partialPivLu().solve(b);

                // Under-relaxation
                double dm = std::max(std::abs(dU3(0) / Ust(0)), std::abs(dU3(1) / Ust(1)));
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

        if (si == 1) {
            std::cout << "  BL march si=" << si << " starting at i=" << i << " (i0=" << i0 << ")" << std::endl;
        }

        while (i < Nsurf) {
            // Initialize with previous state and new ue
            U.col(i) = U.col(i - 1);
            U(3, i) = ue(i);

            // Set ctau at transition and switch to turbulent mode
            if (tran) {
                auto [ct, ct_U] = GetCttr(U.col(i), param);
                std::cout << "  cttr from GetCttr: " << ct << " (Hk will determine this)" << std::endl;
                U(2, i) = ct;
                param.turb = true;  // Set turb BEFORE Newton so we use turbulent residuals
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
                    res = ResidualTransition(state, param, xpair, Upair, Auxpair);
                    if (iNewton == 0) {
                        std::cout << "  Tran Newton start: ctau=" << U(2,i) << " R=("
                                  << res.value(0) << "," << res.value(1) << "," << res.value(2) << ")"
                                  << " |R|=" << res.value.norm() << std::endl;
                    }
                } else {
                    res = ResidualStation(param, xpair, Upair, Auxpair);
                }

                // Debug for first few stations
                if (si == 1 && i <= 3 && iNewton == 0) {
                    std::cout << "  Init si=" << si << " i=" << i << " |R|=" << res.value.norm()
                              << " R=(" << res.value(0) << "," << res.value(1) << "," << res.value(2) << ")"
                              << " th=" << U(0,i) << " ds=" << U(1,i) << std::endl;
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
                    if (tran && iNewton < 5) {
                        std::cout << "  Tran Newton[" << iNewton << "] dU=(" << dU(0) << "," << dU(1) << "," << dU(2) << ")"
                                  << " A(2,2)=" << A(2,2) << " ctau_before=" << U(2,i) << std::endl;
                    }
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
                double dm = std::max(std::abs(dU(0) / std::max(U(0, i - 1), 1e-12)),
                                     std::abs(dU(1) / std::max(U(1, i - 1), 1e-12)));
                if (!direct) dm = std::max(dm, std::abs(dU(3) / std::max(U(3, i - 1), 1e-10)));
                if (param.turb) dm = std::max(dm, std::abs(dU(2) / std::max(U(2, i - 1), 1e-10)));
                else if (direct) dm = std::max(dm, std::abs(dU(2) / 10.0));

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
                    double Hkr = (xi(i) - xi(i - 1)) / std::max(U(0, i - 1), 1e-12);

                    if (param.wake) {
                        double H2 = Hkprev;
                        for (int k = 0; k < 6; ++k) {
                            H2 -= (H2 + 0.03 * Hkr * std::pow(H2 - 1, 3) - Hkprev) /
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
                if (tran && iNewton < 5) {
                    std::cout << "  Tran Newton[" << iNewton << "] ctau_after=" << U(2,i)
                              << " omega=" << omega << " direct=" << direct << std::endl;
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
                } else if (si < 2 && xi(i) > 1e-10) {
                    // Check if we're right after stagnation (tiny xi_prev)
                    if (xi(i - 1) < 1e-8 * xi(Nsurf - 1)) {
                        // Use Thwaites initialization at this station
                        double K_local = ue(i) / std::max(xi(i), 1e-12);
                        auto [th_new, ds_new] = ThwaitesInit(K_local, nu);
                        U(0, i) = th_new;
                        U(1, i) = ds_new;
                    } else {
                        // Extrapolate theta and delta* along surface
                        double xratio = std::sqrt(xi(i) / xi(i - 1));
                        U(0, i) = U(0, i - 1) * xratio;
                        U(1, i) = U(1, i - 1) * xratio;
                    }
                }
            }

            // Check for transition (laminar only)
            if (!param.turb && !tran && U(2, i) > param.ncrit) {
                tran = true;
                std::cout << "  Transition detected: si=" << si << " i=" << i
                          << " Is[i]=" << Is[i] << " n=" << U(2,i) << std::endl;
                continue;  // Redo station with transition
            }

            if (tran) {
                // Store transition location (param.turb already set before Newton)
                state.vsol.Xt(si, 0) = 0.5 * (xi(i - 1) + xi(i));
                tran = false;
                std::cout << "  Transition applied: si=" << si << " i=" << i
                          << " ctau=" << U(2,i) << std::endl;
            }

            ++i;
        }

        // Store states in global array
        for (int j = 0; j < Nsurf; ++j) {
            state.glob.U.col(Is[j]) = U.col(j);
        }
    }
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
        // Protect against small/negative ctau values
        double sa1_safe = std::max(sa1, 1e-10);
        double sa2_safe = std::max(sa2, 1e-10);
        double salog = std::log(sa2_safe / sa1_safe);
        Vector8d salog_U = Vector8d::Zero();
        salog_U(2) = (sa1 > 1e-10) ? -1.0 / sa1 : 0.0;
        salog_U(6) = (sa2 > 1e-10) ? 1.0 / sa2 : 0.0;

        auto [de1, de1_U1] = GetDe(U1, param);
        auto [de2, de2_U2] = GetDe(U2, param);
        auto [de, de_U] = Upwind(0.5, Vector8d::Zero(), de1, de1_U1, de2, de2_U2);

        auto [Us1, Us1_U1] = GetUs(U1, param);
        auto [Us2, Us2_U2] = GetUs(U2, param);
        auto [Us, Us_U] = Upwind(0.5, Vector8d::Zero(), Us1, Us1_U1, Us2, Us2_U2);

        // Hk, upwinded (needed for uq)
        auto [Hk1, Hk1_U1] = GetHk(U1, param);
        auto [Hk2, Hk2_U2] = GetHk(U2, param);
        auto [Hk, Hk_U] = Upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2);

        // Ret, averaged
        auto [Ret1, Ret1_U1] = GetRet(U1, param);
        auto [Ret2, Ret2_U2] = GetRet(U2, param);
        auto [Ret, Ret_U] = Upwind(0.5, Vector8d::Zero(), Ret1, Ret1_U1, Ret2, Ret2_U2);

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
        auto [uq, uq_U] = GetUq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param);

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
        double Cuq = param.Cuq;

        // Full shear lag equation with uq correction
        Rlag = Clag * (cteq - ald * saa) * dx - 2.0 * de * salog +
               2.0 * de * (uq * dx - uelog) * Cuq;
        Rlag_U = Clag_U * (cteq - ald * saa) * dx +
                 Clag * (cteq_U - ald * saa_U) * dx -
                 2.0 * de_U * salog - 2.0 * de * salog_U +
                 2.0 * de_U * (uq * dx - uelog) * Cuq +
                 2.0 * de * (uq_U * dx - uelog_U) * Cuq;
        Rlag_x = Clag * (cteq - ald * saa) * dx_x +
                 2.0 * de * uq * dx_x * Cuq;
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
    constexpr int nNewton = 20;

    double sa1 = U1(2);  // Amplification at point 1

    // Newton iteration to find xt
    for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
        double w2 = (xt - x1) / dx;
        double w1 = 1.0 - w2;

        // State at xt (linear interpolation)
        Vector4d Ut = w1 * U1 + w2 * U2;
        Vector4d Ut_xt = (U2 - U1) / dx;  // d(Ut)/d(xt)
        Ut(2) = ncrit;      // Amplification at transition = ncrit
        Ut_xt(2) = 0.0;     // sa at xt is always ncrit

        // Get damping at U1 and Ut
        Param param = param_base;
        param.turb = false;
        auto [damp1, damp1_U1] = GetDamp(U1, param);
        auto [dampt, dampt_Ut] = GetDamp(Ut, param);
        dampt_Ut(2) = 0.0;  // No dependence on sa at transition

        // Residual for transition location: ncrit - sa1 - 0.5*(xt-x1)*(damp1 + dampt)
        double Rxt = ncrit - sa1 - 0.5 * (xt - x1) * (damp1 + dampt);
        double Rxt_xt = -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);

        double dxt = -Rxt / Rxt_xt;

        // Limit step size
        double dmax = 0.2 * dx * (1.1 - static_cast<double>(iNewton) / nNewton);
        if (std::abs(dxt) > dmax) {
            dxt = dxt * dmax / std::abs(dxt);
        }

        if (std::abs(Rxt) < 1e-10) break;
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

    double Rxt_xt = -0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_xt);

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

    double Rxt_x1 = 0.5 * (damp1 + dampt) - 0.5 * (xt - x1) * dampt_Ut.dot(Ut_x1);
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
    Matrix4d Ut_U1 = w1 * Matrix4d::Identity() + (U2 - U1) * xt_U1.transpose() / dx;
    Matrix4d Ut_U2 = w2 * Matrix4d::Identity() + (U2 - U1) * xt_U2.transpose() / dx;

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
    // R_U1 = Rl_U1 + Rl_Utl * Utl_U1 + Rl_x[:,1] * xt_U1 + Rt_Utt * Utt_U1 + Rt_x[:,0] * xt_U1
    Matrix34d Rl_U1 = Rl.lin_U.block<3, 4>(0, 0);
    Matrix34d Rl_Utl = Rl.lin_U.block<3, 4>(0, 4);
    Matrix34d Rt_Utt = Rt.lin_U.block<3, 4>(0, 0);
    Matrix34d Rt_U2 = Rt.lin_U.block<3, 4>(0, 4);

    Matrix34d R_U1 = Rl_U1 + Rl_Utl * Utl_U1 + Rl.lin_x.col(1) * xt_U1.transpose() +
                     Rt_Utt * Utt_U1 + Rt.lin_x.col(0) * xt_U1.transpose();

    Matrix34d R_U2 = Rl_Utl * Utl_U2 + Rl.lin_x.col(1) * xt_U2.transpose() +
                     Rt_Utt * Utt_U2 + Rt_U2 + Rt.lin_x.col(0) * xt_U2.transpose();

    result.lin_U.block<3, 4>(0, 0) = R_U1;
    result.lin_U.block<3, 4>(0, 4) = R_U2;

    // Combined x linearization
    result.lin_x.col(0) = Rl.lin_x.col(0) + Rl.lin_x.col(1) * xt_x1 + Rt.lin_x.col(0) * xt_x1 +
                          Rl_Utl * Utl_x1 + Rt_Utt * Utt_x1;
    result.lin_x.col(1) = Rt.lin_x.col(1) + Rl.lin_x.col(1) * xt_x2 + Rt.lin_x.col(0) * xt_x2 +
                          Rl_Utl * Utl_x2 + Rt_Utt * Utt_x2;

    return result;
}

}  // namespace cppfoil
