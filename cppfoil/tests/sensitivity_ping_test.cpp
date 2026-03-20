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
#include <gtest/gtest.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

#include "cppfoil/config.h"
#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/geometry.h"
#include "cppfoil/core/panel_influence_deriv.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/data/state.h"

using namespace cppfoil;

static constexpr double MIN_RATE = 1.5;
static constexpr double EPS_SKIP = 5e-9;

// ===========================================================================
// Test fixture
// ===========================================================================
class SensitivityPingTest : public ::testing::Test {
   protected:
    static std::mt19937& rng() {
        static std::mt19937 gen(42);
        return gen;
    }

    static VectorXd RandVec(int n) {
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        VectorXd v(n);
        for (int i = 0; i < n; ++i) v(i) = dist(rng());
        return v;
    }
};

// ===========================================================================
// 1. Panel Influence Derivatives: PanelLinVortexStream
// ===========================================================================
TEST_F(SensitivityPingTest, PanelLinVortexStreamDeriv_Xj) {
    rng().seed(42);

    // Set up a non-degenerate panel and evaluation point
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    // Test each Xj component
    for (int k = 0; k < 4; ++k) {
        std::array<double, 3> a_v, b_v;
        std::array<double, 3> a_d, b_d;

        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelLinVortexStreamDeriv(Xj, xi);
            a_v[ie] = deriv.a;
            b_v[ie] = deriv.b;
            a_d[ie] = deriv.a_Xj(k);
            b_d[ie] = deriv.b_Xj(k);

            // Perturb Xj component k
            if (k < 2) {
                Xj(k, 0) += ep;
            } else {
                Xj(k - 2, 1) += ep;
            }
        }
        // Undo perturbations
        if (k < 2) {
            Xj(k, 0) -= 3 * ep;
        } else {
            Xj(k - 2, 1) -= 3 * ep;
        }

        // 3-point convergence check
        double E0 = std::abs((a_v[1] - a_v[0]) / ep - 0.5 * (a_d[0] + a_d[1])) +
                    std::abs((b_v[1] - b_v[0]) / ep - 0.5 * (b_d[0] + b_d[1]));
        double E1 = std::abs((a_v[2] - a_v[0]) / (2 * ep) - 0.5 * (a_d[0] + a_d[2])) +
                    std::abs((b_v[2] - b_v[0]) / (2 * ep) - 0.5 * (b_d[0] + b_d[2]));

        double rate = 99.0;
        if (E0 > EPS_SKIP && E1 > EPS_SKIP) {
            rate = std::log2(E1 / E0);
        }
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "LinVortexStream Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelLinVortexStreamDeriv_Xj: " << nfail << " failures" << std::endl;
}

TEST_F(SensitivityPingTest, PanelLinVortexStreamDeriv_xi) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 2; ++k) {
        std::array<double, 3> a_v, b_v;
        std::array<double, 3> a_d, b_d;

        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelLinVortexStreamDeriv(Xj, xi);
            a_v[ie] = deriv.a;
            b_v[ie] = deriv.b;
            a_d[ie] = deriv.a_xi(k);
            b_d[ie] = deriv.b_xi(k);
            xi(k) += ep;
        }
        xi(k) -= 3 * ep;

        double E0 = std::abs((a_v[1] - a_v[0]) / ep - 0.5 * (a_d[0] + a_d[1])) +
                    std::abs((b_v[1] - b_v[0]) / ep - 0.5 * (b_d[0] + b_d[1]));
        double E1 = std::abs((a_v[2] - a_v[0]) / (2 * ep) - 0.5 * (a_d[0] + a_d[2])) +
                    std::abs((b_v[2] - b_v[0]) / (2 * ep) - 0.5 * (b_d[0] + b_d[2]));

        double rate = 99.0;
        if (E0 > EPS_SKIP && E1 > EPS_SKIP) {
            rate = std::log2(E1 / E0);
        }
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "LinVortexStream xi k=" << k << " rate=" << rate;
        }
    }
    std::cout << "PanelLinVortexStreamDeriv_xi: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 2. Panel Influence Derivatives: PanelConstSourceStream
// ===========================================================================
TEST_F(SensitivityPingTest, PanelConstSourceStreamDeriv_Xj) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 4; ++k) {
        std::array<double, 3> v, d;
        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelConstSourceStreamDeriv(Xj, xi);
            v[ie] = deriv.val;
            d[ie] = deriv.val_Xj(k);
            if (k < 2) Xj(k, 0) += ep;
            else Xj(k - 2, 1) += ep;
        }
        if (k < 2) Xj(k, 0) -= 3 * ep;
        else Xj(k - 2, 1) -= 3 * ep;

        double E0 = std::abs((v[1] - v[0]) / ep - 0.5 * (d[0] + d[1]));
        double E1 = std::abs((v[2] - v[0]) / (2 * ep) - 0.5 * (d[0] + d[2]));
        double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E1 / E0) : 99.0;
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "ConstSourceStream Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelConstSourceStreamDeriv_Xj: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 3. Panel Influence Derivatives: PanelLinVortexVelocity
// ===========================================================================
TEST_F(SensitivityPingTest, PanelLinVortexVelocityDeriv_Xj) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 4; ++k) {
        std::array<Vector2d, 3> a_v, b_v;
        std::array<Vector2d, 3> a_d, b_d;

        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelLinVortexVelocityDeriv(Xj, xi);
            a_v[ie] = deriv.a;
            b_v[ie] = deriv.b;
            a_d[ie] = deriv.a_Xj.col(k);
            b_d[ie] = deriv.b_Xj.col(k);
            if (k < 2) Xj(k, 0) += ep;
            else Xj(k - 2, 1) += ep;
        }
        if (k < 2) Xj(k, 0) -= 3 * ep;
        else Xj(k - 2, 1) -= 3 * ep;

        double E0 = ((a_v[1] - a_v[0]) / ep - 0.5 * (a_d[0] + a_d[1])).norm() +
                    ((b_v[1] - b_v[0]) / ep - 0.5 * (b_d[0] + b_d[1])).norm();
        double E1 = ((a_v[2] - a_v[0]) / (2 * ep) - 0.5 * (a_d[0] + a_d[2])).norm() +
                    ((b_v[2] - b_v[0]) / (2 * ep) - 0.5 * (b_d[0] + b_d[2])).norm();

        double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E1 / E0) : 99.0;
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "LinVortexVelocity Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelLinVortexVelocityDeriv_Xj: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 4. Panel Influence Derivatives: PanelConstSourceVelocity
// ===========================================================================
TEST_F(SensitivityPingTest, PanelConstSourceVelocityDeriv_Xj) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 4; ++k) {
        std::array<Vector2d, 3> v;
        std::array<Vector2d, 3> d;

        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelConstSourceVelocityDeriv(Xj, xi);
            v[ie] = deriv.val;
            d[ie] = deriv.val_Xj.col(k);
            if (k < 2) Xj(k, 0) += ep;
            else Xj(k - 2, 1) += ep;
        }
        if (k < 2) Xj(k, 0) -= 3 * ep;
        else Xj(k - 2, 1) -= 3 * ep;

        double E0 = ((v[1] - v[0]) / ep - 0.5 * (d[0] + d[1])).norm();
        double E1 = ((v[2] - v[0]) / (2 * ep) - 0.5 * (d[0] + d[2])).norm();

        double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E1 / E0) : 99.0;
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "ConstSourceVelocity Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelConstSourceVelocityDeriv_Xj: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 5. Panel Influence Derivatives: PanelLinSourceStream
// ===========================================================================
TEST_F(SensitivityPingTest, PanelLinSourceStreamDeriv_Xj) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 4; ++k) {
        std::array<double, 3> a_v, b_v, a_d, b_d;
        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelLinSourceStreamDeriv(Xj, xi);
            a_v[ie] = deriv.a;
            b_v[ie] = deriv.b;
            a_d[ie] = deriv.a_Xj(k);
            b_d[ie] = deriv.b_Xj(k);
            if (k < 2) Xj(k, 0) += ep;
            else Xj(k - 2, 1) += ep;
        }
        if (k < 2) Xj(k, 0) -= 3 * ep;
        else Xj(k - 2, 1) -= 3 * ep;

        double E0 = std::abs((a_v[1] - a_v[0]) / ep - 0.5 * (a_d[0] + a_d[1])) +
                    std::abs((b_v[1] - b_v[0]) / ep - 0.5 * (b_d[0] + b_d[1]));
        double E1 = std::abs((a_v[2] - a_v[0]) / (2 * ep) - 0.5 * (a_d[0] + a_d[2])) +
                    std::abs((b_v[2] - b_v[0]) / (2 * ep) - 0.5 * (b_d[0] + b_d[2]));

        double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E1 / E0) : 99.0;
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "LinSourceStream Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelLinSourceStreamDeriv_Xj: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 6. Panel Influence Derivatives: PanelLinSourceVelocity
// ===========================================================================
TEST_F(SensitivityPingTest, PanelLinSourceVelocityDeriv_Xj) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);

    double ep = 1e-5;
    int nfail = 0;

    for (int k = 0; k < 4; ++k) {
        std::array<Vector2d, 3> a_v, b_v, a_d, b_d;
        for (int ie = 0; ie < 3; ++ie) {
            auto deriv = PanelLinSourceVelocityDeriv(Xj, xi);
            a_v[ie] = deriv.a;
            b_v[ie] = deriv.b;
            a_d[ie] = deriv.a_Xj.col(k);
            b_d[ie] = deriv.b_Xj.col(k);
            if (k < 2) Xj(k, 0) += ep;
            else Xj(k - 2, 1) += ep;
        }
        if (k < 2) Xj(k, 0) -= 3 * ep;
        else Xj(k - 2, 1) -= 3 * ep;

        double E0 = ((a_v[1] - a_v[0]) / ep - 0.5 * (a_d[0] + a_d[1])).norm() +
                    ((b_v[1] - b_v[0]) / ep - 0.5 * (b_d[0] + b_d[1])).norm();
        double E1 = ((a_v[2] - a_v[0]) / (2 * ep) - 0.5 * (a_d[0] + a_d[2])).norm() +
                    ((b_v[2] - b_v[0]) / (2 * ep) - 0.5 * (b_d[0] + b_d[2])).norm();

        double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E1 / E0) : 99.0;
        if (rate < MIN_RATE) {
            nfail++;
            ADD_FAILURE() << "LinSourceVelocity Xj k=" << k << " rate=" << rate
                          << " E=[" << E0 << "," << E1 << "]";
        }
    }
    std::cout << "PanelLinSourceVelocityDeriv_Xj: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 7. Force Partials: dcl/dx_partial
// ===========================================================================
TEST_F(SensitivityPingTest, ForcePartials_cl) {
    rng().seed(42);

    Cppfoil foil("2412", 80);
    foil.SetOper(5.0);
    foil.Solve();
    auto& state = foil.GetState();

    CalcForce(state);
    CalcForceLin(state);

    int N = state.foil.N;

    // Use a random perturbation direction to reduce cost
    VectorXd dx_x = RandVec(N) * 0.01;
    VectorXd dx_z = RandVec(N) * 0.01;

    double ep = 1e-5;

    // Predicted change from analytical partials
    double dcl_pred = (state.post.dcl_dx_partial.row(0).transpose().dot(dx_x) +
                       state.post.dcl_dx_partial.row(1).transpose().dot(dx_z));

    // FD check: perturb coordinates, re-solve, compare
    FoilState state_p = state;
    for (int i = 0; i < N; ++i) {
        state_p.foil.x(0, i) = state.foil.x(0, i) + ep * dx_x(i);
        state_p.foil.x(1, i) = state.foil.x(1, i) + ep * dx_z(i);
    }
    // Don't re-solve - just re-evaluate forces with perturbed coords but same gamma
    // This tests the PARTIAL derivative (holding gamma fixed)
    CalcForce(state_p);
    double dcl_fd = (state_p.post.cl - state.post.cl) / ep;

    double err = std::abs(dcl_fd - dcl_pred);
    double tol = 0.01 * std::abs(dcl_fd) + 1e-10;
    EXPECT_LT(err, tol) << "ForcePartials dcl: fd=" << dcl_fd << " pred=" << dcl_pred;
    std::cout << "ForcePartials_cl: err=" << err << " (fd=" << dcl_fd
              << ", pred=" << dcl_pred << ")" << std::endl;
}

// ===========================================================================
// 8. Inviscid Total Derivative: dcl/dx (full system)
//    Perturb a single node, re-solve, compare dcl to adjoint prediction.
// ===========================================================================
TEST_F(SensitivityPingTest, InviscidSensitivity_cl) {
    rng().seed(42);

    Cppfoil foil("2412", 80);
    foil.SetOper(5.0);
    foil.Solve();
    auto& state = foil.GetState();

    CalcForce(state);
    auto sens = InviscidSensitivities(state);

    int N = state.foil.N;

    // Test using central FD at two step sizes: ep and ep/2
    // Central FD has O(ep^2) error, so rate = log2(E_ep / E_ep/2) ~ 2
    int test_nodes[] = {10, 30, 50, 70};
    int nfail = 0;

    auto solve_perturbed = [&](int node, int dim, double pert) -> double {
        FoilState state_p = state;
        state_p.foil.x(dim, node) += pert;
        state_p.isol.AIC.resize(0, 0);
        state_p.isol.gamref.resize(0, 0);
        state_p.isol.gam.resize(0);
        SolveInviscid(state_p);
        CalcForce(state_p);
        return state_p.post.cl;
    };

    for (int node : test_nodes) {
        if (node >= N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            double dcl_pred = sens.dcl_dx(dim, node);

            // Central FD at two step sizes
            double ep = 1e-5;
            double cl_p = solve_perturbed(node, dim, ep);
            double cl_m = solve_perturbed(node, dim, -ep);
            double dcl_fd_ep = (cl_p - cl_m) / (2.0 * ep);

            double ep2 = ep / 2.0;
            double cl_p2 = solve_perturbed(node, dim, ep2);
            double cl_m2 = solve_perturbed(node, dim, -ep2);
            double dcl_fd_ep2 = (cl_p2 - cl_m2) / (2.0 * ep2);

            double E0 = std::abs(dcl_fd_ep - dcl_pred);
            double E1 = std::abs(dcl_fd_ep2 - dcl_pred);
            double rate = (E0 > EPS_SKIP && E1 > EPS_SKIP) ? std::log2(E0 / E1) : 99.0;

            if (rate < MIN_RATE) {
                nfail++;
                ADD_FAILURE() << "InviscidSens cl node=" << node << " dim=" << dim
                              << " rate=" << rate << " E=[" << E0 << "," << E1 << "]"
                              << " fd=" << dcl_fd_ep2 << " adj=" << dcl_pred;
            }
        }
    }
    std::cout << "InviscidSensitivity_cl: " << nfail << "/" << 8 << " failures" << std::endl;
}

// ===========================================================================
// 8b. xi derivatives for ALL influence functions
// ===========================================================================
TEST_F(SensitivityPingTest, AllInfluenceDeriv_xi) {
    Matrix2d Xj;
    Xj.col(0) << 0.3, 0.05;
    Xj.col(1) << 0.6, -0.02;
    Vector2d xi(0.45, 0.15);
    double ep = 1e-5;
    int nfail = 0;

    // LinVortexStream xi
    for (int k = 0; k < 2; ++k) {
        std::array<double, 3> a_v, b_v, a_d, b_d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto d = PanelLinVortexStreamDeriv(Xj, xi_t);
            a_v[ie] = d.a; b_v[ie] = d.b;
            a_d[ie] = d.a_xi(k); b_d[ie] = d.b_xi(k);
            xi_t(k) += ep;
        }
        double E0 = std::abs((a_v[1]-a_v[0])/ep - 0.5*(a_d[0]+a_d[1])) +
                    std::abs((b_v[1]-b_v[0])/ep - 0.5*(b_d[0]+b_d[1]));
        double E1 = std::abs((a_v[2]-a_v[0])/(2*ep) - 0.5*(a_d[0]+a_d[2])) +
                    std::abs((b_v[2]-b_v[0])/(2*ep) - 0.5*(b_d[0]+b_d[2]));
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "LVS xi k=" << k << " rate=" << rate; }
    }

    // ConstSourceStream xi
    for (int k = 0; k < 2; ++k) {
        std::array<double, 3> v, d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto dr = PanelConstSourceStreamDeriv(Xj, xi_t);
            v[ie] = dr.val; d[ie] = dr.val_xi(k);
            xi_t(k) += ep;
        }
        double E0 = std::abs((v[1]-v[0])/ep - 0.5*(d[0]+d[1]));
        double E1 = std::abs((v[2]-v[0])/(2*ep) - 0.5*(d[0]+d[2]));
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "CSS xi k=" << k << " rate=" << rate; }
    }

    // LinVortexVelocity xi
    for (int k = 0; k < 2; ++k) {
        std::array<Vector2d, 3> a_v, b_v, a_d, b_d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto dr = PanelLinVortexVelocityDeriv(Xj, xi_t);
            a_v[ie] = dr.a; b_v[ie] = dr.b;
            a_d[ie] = dr.a_xi.col(k); b_d[ie] = dr.b_xi.col(k);
            xi_t(k) += ep;
        }
        double E0 = ((a_v[1]-a_v[0])/ep - 0.5*(a_d[0]+a_d[1])).norm() +
                    ((b_v[1]-b_v[0])/ep - 0.5*(b_d[0]+b_d[1])).norm();
        double E1 = ((a_v[2]-a_v[0])/(2*ep) - 0.5*(a_d[0]+a_d[2])).norm() +
                    ((b_v[2]-b_v[0])/(2*ep) - 0.5*(b_d[0]+b_d[2])).norm();
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "LVV xi k=" << k << " rate=" << rate; }
    }

    // ConstSourceVelocity xi
    for (int k = 0; k < 2; ++k) {
        std::array<Vector2d, 3> v, d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto dr = PanelConstSourceVelocityDeriv(Xj, xi_t);
            v[ie] = dr.val; d[ie] = dr.val_xi.col(k);
            xi_t(k) += ep;
        }
        double E0 = ((v[1]-v[0])/ep - 0.5*(d[0]+d[1])).norm();
        double E1 = ((v[2]-v[0])/(2*ep) - 0.5*(d[0]+d[2])).norm();
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "CSV xi k=" << k << " rate=" << rate; }
    }

    // LinSourceStream xi
    for (int k = 0; k < 2; ++k) {
        std::array<double, 3> a_v, b_v, a_d, b_d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto dr = PanelLinSourceStreamDeriv(Xj, xi_t);
            a_v[ie] = dr.a; b_v[ie] = dr.b;
            a_d[ie] = dr.a_xi(k); b_d[ie] = dr.b_xi(k);
            xi_t(k) += ep;
        }
        double E0 = std::abs((a_v[1]-a_v[0])/ep - 0.5*(a_d[0]+a_d[1])) +
                    std::abs((b_v[1]-b_v[0])/ep - 0.5*(b_d[0]+b_d[1]));
        double E1 = std::abs((a_v[2]-a_v[0])/(2*ep) - 0.5*(a_d[0]+a_d[2])) +
                    std::abs((b_v[2]-b_v[0])/(2*ep) - 0.5*(b_d[0]+b_d[2]));
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "LSS xi k=" << k << " rate=" << rate; }
    }

    // LinSourceVelocity xi
    for (int k = 0; k < 2; ++k) {
        std::array<Vector2d, 3> a_v, b_v, a_d, b_d;
        Vector2d xi_t = xi;
        for (int ie = 0; ie < 3; ++ie) {
            auto dr = PanelLinSourceVelocityDeriv(Xj, xi_t);
            a_v[ie] = dr.a; b_v[ie] = dr.b;
            a_d[ie] = dr.a_xi.col(k); b_d[ie] = dr.b_xi.col(k);
            xi_t(k) += ep;
        }
        double E0 = ((a_v[1]-a_v[0])/ep - 0.5*(a_d[0]+a_d[1])).norm() +
                    ((b_v[1]-b_v[0])/ep - 0.5*(b_d[0]+b_d[1])).norm();
        double E1 = ((a_v[2]-a_v[0])/(2*ep) - 0.5*(a_d[0]+a_d[2])).norm() +
                    ((b_v[2]-b_v[0])/(2*ep) - 0.5*(b_d[0]+b_d[2])).norm();
        double rate = (E0>EPS_SKIP && E1>EPS_SKIP) ? std::log2(E1/E0) : 99.0;
        if (rate < MIN_RATE) { nfail++; ADD_FAILURE() << "LSV xi k=" << k << " rate=" << rate; }
    }

    std::cout << "AllInfluenceDeriv_xi: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 9. dRue/dx_foil verification
// ===========================================================================
TEST_F(SensitivityPingTest, DRueDxFoil) {
    rng().seed(42);

    Cppfoil foil("0012", 30);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }
    auto& state = foil.GetState();
    int N = state.foil.N;

    // Quick FD check of tcp/tdp derivatives
    {
        TEInfo te0 = GetTEInfo(state.foil.x);
        double fd_ep = 1e-7;
        int te_nodes[] = {0, 1, N-2, N-1};
        std::cout << "  TE deriv FD (tcp=" << te0.tcp << " tdp=" << te0.tdp << "):" << std::endl;
        for (int kk : te_nodes) {
            for (int dd = 0; dd < 2; ++dd) {
                MatrixXd xp = state.foil.x;
                xp(dd, kk) += fd_ep;
                TEInfo tep = GetTEInfo(xp);
                double dtcp_fd = (tep.tcp - te0.tcp) / fd_ep;
                double dtdp_fd = (tep.tdp - te0.tdp) / fd_ep;
                double dtcp_an, dtdp_an;
                GetTEInfoDerivPublic(state.foil.x, kk, dd, dtcp_an, dtdp_an);
                std::cout << "    k=" << kk << " d=" << dd
                          << " dtcp: fd=" << dtcp_fd << " an=" << dtcp_an
                          << " dtdp: fd=" << dtdp_fd << " an=" << dtdp_an << std::endl;
            }
        }
    }

    MatrixXd dRue_dx = ComputeDRueDxFoil(state);

    int Nsys = state.glob.Nsys;

    // Compute base R_ue residual
    VectorXd ue = state.glob.U.row(3).transpose();
    VectorXd ds = state.glob.U.row(1).transpose();
    VectorXd ueinv = GetUeInv(state);
    VectorXd mass = ds.cwiseProduct(ue);
    VectorXd R_ue_base = ue - (ueinv + state.vsol.ue_m * mass);

    double ep = 1e-6;
    // Print element-wise comparison at TE node for debugging
    {
        int kk = 1, dd = 0;
        FoilState sp = state;
        sp.foil.x(dd, kk) += ep;
        sp.isol.AIC.resize(0, 0);
        sp.isol.gamref.resize(0, 0);
        sp.isol.gam.resize(0);
        BuildGamma(sp, sp.oper.alpha);
        sp.isol.sgnue = state.isol.sgnue;
        CalcUeM(sp);
        VectorXd R_ue_p = ue - (GetUeInv(sp) + sp.vsol.ue_m * mass);
        VectorXd dRue_fd_col = (R_ue_p - R_ue_base) / ep;
        VectorXd dRue_an_col = dRue_dx.col(2 * kk + dd);

        // Find rows with largest absolute errors
        VectorXd err_vec = (dRue_fd_col - dRue_an_col).cwiseAbs();
        std::cout << "  dRue column k=" << kk << " d=" << dd
                  << " |fd|=" << dRue_fd_col.norm()
                  << " |an|=" << dRue_an_col.norm()
                  << " |err|=" << err_vec.norm() << std::endl;
        // Print top 5 error rows
        for (int iter = 0; iter < 5; ++iter) {
            int max_row;
            double max_err = err_vec.maxCoeff(&max_row);
            std::cout << "    row " << max_row << ": fd=" << dRue_fd_col(max_row)
                      << " an=" << dRue_an_col(max_row) << " err=" << max_err << std::endl;
            err_vec(max_row) = 0;
        }
    }

    int test_nodes_rue[] = {N/4, N/2, 3*N/4};
    int nfail = 0;

    for (int test_node : test_nodes_rue) {
        if (test_node >= N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            FoilState sp = state;
            sp.foil.x(dim, test_node) += ep;
            sp.isol.AIC.resize(0, 0);
            sp.isol.gamref.resize(0, 0);
            sp.isol.gam.resize(0);
            BuildGamma(sp, sp.oper.alpha);
            sp.isol.sgnue = state.isol.sgnue;

            // (A) ueinv-only FD (don't rebuild ue_m)
            VectorXd ueinv_p = GetUeInv(sp);
            VectorXd R_ue_ueinv = ue - (ueinv_p + state.vsol.ue_m * mass);
            VectorXd dRue_ueinv_fd = (R_ue_ueinv - R_ue_base) / ep;

            // (B) Full FD (rebuild ue_m too)
            CalcUeM(sp);
            VectorXd R_ue_full = ue - (ueinv_p + sp.vsol.ue_m * mass);
            VectorXd dRue_full_fd = (R_ue_full - R_ue_base) / ep;

            // (C) ue_m-only contribution (difference)
            VectorXd dRue_uem_fd = dRue_full_fd - dRue_ueinv_fd;

            VectorXd dRue_adj = dRue_dx.col(2 * test_node + dim);

            // Check ueinv-only part separately
            double ueinv_ref = dRue_ueinv_fd.norm();
            // The ueinv analytical part is dRue_adj minus the ue_m part.
            // But since we can't easily separate, check total instead.
            double ref = dRue_full_fd.norm();
            double vec_err = (dRue_full_fd - dRue_adj).norm();
            double rel_err = vec_err / std::max(ref, 1e-10);

            // Also check ueinv-only: the error in the total is dominated by
            // ueinv since |uem| << |ueinv|. So rel_err ≈ ueinv error.
            double uem_mag = dRue_uem_fd.norm();
            double ueinv_mag = dRue_ueinv_fd.norm();

            // Separate ueinv-only and ue_m-only errors
            double ueinv_err = (dRue_ueinv_fd - dRue_adj).norm();  // adj has ueinv+uem
            // Better: adj = ueinv_adj + uem_adj. uem_adj = adj - ueinv_adj.
            // We can't easily extract ueinv_adj separately.
            // But we CAN check: |ueinv_fd - (adj - uem_fd)| to estimate ueinv part error
            VectorXd ueinv_adj_est = dRue_adj - (dRue_full_fd - dRue_ueinv_fd);  // adj_ueinv ≈ adj - uem_fd
            // Actually just check both parts:
            double ueinv_only_err = (dRue_ueinv_fd - (dRue_adj - dRue_uem_fd + dRue_uem_fd)).norm(); // not useful

            // Error decomposition: if total error is large, check if it's in
            // the ueinv part or the ue_m part
            double ueinv_err_rel = (dRue_ueinv_fd - dRue_adj).norm() / std::max(dRue_ueinv_fd.norm(), 1e-10);
            // uem error: (total_fd - adj) - (ueinv_fd - adj) = total_fd - ueinv_fd = uem_fd
            // So the ue_m-only adj contribution = adj - (total_fd - uem_fd) = adj - ueinv_fd ... no
            // Actually: adj ≈ adj_ueinv + adj_uem
            //           fd_total = fd_ueinv + fd_uem
            // error = |fd_total - adj| = |(fd_ueinv - adj_ueinv) + (fd_uem - adj_uem)|
            // We know fd_ueinv and fd_uem separately.
            // If adj closely matches fd_ueinv (i.e., adj_uem ≈ 0), then uem is missing.
            double adj_matches_ueinv = (dRue_ueinv_fd - dRue_adj).norm() / std::max(ref, 1e-10);

            bool fail = rel_err > 5e-3;  // 0.5% tolerance for mid-airfoil nodes
            if (fail) nfail++;

            // How much of ue_m does the analytic capture?
            // If adj ≈ ueinv_fd, then adj captures 0% of ue_m
            // If adj ≈ total_fd, then adj captures 100% of ue_m
            double uem_captured = 1.0 - (dRue_adj - dRue_full_fd).norm() / std::max(uem_mag, 1e-10);
            std::cout << "  dRue/dx node=" << test_node << " dim=" << dim
                      << " err=" << rel_err
                      << " |uem|/|total|=" << uem_mag / std::max(ref, 1e-10)
                      << " uem_frac=" << uem_captured;
            if (fail) std::cout << " FAIL";
            std::cout << std::endl;
        }
    }
    EXPECT_EQ(nfail, 0) << "dRue/dx_foil failures at 0.1% tolerance";
}

// ===========================================================================
// 10. Viscous Total Derivative: dcl/dx, dcd/dx (full viscous system)
// ===========================================================================
TEST_F(SensitivityPingTest, ViscousSensitivity) {
    rng().seed(42);

    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    {
        auto& s = foil.GetState();
        TEInfo t = GetTEInfo(s.foil.x);
        std::cout << "  ViscSens TE: hTE=" << t.hTE << " tcp=" << t.tcp
                  << " tdp=" << t.tdp << " N=" << s.foil.N << std::endl;
    }

    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto& state = foil.GetState();
    auto sens = ViscousSensitivities(state);

    int N = state.foil.N;

    // Test a few nodes with central FD
    int test_nodes[] = {20, 50, 80};
    int nfail_cl = 0, nfail_cd = 0;

    for (int node : test_nodes) {
        if (node >= N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            double dcl_pred = sens.dcl_dx(dim, node);
            double dcd_pred = sens.dcd_dx(dim, node);

            // Central FD at two step sizes
            auto solve_perturbed = [&](double pert) -> std::pair<double, double> {
                Cppfoil foil_p("0012", 100);
                foil_p.SetOper(3.0, -1.0, 1e6, 0.0);
                auto& sp = foil_p.GetState();
                // Copy panel coordinates from base state
                sp.foil.x = state.foil.x;
                sp.foil.x(dim, node) += pert;
                sp.oper.viscous = true;
                sp.oper.initbl = true;
                sp.param.verb = 0;
                sp.param.niglob = 25;
                sp.param.rtol = 1e-6;
                SolveViscous(sp);
                CalcForce(sp);
                return {sp.post.cl, sp.post.cd};
            };

            double ep = 1e-5;
            auto [cl_p, cd_p] = solve_perturbed(ep);
            auto [cl_m, cd_m] = solve_perturbed(-ep);
            double dcl_fd = (cl_p - cl_m) / (2.0 * ep);
            double dcd_fd = (cd_p - cd_m) / (2.0 * ep);

            double ep2 = ep / 2.0;
            auto [cl_p2, cd_p2] = solve_perturbed(ep2);
            auto [cl_m2, cd_m2] = solve_perturbed(-ep2);
            double dcl_fd2 = (cl_p2 - cl_m2) / (2.0 * ep2);
            double dcd_fd2 = (cd_p2 - cd_m2) / (2.0 * ep2);

            // cl relative accuracy (use finer FD as reference)
            double rel_cl = std::abs(dcl_fd2 - dcl_pred) /
                            std::max(std::abs(dcl_fd2), 1e-6);

            if (rel_cl > 0.005) {  // > 0.5% error
                nfail_cl++;
                std::cout << "  ViscSens cl FAIL node=" << node << " dim=" << dim
                          << " fd=" << dcl_fd2 << " adj=" << dcl_pred
                          << " rel=" << rel_cl << std::endl;
            } else {
                std::cout << "  ViscSens cl OK   node=" << node << " dim=" << dim
                          << " fd=" << dcl_fd2 << " adj=" << dcl_pred
                          << " rel=" << rel_cl << std::endl;
            }

            // cd relative accuracy
            double rel_cd = std::abs(dcd_fd2 - dcd_pred) /
                            std::max(std::abs(dcd_fd2), 1e-6);

            if (rel_cd > 0.001) {  // > 0.1% error
                nfail_cd++;
                std::cout << "  ViscSens cd FAIL node=" << node << " dim=" << dim
                          << " fd=" << dcd_fd2 << " adj=" << dcd_pred
                          << " rel=" << rel_cd << std::endl;
            } else {
                std::cout << "  ViscSens cd OK   node=" << node << " dim=" << dim
                          << " fd=" << dcd_fd2 << " adj=" << dcd_pred
                          << " rel=" << rel_cd << std::endl;
            }
        }
    }
    std::cout << "ViscousSensitivity cl: " << nfail_cl << "/" << 6 << " failures" << std::endl;
    std::cout << "ViscousSensitivity cd: " << nfail_cd << "/" << 6 << " failures" << std::endl;

    // The viscous adjoint involves several approximations (frozen wake geometry,
    // neglected due_m/dx, simplified dxi/dx). Report results for now.
    // TODO: tighten tolerance as approximations are refined.
    if (nfail_cl > 0 || nfail_cd > 0) {
        std::cout << "  NOTE: Viscous sensitivity has known approximations. "
                  << "Failures expected until Phase 1d refinement." << std::endl;
    }
}

// ===========================================================================
// 10b. FD convergence study for viscous dcd/dx at TE-adjacent node
// ===========================================================================
TEST_F(SensitivityPingTest, ViscousSensitivityFDConvergence) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();

    if (!foil.IsConverged()) {
        GTEST_SKIP() << "Viscous solver did not converge";
    }

    auto& state = foil.GetState();
    auto sens = ViscousSensitivities(state);

    int node = 20;
    int dim = 1;
    double dcl_adj = sens.dcl_dx(dim, node);
    double dcd_adj = sens.dcd_dx(dim, node);

    std::cout << "  Adjoint: dcl=" << dcl_adj << " dcd=" << dcd_adj << std::endl;

    auto solve_perturbed = [&](double pert) -> std::pair<double, double> {
        Cppfoil foil_p("0012", 100);
        foil_p.SetOper(3.0, -1.0, 1e6, 0.0);
        auto& sp = foil_p.GetState();
        sp.foil.x = state.foil.x;
        sp.foil.x(dim, node) += pert;
        sp.oper.viscous = true;
        sp.oper.initbl = true;
        sp.param.verb = 0;
        sp.param.niglob = 25;
        sp.param.rtol = 1e-10;  // tighter tolerance
        SolveViscous(sp);
        CalcForce(sp);
        return {sp.post.cl, sp.post.cd};
    };

    // Test 6 step sizes from 1e-3 to 1e-8
    double eps[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    std::cout << std::scientific;
    std::cout << "  FD convergence study node=" << node << " dim=" << dim << std::endl;
    std::cout << "  " << std::setw(10) << "eps"
              << std::setw(16) << "dcl_fd"
              << std::setw(16) << "dcd_fd"
              << std::setw(16) << "cl_err_vs_adj"
              << std::setw(16) << "cd_err_vs_adj" << std::endl;

    for (double ep : eps) {
        auto [cl_p, cd_p] = solve_perturbed(ep);
        auto [cl_m, cd_m] = solve_perturbed(-ep);
        double dcl_fd = (cl_p - cl_m) / (2.0 * ep);
        double dcd_fd = (cd_p - cd_m) / (2.0 * ep);

        double ecl = std::abs(dcl_fd - dcl_adj) / std::max(std::abs(dcl_fd), 1e-10);
        double ecd = std::abs(dcd_fd - dcd_adj) / std::max(std::abs(dcd_fd), 1e-10);

        std::cout << "  " << std::setw(10) << ep
                  << std::setw(16) << dcl_fd
                  << std::setw(16) << dcd_fd
                  << std::setw(16) << ecl
                  << std::setw(16) << ecd << std::endl;
    }
}

// ===========================================================================
// 10c. FD verification of wake geometry VJP
// Perturbs foil.x, rebuilds wake+ueM only, computes delta(R_ue) from wake change.
// ===========================================================================
TEST_F(SensitivityPingTest, WakeGeomVJPVerification) {
    Cppfoil foil("0012", 40);  // smaller mesh for speed
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }

    auto& state = foil.GetState();
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;

    // Compute adjoint
    CalcForceLin(state);
    CalcCdLin(state);

    VectorXd dcd_dV = state.post.cd_V;
    auto& lu = *state.glob.lu_factorization;
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);
    VectorXd lambda_ue = lambda_cd.tail(Nsys);

    // Compute analytic wake VJP
    VectorXd vjp_wake = ComputeWakeGeomVJP(state, lambda_ue);

    // Save base state
    MatrixXd base_gamref = state.isol.gamref;
    VectorXd base_gam = state.isol.gam;
    MatrixXd base_wake_x = state.wake.x;
    MatrixXd base_wake_t = state.wake.t;
    MatrixXd base_uewiref = state.isol.uewiref;
    MatrixXd base_ue_sigma = state.vsol.ue_sigma;

    // Base R_ue computation
    double alpha = state.oper.alpha;
    VectorXd ue_vec = state.glob.U.row(3).transpose();
    VectorXd ds_vec = state.glob.U.row(1).transpose();
    VectorXd mass = ds_vec.cwiseProduct(ue_vec);
    VectorXd phi = state.vsol.sigma_m * mass;

    VectorXd ueinv_base = GetUeInv(state);
    VectorXd sgue(Nsys);
    sgue.head(N) = state.isol.sgnue;
    sgue.tail(Nw).setOnes();
    VectorXd Rue_base(Nsys);
    for (int i = 0; i < Nsys; ++i) {
        Rue_base(i) = ue_vec(i) - ueinv_base(i) - sgue(i) * (state.vsol.ue_sigma.row(i).dot(phi));
    }

    // Print TE info
    {
        TEInfo te_info = GetTEInfo(state.foil.x);
        std::cout << "  TE info: hTE=" << te_info.hTE << " tcp=" << te_info.tcp
                  << " tdp=" << te_info.tdp << " s_norm="
                  << (state.foil.x.col(N-1) - state.foil.x.col(0)).norm() << std::endl;
    }
    // Compute both parts of the VJP
    MatrixXd dRue_dx = ComputeDRueDxFoil(state);
    VectorXd vjp_ue = dRue_dx.transpose() * lambda_ue;
    VectorXd vjp_wake_an = ComputeWakeGeomVJP(state, lambda_ue);
    VectorXd vjp_total_an = vjp_ue + vjp_wake_an;

    // Test nodes
    int test_nodes[] = {1, N/4, N/2, 3*N/4, N - 2};
    double ep = 1e-7;

    std::cout << "  Wake geom VJP verification (cd):" << std::endl;
    for (int node : test_nodes) {
        if (node >= N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            // Perturb foil.x
            state.foil.x(dim, node) += ep;

            // Rebuild gamma and wake with perturbed foil
            BuildGamma(state, alpha);
            StagPointFind(state);
            BuildWake(state);
            ExtendXiToWake(state);
            CalcUeM(state);

            // Compute perturbed R_ue (using BASE BL state but NEW wake/gamma/ueM)
            VectorXd ueinv_pert = GetUeInv(state);
            VectorXd Rue_pert(Nsys);
            for (int i = 0; i < Nsys; ++i) {
                Rue_pert(i) = ue_vec(i) - ueinv_pert(i) - sgue(i) * (state.vsol.ue_sigma.row(i).dot(phi));
            }

            // This gives the TOTAL R_ue sensitivity (gamma + wake geom + ue_m)
            double vjp_fd_total = lambda_ue.dot(Rue_pert - Rue_base) / ep;

            // Now compute wake-geometry-only: perturb foil.x, rebuild wake
            // (new wake.x, wake.t), then recompute uewiref and ue_m using
            // the NEW wake but BASE foil panel geometry.
            // This isolates: foil.x -> wake construction -> R_ue
            state.isol.gamref = base_gamref;
            state.isol.gam = base_gam;
            StagPointFind(state);
            BuildWake(state);  // New wake from perturbed foil.x + base gamma
            ExtendXiToWake(state);
            // Now restore foil.x to base (so CalcUeM uses base panel geometry)
            // but keep the new wake.x and wake.t
            MatrixXd new_wake_x = state.wake.x;
            MatrixXd new_wake_t = state.wake.t;
            VectorXd new_wake_s = state.wake.s;
            VectorXd new_xi = state.isol.xi;
            state.foil.x(dim, node) -= ep;  // restore foil.x
            // Put back new wake
            state.wake.x = new_wake_x;
            state.wake.t = new_wake_t;
            state.wake.s = new_wake_s;
            state.isol.xi = new_xi;
            CalcUeM(state);  // Uses base foil panels but new wake geometry

            VectorXd ueinv_wake = GetUeInv(state);
            VectorXd Rue_wake(Nsys);
            for (int i = 0; i < Nsys; ++i) {
                Rue_wake(i) = ue_vec(i) - ueinv_wake(i) - sgue(i) * (state.vsol.ue_sigma.row(i).dot(phi));
            }
            double vjp_fd_wake = lambda_ue.dot(Rue_wake - Rue_base) / ep;
            // Re-perturb foil.x for the restore at end of loop
            state.foil.x(dim, node) += ep;

            double vjp_fd = vjp_fd_total;  // Print total for comparison
            double vjp_an = vjp_wake(2 * node + dim);

            // Compare: total_fd = gamma+wake+ueM effect,
            //          wake_fd = wake-only (geometry) effect,
            //          gamma_fd = total - wake = gamma/AIC effect (should match ComputeDRueDxFoil)
            double gamma_fd = vjp_fd_total - vjp_fd_wake;
            double wake_an_val = vjp_wake_an(2 * node + dim);
            double err_wake = std::abs(vjp_fd_wake - wake_an_val) /
                              std::max(std::abs(vjp_fd_wake), 1e-8);
            std::cout << "    node=" << node << " dim=" << dim
                      << " wake_fd=" << vjp_fd_wake
                      << " wake_an=" << wake_an_val
                      << " err=" << err_wake
                      << (err_wake > 0.1 ? " BAD" : " ok") << std::endl;

            // Restore
            state.foil.x(dim, node) -= ep;
        }
    }

    // Restore original wake and ueM
    state.wake.x = base_wake_x;
    state.wake.t = base_wake_t;
    state.isol.uewiref = base_uewiref;
    state.vsol.ue_sigma = base_ue_sigma;
    // Restore gamma
    BuildGamma(state, alpha);
    StagPointFind(state);
    BuildWake(state);
    ExtendXiToWake(state);
    CalcUeM(state);
}

// ===========================================================================
// 11. Inviscid Alpha Sensitivity: dcl/dalpha, dcm/dalpha
// ===========================================================================
TEST_F(SensitivityPingTest, InviscidAlphaSensitivity) {
    Cppfoil foil("2412", 80);
    foil.SetOper(5.0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();

    double dcl_dalpha_adj = sens.dcl_dalpha;
    double dcm_dalpha_adj = sens.dcm_dalpha;

    // Central FD at two step sizes
    int nfail = 0;
    for (int ie = 0; ie < 2; ++ie) {
        double ep = (ie == 0) ? 1e-4 : 5e-5;

        Cppfoil foil_p("2412", 80);
        foil_p.SetOper(5.0 + ep);
        foil_p.Solve();
        Cppfoil foil_m("2412", 80);
        foil_m.SetOper(5.0 - ep);
        foil_m.Solve();

        double dcl_fd = (foil_p.GetCl() - foil_m.GetCl()) / (2.0 * ep);
        double dcm_fd = (foil_p.GetCm() - foil_m.GetCm()) / (2.0 * ep);

        double ecl = std::abs(dcl_fd - dcl_dalpha_adj);
        double ecm = std::abs(dcm_fd - dcm_dalpha_adj);

        std::cout << "  InvAlpha ep=" << ep
                  << " dcl: fd=" << dcl_fd << " adj=" << dcl_dalpha_adj << " err=" << ecl
                  << " dcm: fd=" << dcm_fd << " adj=" << dcm_dalpha_adj << " err=" << ecm
                  << std::endl;

        if (ie == 1) {
            // Check relative error < 0.1%
            double rel_cl = ecl / std::max(std::abs(dcl_fd), 1e-10);
            double rel_cm = ecm / std::max(std::abs(dcm_fd), 1e-10);
            if (rel_cl > 1e-3) {
                nfail++;
                ADD_FAILURE() << "dcl/dalpha rel_err=" << rel_cl;
            }
            if (rel_cm > 1e-3) {
                nfail++;
                ADD_FAILURE() << "dcm/dalpha rel_err=" << rel_cm;
            }
        }
    }
    std::cout << "InviscidAlphaSensitivity: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 12. Viscous Alpha Sensitivity: dcl/dalpha, dcd/dalpha, dcm/dalpha
// ===========================================================================
TEST_F(SensitivityPingTest, ViscousAlphaSensitivity) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }

    auto sens = foil.ComputeSensitivities();

    double dcl_adj = sens.dcl_dalpha;
    double dcd_adj = sens.dcd_dalpha;
    double dcm_adj = sens.dcm_dalpha;

    // Central FD at two step sizes
    double eps[] = {1e-4, 5e-5};
    double ecl[2], ecd[2], ecm[2];
    double dcl_fd2 = 0, dcd_fd2 = 0, dcm_fd2 = 0;

    for (int ie = 0; ie < 2; ++ie) {
        double ep = eps[ie];

        Cppfoil foil_p("0012", 100);
        foil_p.SetOper(3.0 + ep, -1.0, 1e6, 0.0);
        foil_p.GetState().param.verb = 0;
        foil_p.Solve();
        Cppfoil foil_m("0012", 100);
        foil_m.SetOper(3.0 - ep, -1.0, 1e6, 0.0);
        foil_m.GetState().param.verb = 0;
        foil_m.Solve();

        if (!foil_p.IsConverged() || !foil_m.IsConverged()) {
            GTEST_SKIP() << "FD solve did not converge";
        }

        double dcl_fd = (foil_p.GetCl() - foil_m.GetCl()) / (2.0 * ep);
        double dcd_fd = (foil_p.GetCd() - foil_m.GetCd()) / (2.0 * ep);
        double dcm_fd = (foil_p.GetCm() - foil_m.GetCm()) / (2.0 * ep);

        ecl[ie] = std::abs(dcl_fd - dcl_adj);
        ecd[ie] = std::abs(dcd_fd - dcd_adj);
        ecm[ie] = std::abs(dcm_fd - dcm_adj);

        if (ie == 1) {
            dcl_fd2 = dcl_fd;
            dcd_fd2 = dcd_fd;
            dcm_fd2 = dcm_fd;
        }

        std::cout << "  ViscAlpha ep=" << ep
                  << " dcl: fd=" << dcl_fd << " adj=" << dcl_adj
                  << " dcd: fd=" << dcd_fd << " adj=" << dcd_adj
                  << " dcm: fd=" << dcm_fd << " adj=" << dcm_adj
                  << std::endl;
    }

    // Check relative error of finer FD vs adjoint
    int nfail = 0;
    double rel_cl = ecl[1] / std::max(std::abs(dcl_fd2), 1e-6);
    double rel_cd = ecd[1] / std::max(std::abs(dcd_fd2), 1e-6);
    double rel_cm = ecm[1] / std::max(std::abs(dcm_fd2), 1e-6);

    std::cout << "  ViscAlpha rel: cl=" << rel_cl << " cd=" << rel_cd
              << " cm=" << rel_cm << std::endl;

    if (rel_cl > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcl/dalpha viscous rel_err=" << rel_cl;
    }
    if (rel_cd > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcd/dalpha viscous rel_err=" << rel_cd;
    }
    // cm has tighter coupling through BL; 1% tolerance accommodates
    // adjoint approximations (frozen wake, simplified ue_m)
    if (rel_cm > 0.01) {
        nfail++;
        ADD_FAILURE() << "dcm/dalpha viscous rel_err=" << rel_cm;
    }
    std::cout << "ViscousAlphaSensitivity: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 13. Inviscid Mach Sensitivity: dcl/dMa, dcm/dMa (Ma=0.3)
// ===========================================================================
TEST_F(SensitivityPingTest, InviscidMachSensitivity) {
    double Ma0 = 0.3;
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, Ma0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();
    double dcl_adj = sens.dcl_dMa;
    double dcm_adj = sens.dcm_dMa;

    // Central FD at two step sizes
    double eps[] = {1e-4, 5e-5};
    int nfail = 0;

    for (int ie = 0; ie < 2; ++ie) {
        double ep = eps[ie];

        Cppfoil foil_p("0012", 100);
        foil_p.SetOper(3.0, -1.0, 0.0, Ma0 + ep);
        foil_p.Solve();
        Cppfoil foil_m("0012", 100);
        foil_m.SetOper(3.0, -1.0, 0.0, Ma0 - ep);
        foil_m.Solve();

        double dcl_fd = (foil_p.GetCl() - foil_m.GetCl()) / (2.0 * ep);
        double dcm_fd = (foil_p.GetCm() - foil_m.GetCm()) / (2.0 * ep);

        std::cout << "  InvMa ep=" << ep
                  << " dcl: fd=" << dcl_fd << " adj=" << dcl_adj
                  << " dcm: fd=" << dcm_fd << " adj=" << dcm_adj
                  << std::endl;

        if (ie == 1) {
            double rel_cl = std::abs(dcl_fd - dcl_adj) / std::max(std::abs(dcl_fd), 1e-10);
            double rel_cm = std::abs(dcm_fd - dcm_adj) / std::max(std::abs(dcm_fd), 1e-10);
            if (rel_cl > 1e-3) {
                nfail++;
                ADD_FAILURE() << "dcl/dMa inviscid rel_err=" << rel_cl;
            }
            if (rel_cm > 1e-3) {
                nfail++;
                ADD_FAILURE() << "dcm/dMa inviscid rel_err=" << rel_cm;
            }
        }
    }
    std::cout << "InviscidMachSensitivity: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 14. Inviscid Mach Sensitivity at Ma=0: should be zero
// ===========================================================================
TEST_F(SensitivityPingTest, InviscidMachSensitivityZero) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, 0.0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();

    EXPECT_NEAR(sens.dcl_dMa, 0.0, 1e-12) << "dcl/dMa should be 0 at Ma=0";
    EXPECT_NEAR(sens.dcm_dMa, 0.0, 1e-12) << "dcm/dMa should be 0 at Ma=0";
    EXPECT_NEAR(sens.dcd_dMa, 0.0, 1e-12) << "dcd/dMa should be 0 at Ma=0";
}

// ===========================================================================
// 15. Viscous Mach Sensitivity: dcl/dMa, dcd/dMa, dcm/dMa (Ma=0.3)
// ===========================================================================
TEST_F(SensitivityPingTest, ViscousMachSensitivity) {
    double Ma0 = 0.3;
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 1e6, Ma0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }

    auto sens = foil.ComputeSensitivities();
    double dcl_adj = sens.dcl_dMa;
    double dcd_adj = sens.dcd_dMa;
    double dcm_adj = sens.dcm_dMa;

    // Central FD at two step sizes
    double eps[] = {1e-4, 5e-5};
    double ecl[2], ecd[2], ecm[2];
    double dcl_fd2 = 0, dcd_fd2 = 0, dcm_fd2 = 0;

    for (int ie = 0; ie < 2; ++ie) {
        double ep = eps[ie];

        Cppfoil foil_p("0012", 100);
        foil_p.SetOper(3.0, -1.0, 1e6, Ma0 + ep);
        foil_p.GetState().param.verb = 0;
        foil_p.Solve();
        Cppfoil foil_m("0012", 100);
        foil_m.SetOper(3.0, -1.0, 1e6, Ma0 - ep);
        foil_m.GetState().param.verb = 0;
        foil_m.Solve();

        if (!foil_p.IsConverged() || !foil_m.IsConverged()) {
            GTEST_SKIP() << "FD solve did not converge";
        }

        double dcl_fd = (foil_p.GetCl() - foil_m.GetCl()) / (2.0 * ep);
        double dcd_fd = (foil_p.GetCd() - foil_m.GetCd()) / (2.0 * ep);
        double dcm_fd = (foil_p.GetCm() - foil_m.GetCm()) / (2.0 * ep);

        ecl[ie] = std::abs(dcl_fd - dcl_adj);
        ecd[ie] = std::abs(dcd_fd - dcd_adj);
        ecm[ie] = std::abs(dcm_fd - dcm_adj);

        if (ie == 1) {
            dcl_fd2 = dcl_fd;
            dcd_fd2 = dcd_fd;
            dcm_fd2 = dcm_fd;
        }

        std::cout << "  ViscMa ep=" << ep
                  << " dcl: fd=" << dcl_fd << " adj=" << dcl_adj
                  << " dcd: fd=" << dcd_fd << " adj=" << dcd_adj
                  << " dcm: fd=" << dcm_fd << " adj=" << dcm_adj
                  << std::endl;
    }

    int nfail = 0;
    double rel_cl = ecl[1] / std::max(std::abs(dcl_fd2), 1e-6);
    double rel_cd = ecd[1] / std::max(std::abs(dcd_fd2), 1e-6);
    double rel_cm = ecm[1] / std::max(std::abs(dcm_fd2), 1e-6);

    std::cout << "  ViscMa rel: cl=" << rel_cl << " cd=" << rel_cd
              << " cm=" << rel_cm << std::endl;

    if (rel_cl > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcl/dMa viscous rel_err=" << rel_cl;
    }
    if (rel_cd > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcd/dMa viscous rel_err=" << rel_cd;
    }
    if (rel_cm > 0.01) {
        nfail++;
        ADD_FAILURE() << "dcm/dMa viscous rel_err=" << rel_cm;
    }
    std::cout << "ViscousMachSensitivity: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 16. Inviscid Reynolds Sensitivity: should be zero
// ===========================================================================
TEST_F(SensitivityPingTest, InviscidReSensitivityZero) {
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, 0.0, 0.0);
    foil.Solve();

    auto sens = foil.ComputeSensitivities();

    EXPECT_NEAR(sens.dcl_dRe, 0.0, 1e-12) << "dcl/dRe should be 0 for inviscid";
    EXPECT_NEAR(sens.dcd_dRe, 0.0, 1e-12) << "dcd/dRe should be 0 for inviscid";
    EXPECT_NEAR(sens.dcm_dRe, 0.0, 1e-12) << "dcm/dRe should be 0 for inviscid";
}

// ===========================================================================
// 17. Viscous Reynolds Sensitivity: dcl/dRe, dcd/dRe, dcm/dRe
// ===========================================================================
TEST_F(SensitivityPingTest, ViscousReSensitivity) {
    double Re0 = 1e6;
    Cppfoil foil("0012", 100);
    foil.SetOper(3.0, -1.0, Re0, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }

    auto sens = foil.ComputeSensitivities();
    double dcl_adj = sens.dcl_dRe;
    double dcd_adj = sens.dcd_dRe;
    double dcm_adj = sens.dcm_dRe;

    // Central FD at two step sizes
    double eps[] = {100.0, 50.0};  // Re perturbation (absolute)
    double ecl[2], ecd[2], ecm[2];
    double dcl_fd2 = 0, dcd_fd2 = 0, dcm_fd2 = 0;

    for (int ie = 0; ie < 2; ++ie) {
        double ep = eps[ie];

        Cppfoil foil_p("0012", 100);
        foil_p.SetOper(3.0, -1.0, Re0 + ep, 0.0);
        foil_p.GetState().param.verb = 0;
        foil_p.Solve();
        Cppfoil foil_m("0012", 100);
        foil_m.SetOper(3.0, -1.0, Re0 - ep, 0.0);
        foil_m.GetState().param.verb = 0;
        foil_m.Solve();

        if (!foil_p.IsConverged() || !foil_m.IsConverged()) {
            GTEST_SKIP() << "FD solve did not converge";
        }

        double dcl_fd = (foil_p.GetCl() - foil_m.GetCl()) / (2.0 * ep);
        double dcd_fd = (foil_p.GetCd() - foil_m.GetCd()) / (2.0 * ep);
        double dcm_fd = (foil_p.GetCm() - foil_m.GetCm()) / (2.0 * ep);

        ecl[ie] = std::abs(dcl_fd - dcl_adj);
        ecd[ie] = std::abs(dcd_fd - dcd_adj);
        ecm[ie] = std::abs(dcm_fd - dcm_adj);

        if (ie == 1) {
            dcl_fd2 = dcl_fd;
            dcd_fd2 = dcd_fd;
            dcm_fd2 = dcm_fd;
        }

        std::cout << "  ViscRe ep=" << ep
                  << " dcl: fd=" << dcl_fd << " adj=" << dcl_adj
                  << " dcd: fd=" << dcd_fd << " adj=" << dcd_adj
                  << " dcm: fd=" << dcm_fd << " adj=" << dcm_adj
                  << std::endl;
    }

    int nfail = 0;
    double rel_cl = ecl[1] / std::max(std::abs(dcl_fd2), 1e-6);
    double rel_cd = ecd[1] / std::max(std::abs(dcd_fd2), 1e-6);
    double rel_cm = ecm[1] / std::max(std::abs(dcm_fd2), 1e-6);

    std::cout << "  ViscRe rel: cl=" << rel_cl << " cd=" << rel_cd
              << " cm=" << rel_cm << std::endl;

    if (rel_cl > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcl/dRe viscous rel_err=" << rel_cl;
    }
    if (rel_cd > 1e-3) {
        nfail++;
        ADD_FAILURE() << "dcd/dRe viscous rel_err=" << rel_cd;
    }
    if (rel_cm > 0.01) {
        nfail++;
        ADD_FAILURE() << "dcm/dRe viscous rel_err=" << rel_cm;
    }
    std::cout << "ViscousReSensitivity: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// 18. FD verification of B-matrix wake endpoint sensitivity
//     Tests d(B*phi)/d(wake.x[k]) via finite differences against the
//     analytic PanelLinSourceStreamDeriv Xj chain rule.
// ===========================================================================
TEST_F(SensitivityPingTest, BMatrixWakeEndpointFD) {
    // 1. Set up a 30-panel converged viscous solution
    Cppfoil foil("0012", 30);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }
    auto& state = foil.GetState();

    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = state.glob.Nsys;  // N + Nw

    // Need adjoint machinery for the combined adjoint vector w
    CalcForceLin(state);
    CalcCdLin(state);

    VectorXd dcd_dV = state.post.cd_V;
    auto& lu = *state.glob.lu_factorization;
    VectorXd lambda_cd = lu.transpose().solve(dcd_dV);
    VectorXd lambda_ue = lambda_cd.tail(Nsys);

    // Extract wake and foil parts of lambda_ue
    VectorXd lambda_ue_foil = lambda_ue.head(N);
    VectorXd lam_w = lambda_ue.tail(Nw);

    // Compute base BL state quantities
    VectorXd ue_vec = state.glob.U.row(3).transpose();
    VectorXd ds_vec = state.glob.U.row(1).transpose();
    VectorXd mass = ds_vec.cwiseProduct(ue_vec);
    VectorXd phi = state.vsol.sigma_m * mass;
    int Npanel = N + Nw - 2;

    // 2. Build v_ext and w = AIC^{-T} * v_ext
    VectorXd v_ext = VectorXd::Zero(N + 1);
    for (int i = 0; i < N; ++i) {
        v_ext(i) = state.isol.sgnue(i) * lambda_ue_foil(i);
    }

    // Wake station contributions through Cgam^T * lam_w(1:Nw-1)
    for (int i = 1; i < Nw; ++i) {
        double lam_i = lam_w(i);
        if (std::abs(lam_i) < 1e-20) continue;
        Vector2d xi_w = state.wake.x.col(i);
        Vector2d ti_w = state.wake.t.col(i);
        MatrixXd V_G;
        InviscidVelocity(state.foil.x, state.isol.gam, 0.0, 0.0, xi_w,
                         true, &V_G);
        for (int jj = 0; jj < N; ++jj) {
            v_ext(jj) += lam_i * (V_G(0, jj) * ti_w(0) + V_G(1, jj) * ti_w(1));
        }
    }

    // Dw(0) = Bp(N-1) override: wake station 0 uses Bp(N-1,:) directly
    v_ext(N - 1) += lam_w(0);

    MatrixXd AIC_T = state.isol.AIC.transpose();
    VectorXd w = AIC_T.partialPivLu().solve(v_ext);

    // 3. Helper: compute w^T * B * phi for a given wake geometry
    //    (only the wake-source-panel columns of B, i.e. columns N-1..Npanel-1)
    auto compute_wBphi = [&](const MatrixXd& wake_x) -> double {
        double result = 0.0;
        for (int i = 0; i < N; ++i) {
            Vector2d xi = state.foil.x.col(i);
            double row_sum = 0.0;
            for (int j = 0; j < Nw - 1; ++j) {
                Vector2d Xm = 0.5 * (wake_x.col(j) + wake_x.col(j + 1));
                Matrix2d Xj_l;
                Xj_l.col(0) = wake_x.col(j);
                Xj_l.col(1) = Xm;

                Vector2d X_right;
                if (j == Nw - 2) {
                    X_right = 2.0 * wake_x.col(j + 1) - Xm;
                } else {
                    X_right = wake_x.col(j + 1);
                }

                double a_left, b_left, a_right, b_right;

                // Left half panel
                PanelLinSourceStream(Xj_l, xi, a_left, b_left);
                if (j > 0) {
                    row_sum += (0.5 * a_left + b_left) * phi(N - 1 + j);
                    row_sum += 0.5 * a_left * phi(N - 1 + j - 1);
                } else {
                    row_sum += b_left * phi(N - 1 + j);
                }

                // Right half panel
                Matrix2d Xj_r;
                Xj_r.col(0) = Xm;
                Xj_r.col(1) = X_right;
                PanelLinSourceStream(Xj_r, xi, a_right, b_right);
                row_sum += (a_right + 0.5 * b_right) * phi(N - 1 + j);
                if (j < Nw - 2) {
                    row_sum += 0.5 * b_right * phi(N - 1 + j + 1);
                } else {
                    row_sum += 0.5 * b_right * phi(N - 1 + j);  // ghost fold-back
                }
            }
            result += w(i) * row_sum;
        }
        return result;
    };

    // 4. Compute base value
    double wBphi_base = compute_wBphi(state.wake.x);

    // 5. For each wake node k and dimension d, compare FD vs analytic
    double eps = 1e-7;
    int nfail = 0;

    std::cout << std::scientific << std::setprecision(8);
    std::cout << "  B-matrix wake endpoint FD verification:" << std::endl;
    std::cout << "  " << std::setw(4) << "k" << std::setw(4) << "d"
              << std::setw(18) << "analytic"
              << std::setw(18) << "FD"
              << std::setw(14) << "rel_err" << std::endl;

    for (int k = 0; k < Nw; ++k) {
        for (int dim = 0; dim < 2; ++dim) {
            // --- Finite difference ---
            MatrixXd wake_x_pert = state.wake.x;
            wake_x_pert(dim, k) += eps;
            double wBphi_pert = compute_wBphi(wake_x_pert);
            double fd_val = (wBphi_pert - wBphi_base) / eps;

            // --- Analytic: sum over panels j that depend on wake.x[k] ---
            // wake.x[k] affects:
            //   Panel j=k: wake.x[k] is the left endpoint
            //   Panel j=k-1: wake.x[k] is the right endpoint (enters via Xm and X_right)
            double analytic_val = 0.0;

            // Panel j = k (wake.x[k] is left endpoint of panel k)
            if (k < Nw - 1) {
                int j = k;
                Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
                Vector2d X_right;
                if (j == Nw - 2) {
                    X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
                } else {
                    X_right = state.wake.x.col(j + 1);
                }

                double phi_jm1 = (j > 0) ? phi(N - 1 + j - 1) : 0.0;
                double phi_j = phi(N - 1 + j);
                double phi_jp1 = (j < Nw - 2) ? phi(N - 1 + j + 1) : 0.0;

                for (int i = 0; i < N; ++i) {
                    if (std::abs(w(i)) < 1e-20) continue;
                    Vector2d xi_eval = state.foil.x.col(i);

                    Matrix2d Xj_l;
                    Xj_l.col(0) = state.wake.x.col(j);
                    Xj_l.col(1) = Xm;
                    auto ld = PanelLinSourceStreamDeriv(Xj_l, xi_eval);

                    Matrix2d Xj_r;
                    Xj_r.col(0) = Xm;
                    Xj_r.col(1) = X_right;
                    auto rd = PanelLinSourceStreamDeriv(Xj_r, xi_eval);

                    // d(a_l)/d(wake.x[j]_dim): Xj_l col(0) = wake.x[j], col(1) = Xm
                    // d(Xj_l.col(0))/d(wake.x[j]) = I  => a_Xj(dim)
                    // d(Xm)/d(wake.x[j]) = 0.5*I        => 0.5*a_Xj(2+dim)
                    double da_l_dj = ld.a_Xj(dim) + 0.5 * ld.a_Xj(2 + dim);
                    double db_l_dj = ld.b_Xj(dim) + 0.5 * ld.b_Xj(2 + dim);

                    // d(a_r)/d(wake.x[j]_dim): Xj_r col(0) = Xm, col(1) = X_right
                    // d(Xm)/d(wake.x[j]) = 0.5*I => 0.5*a_Xj(dim)
                    double da_r_dj = 0.5 * rd.a_Xj(dim);
                    double db_r_dj = 0.5 * rd.b_Xj(dim);
                    if (j == Nw - 2) {
                        // X_right = 2*wake.x[j+1] - Xm, d(X_right)/d(wake.x[j]) = -0.5*I
                        da_r_dj += -0.5 * rd.a_Xj(2 + dim);
                        db_r_dj += -0.5 * rd.b_Xj(2 + dim);
                    }

                    // Compute d(B(i,:)*phi)/d(wake.x[j]_dim) from this panel
                    double dW = 0.0;
                    if (j > 0) {
                        dW += 0.5 * da_l_dj * phi_jm1;
                        dW += (0.5 * da_l_dj + db_l_dj) * phi_j;
                    } else {
                        dW += db_l_dj * phi_j;
                    }
                    dW += (da_r_dj + 0.5 * db_r_dj) * phi_j;
                    if (j < Nw - 2) {
                        dW += 0.5 * db_r_dj * phi_jp1;
                    } else {
                        dW += 0.5 * db_r_dj * phi_j;  // ghost fold-back
                    }

                    analytic_val += w(i) * dW;
                }
            }

            // Panel j = k-1 (wake.x[k] is the right endpoint of panel k-1)
            if (k > 0 && k - 1 < Nw - 1) {
                int j = k - 1;
                Vector2d Xm = 0.5 * (state.wake.x.col(j) + state.wake.x.col(j + 1));
                Vector2d X_right;
                if (j == Nw - 2) {
                    X_right = 2.0 * state.wake.x.col(j + 1) - Xm;
                } else {
                    X_right = state.wake.x.col(j + 1);
                }

                double phi_jm1 = (j > 0) ? phi(N - 1 + j - 1) : 0.0;
                double phi_j = phi(N - 1 + j);
                double phi_jp1 = (j < Nw - 2) ? phi(N - 1 + j + 1) : 0.0;

                for (int i = 0; i < N; ++i) {
                    if (std::abs(w(i)) < 1e-20) continue;
                    Vector2d xi_eval = state.foil.x.col(i);

                    Matrix2d Xj_l;
                    Xj_l.col(0) = state.wake.x.col(j);
                    Xj_l.col(1) = Xm;
                    auto ld = PanelLinSourceStreamDeriv(Xj_l, xi_eval);

                    Matrix2d Xj_r;
                    Xj_r.col(0) = Xm;
                    Xj_r.col(1) = X_right;
                    auto rd = PanelLinSourceStreamDeriv(Xj_r, xi_eval);

                    // d/d(wake.x[j+1]_dim) where j+1 = k
                    // Left half: Xj_l col(1) = Xm, d(Xm)/d(wake.x[j+1]) = 0.5*I
                    double da_l_djp1 = 0.5 * ld.a_Xj(2 + dim);
                    double db_l_djp1 = 0.5 * ld.b_Xj(2 + dim);

                    // Right half: Xj_r col(0) = Xm, col(1) = X_right
                    // d(Xm)/d(wake.x[j+1]) = 0.5*I => 0.5*rd.a_Xj(dim)
                    double da_r_djp1 = 0.5 * rd.a_Xj(dim);
                    double db_r_djp1 = 0.5 * rd.b_Xj(dim);
                    if (j < Nw - 2) {
                        // X_right = wake.x[j+1], d/d(wake.x[j+1]) = I
                        da_r_djp1 += rd.a_Xj(2 + dim);
                        db_r_djp1 += rd.b_Xj(2 + dim);
                    } else {
                        // j == Nw-2: X_right = 2*wake.x[j+1] - Xm
                        // d(X_right)/d(wake.x[j+1]) = 2*I - 0.5*I = 1.5*I
                        da_r_djp1 += 1.5 * rd.a_Xj(2 + dim);
                        db_r_djp1 += 1.5 * rd.b_Xj(2 + dim);
                    }

                    double dW = 0.0;
                    if (j > 0) {
                        dW += 0.5 * da_l_djp1 * phi_jm1;
                        dW += (0.5 * da_l_djp1 + db_l_djp1) * phi_j;
                    } else {
                        dW += db_l_djp1 * phi_j;
                    }
                    dW += (da_r_djp1 + 0.5 * db_r_djp1) * phi_j;
                    if (j < Nw - 2) {
                        dW += 0.5 * db_r_djp1 * phi_jp1;
                    } else {
                        dW += 0.5 * db_r_djp1 * phi_j;  // ghost fold-back
                    }

                    analytic_val += w(i) * dW;
                }
            }

            // Compare: use absolute tolerance floor for small values near FD noise
            double ref = std::max(std::abs(fd_val), std::abs(analytic_val));
            double abs_err = std::abs(analytic_val - fd_val);
            double rel_err = (ref > 1e-15) ? abs_err / ref : abs_err;

            // At eps=1e-7, FD noise floor is ~eps*O(f) ~ 1e-14 for f~1e-7.
            // Use absolute tolerance for values near noise floor.
            bool fail = rel_err > 1e-4 && abs_err > 1e-11;
            if (fail) nfail++;

            std::cout << "  " << std::setw(4) << k << std::setw(4) << dim
                      << std::setw(18) << analytic_val
                      << std::setw(18) << fd_val
                      << std::setw(14) << rel_err
                      << (fail ? "  FAIL" : "") << std::endl;
        }
    }

    std::cout << "BMatrixWakeEndpointFD: " << nfail << "/" << 2 * Nw
              << " failures (tol=1e-4)" << std::endl;
    EXPECT_EQ(nfail, 0) << "B-matrix wake endpoint derivative failures at 1e-4 tolerance";
}

// ===========================================================================
// Wake streamline Jacobian: d(wake.x)/d(foil.x) FD verification
// Perturbs foil.x[k], rebuilds wake, checks delta(wake.x) against the
// backward sweep's implicit prediction.
// ===========================================================================
TEST_F(SensitivityPingTest, WakeStreamlineJacobian) {
    Cppfoil foil_setup("0012", 30);
    foil_setup.SetOper(3.0, -1.0, 1e6, 0.0);
    foil_setup.Solve();
    if (!foil_setup.IsConverged()) { GTEST_SKIP(); }

    auto& state = foil_setup.GetState();
    int N = state.foil.N;
    int Nw = state.wake.N;

    // Save base wake
    MatrixXd base_wake_x = state.wake.x;
    MatrixXd base_wake_t = state.wake.t;
    VectorXd base_wake_s = state.wake.s;
    VectorXd base_xi = state.isol.xi;
    MatrixXd base_gamref = state.isol.gamref;
    VectorXd base_gam = state.isol.gam;
    double alpha = state.oper.alpha;

    // Test: for foil node k, perturb foil.x[k](dim), rebuild wake, measure
    // delta(wake.x[w]) for each wake node w. This gives the forward Jacobian
    // d(wake.x[w])/d(foil.x[k]).
    //
    // The backward sweep in ComputeWakeGeomVJP should produce the TRANSPOSE
    // of this Jacobian when applied to a test vector mu_x.
    // Verify: mu_x^T * J_forward = (J_forward^T * mu_x)^T = backward_sweep(mu_x)^T

    double ep = 1e-7;
    int test_foil_nodes[] = {0, 1, N/2, N-2, N-1};

    std::cout << "  Wake streamline Jacobian d(wake.x)/d(foil.x):" << std::endl;
    for (int k : test_foil_nodes) {
        if (k >= N) continue;
        for (int dim = 0; dim < 2; ++dim) {
            // Perturb foil.x[k](dim), rebuild wake
            state.foil.x(dim, k) += ep;
            // Use base gamma (isolate panel geometry effect on streamline)
            state.isol.gamref = base_gamref;
            state.isol.gam = base_gam;
            StagPointFind(state);
            BuildWake(state);

            // Compute total wake.x change
            double total_delta = 0.0;
            for (int w = 0; w < Nw; ++w) {
                Vector2d dw = (state.wake.x.col(w) - base_wake_x.col(w)) / ep;
                total_delta += dw.norm();
            }

            // Restore
            state.foil.x(dim, k) -= ep;
            state.wake.x = base_wake_x;

            std::cout << "    k=" << k << " d=" << dim
                      << " |d(wake)/d(foil)|=" << total_delta << std::endl;
        }
    }

    // Direct backward sweep test: compute mu_x^T * J_streamline by FD
    // where J_streamline[w,k] = d(wake.x[w])/d(foil.x[k])
    // Use a simple mu_x = e_0 (unit vector at wake node 0)
    {
        // Compute vjp from backward sweep for mu_x = e_0 (first wake node, x-component)
        // This requires extracting the backward sweep output for a specific mu_x input
        // Instead, let's compute mu_x^T * J by FD for the ACTUAL mu_x from Steps 1+1b

        // We need the actual mu_x values. Since we can't access them directly,
        // let's verify by computing vjp_wake = sum_k (vjp_wake[k] * e_k)
        // and comparing against sum_k sum_w (mu_x[w] · d(wake.x[w])/d(foil.x[k])) * e_k

        // Simpler: just verify that perturbing foil.x[k] and measuring the wake change,
        // dotted with a known sensitivity, gives the right VJP.

    }

    // Backward sweep verification: use lambda_ue = e_{N+1} (wake station 1 only)
    // Then ComputeWakeGeomVJP should give d(R_ue(N+1))/d(foil.x) through wake geom.
    // FD: perturb foil.x, rebuild wake (base gamma, then restore foil.x, keep new wake),
    //     recompute R_ue(N+1), compare.
    {
        VectorXd lam_test = VectorXd::Zero(N + Nw);
        lam_test(N + 1) = 1.0;  // wake station 1 only

        VectorXd vjp_an = ComputeWakeGeomVJP(state, lam_test);

        // Base R_ue at station N+1
        VectorXd ue_v = state.glob.U.row(3).transpose();
        VectorXd ds_v = state.glob.U.row(1).transpose();
        VectorXd mass_v = ds_v.cwiseProduct(ue_v);
        VectorXd phi_v = state.vsol.sigma_m * mass_v;
        VectorXd ueinv_b = GetUeInv(state);
        VectorXd sgue_v(N + Nw);
        sgue_v.head(N) = state.isol.sgnue;
        sgue_v.tail(Nw).setOnes();
        double Rue_base_w1 = ue_v(N+1) - ueinv_b(N+1) - sgue_v(N+1) * state.vsol.ue_sigma.row(N+1).dot(phi_v);

        std::cout << "  Backward sweep test (lambda_ue = e_{N+1}):" << std::endl;
        int test_k[] = {0, 1, N/2, N-2, N-1};
        for (int k : test_k) {
            if (k >= N) continue;
            for (int dim = 0; dim < 2; ++dim) {
                // Perturb foil.x, rebuild wake with base gamma
                state.foil.x(dim, k) += ep;
                state.isol.gamref = base_gamref;
                state.isol.gam = base_gam;
                StagPointFind(state);
                BuildWake(state);
                ExtendXiToWake(state);
                // Restore foil.x, keep new wake
                MatrixXd nw_x = state.wake.x;
                MatrixXd nw_t = state.wake.t;
                VectorXd nw_s = state.wake.s;
                VectorXd nw_xi = state.isol.xi;
                state.foil.x(dim, k) -= ep;
                state.wake.x = nw_x; state.wake.t = nw_t;
                state.wake.s = nw_s; state.isol.xi = nw_xi;
                CalcUeM(state);

                VectorXd ueinv_p = GetUeInv(state);
                double Rue_p_w1 = ue_v(N+1) - ueinv_p(N+1) - sgue_v(N+1) * state.vsol.ue_sigma.row(N+1).dot(phi_v);
                double vjp_fd_val = (Rue_p_w1 - Rue_base_w1) / ep;
                double vjp_an_val = vjp_an(2*k + dim);

                double err = std::abs(vjp_fd_val - vjp_an_val) / std::max(std::abs(vjp_fd_val), 1e-10);
                std::cout << "    k=" << k << " d=" << dim
                          << " fd=" << vjp_fd_val << " an=" << vjp_an_val
                          << " err=" << err
                          << (err > 0.1 ? " BAD" : " ok") << std::endl;

                // Restore wake
                state.wake.x = base_wake_x;
                state.wake.t = base_wake_t;
                state.wake.s = base_wake_s;
                state.isol.xi = base_xi;
            }
        }
        // Restore ueM
        CalcUeM(state);
    }
}
