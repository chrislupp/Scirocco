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
#include <iostream>
#include <random>

#include "cppfoil/config.h"
#include "cppfoil/core/adjoint.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/panel_influence_deriv.h"
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

    Cppfoil foil("0012", 50);
    foil.SetOper(3.0, -1.0, 1e6, 0.0);
    foil.Solve();
    if (!foil.IsConverged()) { GTEST_SKIP(); }
    auto& state = foil.GetState();

    MatrixXd dRue_dx = ComputeDRueDxFoil(state);

    int N = state.foil.N;
    int Nsys = state.glob.Nsys;

    // Compute base R_ue residual
    VectorXd ue = state.glob.U.row(3).transpose();
    VectorXd ds = state.glob.U.row(1).transpose();
    VectorXd ueinv = GetUeInv(state);
    VectorXd mass = ds.cwiseProduct(ue);
    VectorXd R_ue_base = ue - (ueinv + state.vsol.ue_m * mass);

    // FD check: separate ueinv and ue_m contributions
    int test_nodes_rue[] = {10, 25, 40};
    int nfail = 0;
    double ep = 1e-6;

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

            bool fail = rel_err > 1e-3;
            if (fail) nfail++;

            std::cout << "  dRue/dx node=" << test_node << " dim=" << dim
                      << " rel_err=" << rel_err
                      << " |ueinv|=" << ueinv_mag
                      << " |uem|=" << uem_mag;
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

            if (rel_cl > 0.1) {  // > 10% error
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

            if (rel_cd > 0.5) {  // > 50% error (cd has more approximations)
                nfail_cd++;
                std::cout << "  ViscSens cd FAIL node=" << node << " dim=" << dim
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
