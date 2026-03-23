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
// 8. xi derivatives for ALL influence functions
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
// 15. Inviscid Reynolds Sensitivity: should be zero
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
