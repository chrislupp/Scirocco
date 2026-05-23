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
#include "cppfoil/core/boundary_layer.h"
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/forces.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/wake.h"
#include "cppfoil/cppfoil.h"
#include "cppfoil/data/param.h"
#include "cppfoil/data/state.h"
#include "cppfoil/properties/bl_properties.h"
#include "cppfoil/util/linearization.h"

using namespace cppfoil;

// Minimum convergence rate (mfoil uses 1.5 as alarm threshold)
static constexpr double MIN_RATE = 1.5;
// Machine-precision threshold below which we skip rate check
// Raised from mfoil's 5e-11 to account for linear functions (e.g. GetRet w.r.t.
// th) where FD errors are O(eps_mach * |f|/ep) ~ 1e-9 and don't converge
static constexpr double EPS_SKIP = 5e-9;

// ---------------------------------------------------------------------------
// check_ping: 3-evaluation convergence check matching mfoil check_ping
// v[0..2] at 0, +ep, +2*ep;  rate = log2(E[1]/E[0])
// ---------------------------------------------------------------------------
static double check_ping(double ep, const Eigen::VectorXd& v0,
                         const Eigen::VectorXd& v_u0, const Eigen::VectorXd& v1,
                         const Eigen::VectorXd& v_u1, const Eigen::VectorXd& v2,
                         const Eigen::VectorXd& v_u2) {
    double E0 = ((v1 - v0) / ep - 0.5 * (v_u0 + v_u1)).norm();
    double E1 = ((v2 - v0) / (2.0 * ep) - 0.5 * (v_u0 + v_u2)).norm();
    if (E0 < EPS_SKIP && E1 < EPS_SKIP) return 99.0;
    double rate = std::log2(E1 / E0);
    return rate;
}

// ---------------------------------------------------------------------------
// Helper: make a Param for standalone BL property tests
// ---------------------------------------------------------------------------
static Param MakeParam(double Re, double Ma, bool turb, bool wake) {
    Cppfoil foil("0012", 130);
    foil.SetOper(3.0, -1.0, Re, Ma);
    auto& state = foil.GetState();
    InitThermo(state);
    Param p = BuildParam(state, 1);
    p.turb = turb;
    p.wake = wake;
    return p;
}

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------
class PingTest : public ::testing::Test {
   protected:
    static std::unique_ptr<Cppfoil> foil_;

    static void SetUpTestSuite() {
        foil_ = std::make_unique<Cppfoil>("0012", 130);
        foil_->SetOper(3.0, -1.0, 1e6, 0.0);
        auto& s = foil_->GetState();
        s.oper.viscous = true;
        s.oper.initbl = true;
        s.param.verb = 0;
        s.param.niglob = 25;
        s.param.rtol = 1e-6;
        SolveViscous(s);
        CalcForce(s);
    }

    FoilState CopyState() { return foil_->GetState(); }

    static std::mt19937& rng() {
        static std::mt19937 gen(42);
        return gen;
    }

    static VectorXd RandVec(int n) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        VectorXd v(n);
        for (int i = 0; i < n; ++i) v(i) = dist(rng());
        return v;
    }

    static MatrixXd RandMat(int rows, int cols) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        MatrixXd m(rows, cols);
        for (int j = 0; j < cols; ++j)
            for (int i = 0; i < rows; ++i) m(i, j) = dist(rng());
        return m;
    }
};

std::unique_ptr<Cppfoil> PingTest::foil_;

// ===========================================================================
// A. BL Property Function Derivatives
// Matches mfoil exactly: ep and ep/2, rate = log2(E[0]/E[1])
// ===========================================================================

TEST_F(PingTest, BLPropertyDerivatives) {
    double Rev[] = {2e3, 1e5};
    Vector4d Ulam;
    Ulam << 0.01, 0.02, 8.4, 0.9;
    Vector4d Uturb;
    Uturb << 0.023, 0.05, 0.031, 1.1;

    struct RegimeInfo {
        bool turb, wake;
        const char* name;
    };
    RegimeInfo regimes[] = {
        {false, false, "lam"},
        {true, false, "turb"},
        {true, true, "wake"},
    };

    using BLFunc = std::function<ScalarLin4(const Vector4d&, const Param&)>;
    struct FuncEntry {
        const char* name;
        BLFunc func;
    };
    FuncEntry funcs[] = {
        {"GetHk",
         [](const Vector4d& U, const Param& p) { return GetHk(U, p); }},
        {"GetRet",
         [](const Vector4d& U, const Param& p) { return GetRet(U, p); }},
        {"GetCf",
         [](const Vector4d& U, const Param& p) { return GetCf(U, p); }},
        {"GetCDi",
         [](const Vector4d& U, const Param& p) { return GetCDi(U, p); }},
        {"GetHs",
         [](const Vector4d& U, const Param& p) { return GetHs(U, p); }},
        {"GetHss",
         [](const Vector4d& U, const Param& p) { return GetHss(U, p); }},
        {"GetUs",
         [](const Vector4d& U, const Param& p) { return GetUs(U, p); }},
        {"GetCteq",
         [](const Vector4d& U, const Param& p) { return GetCteq(U, p); }},
        {"GetCttr",
         [](const Vector4d& U, const Param& p) { return GetCttr(U, p); }},
        {"GetDe",
         [](const Vector4d& U, const Param& p) { return GetDe(U, p); }},
        {"GetDamp",
         [](const Vector4d& U, const Param& p) { return GetDamp(U, p); }},
        {"GetMach2",
         [](const Vector4d& U, const Param& p) { return GetMach2(U, p); }},
        {"GetRho",
         [](const Vector4d& U, const Param& p) { return GetRho(U, p); }},
    };

    int nfail = 0;
    for (double Re : Rev) {
        for (const auto& regime : regimes) {
            Param param = MakeParam(Re, 0.4, regime.turb, regime.wake);
            Vector4d Ubase = regime.turb ? Uturb : Ulam;

            for (const auto& fe : funcs) {
                std::string srates, serr;
                bool any_bad = false;

                for (int k = 0; k < 4; ++k) {
                    Vector4d U0 = Ubase;
                    double ep = 1e-2 * Ubase(k);
                    double E[2];

                    // Evaluate at base point
                    auto [v0, v_U0] = fe.func(U0, param);

                    // Two epsilons: ep and ep/2 (matching mfoil exactly)
                    double epi = ep;
                    for (int iep = 0; iep < 2; ++iep) {
                        U0(k) = Ubase(k) + epi;
                        auto [v1, v_U1] = fe.func(U0, param);
                        E[iep] = std::abs((v1 - v0) / epi -
                                          0.5 * (v_U1(k) + v_U0(k)));
                        epi /= 2.0;
                    }
                    U0(k) = Ubase(k);

                    // rate = log2(E[0]/E[1]) — larger ep first
                    double rate = 99.0;
                    if (E[0] > EPS_SKIP && E[1] > EPS_SKIP) {
                        rate = std::log2(E[0] / E[1]);
                    }

                    if (rate < MIN_RATE) {
                        any_bad = true;
                        nfail++;
                        ADD_FAILURE()
                            << fe.name << " " << regime.name << " Re=" << Re
                            << " k=" << k << " rate=" << rate << " E=[" << E[0]
                            << "," << E[1] << "]";
                    }
                }
            }
        }
    }
    std::cout << "BLPropertyDerivatives: " << nfail << " failures" << std::endl;
}

// ===========================================================================
// A2. ResidualStation derivatives (non-stagnation)
// Uses ep and ep/2 approach matching mfoil
// ===========================================================================

TEST_F(PingTest, ResidualStationDerivatives) {
    double Rev[] = {2e3, 1e5};
    Vector4d Ulam;
    Ulam << 0.01, 0.02, 8.4, 0.9;
    Vector4d Uturb;
    Uturb << 0.023, 0.05, 0.031, 1.1;

    struct RegimeInfo {
        bool turb, wake;
        const char* name;
    };
    RegimeInfo regimes[] = {
        {false, false, "lam"},
        {true, false, "turb"},
        {true, true, "wake"},
    };

    Vector2d xi_base;
    xi_base << 0.7, 0.8;
    Vector2d Aux;
    Aux << 0.002, 0.0018;
    Vector2d dx;
    dx << -0.2, 0.3;

    int nfail = 0;
    for (double Re : Rev) {
        for (const auto& regime : regimes) {
            Param param = MakeParam(Re, 0.4, regime.turb, regime.wake);
            Vector4d U1base = regime.turb ? Uturb : Ulam;
            Vector4d U2base;
            U2base << U1base(0) * 1.1, U1base(1) * 0.8, U1base(2) * 0.9,
                U1base(3) * 1.2;

            for (int k = 0; k < 8; ++k) {
                double ep_val = (k < 4) ? U1base(k % 4) : U2base(k % 4);
                double ep = 1e-2 * ep_val;
                double E[2];

                // Evaluate at base
                Matrix42d Umat;
                Umat.col(0) = U1base;
                Umat.col(1) = U2base;
                auto [R0, R_U0, R_x0, R_Ret0] =
                    ResidualStation(param, xi_base, Umat, Aux);
                (void)R_Ret0;

                // Two epsilons: ep and ep/2
                double epi = ep;
                for (int iep = 0; iep < 2; ++iep) {
                    Matrix42d Umat1 = Umat;
                    Umat1(k % 4, k / 4) += epi;
                    Vector2d xi1 = xi_base + epi * dx;
                    auto [R1, R_U1, R_x1, R_Ret1] =
                        ResidualStation(param, xi1, Umat1, Aux);
                    (void)R_Ret1;

                    // Error: FD vs averaged analytical derivative
                    Vector3d v_u0 = R_U0.col(k) + R_x0 * dx;
                    Vector3d v_u1 = R_U1.col(k) + R_x1 * dx;
                    E[iep] = ((R1 - R0) / epi - 0.5 * (v_u0 + v_u1)).norm();
                    epi /= 2.0;
                }

                double rate = 99.0;
                if (E[0] > EPS_SKIP && E[1] > EPS_SKIP) {
                    rate = std::log2(E[0] / E[1]);
                }

                if (rate < MIN_RATE) {
                    nfail++;
                    ADD_FAILURE()
                        << "ResidualStation " << regime.name << " Re=" << Re
                        << " k=" << k << " rate=" << rate << " E=[" << E[0]
                        << "," << E[1] << "]";
                }
            }
        }
    }
    std::cout << "ResidualStationDerivatives: " << nfail << " failures"
              << std::endl;
}

// ===========================================================================
// B. Transition Residual Derivatives
// ===========================================================================

TEST_F(PingTest, ResidualTransitionDerivatives) {
    FoilState state = CopyState();
    state.oper.Re = 2e6;
    state.oper.Ma = 0.0;
    InitThermo(state);
    Param param = BuildParam(state, 1);

    Matrix42d U;
    U.col(0) << 0.01, 0.02, 8.95, 0.9;
    U.col(1) << 0.013, 0.023, 0.028, 0.85;
    Vector2d x;
    x << 0.7, 0.8;
    Vector2d Aux = Vector2d::Zero();

    rng().seed(42);
    Eigen::Matrix<double, 4, 2> dU_mat;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 4; ++i)
            dU_mat(i, j) =
                std::uniform_real_distribution<double>(0.0, 1.0)(rng());
    Vector2d dx;
    dx(0) = std::uniform_real_distribution<double>(0.0, 1.0)(rng());
    dx(1) = std::uniform_real_distribution<double>(0.0, 1.0)(rng());

    double ep = 1e-4;

    std::array<Vector3d, 3> v;
    std::array<Vector3d, 3> v_u;

    for (int ie = 0; ie < 3; ++ie) {
        auto [R, R_U, R_x, R_Ret_unused] =
            ResidualTransition(state, param, x, U, Aux);
        (void)R_Ret_unused;
        v[ie] = R;
        Eigen::Matrix<double, 8, 1> dU_flat;
        dU_flat << dU_mat.col(0), dU_mat.col(1);
        v_u[ie] = R_U * dU_flat + R_x * dx;
        U += ep * dU_mat;
        x += ep * dx;
    }

    double rate = check_ping(ep, v[0], v_u[0], v[1], v_u[1], v[2], v_u[2]);
    EXPECT_GT(rate, MIN_RATE) << "Transition residual rate=" << rate;
    std::cout << "ResidualTransition: rate=" << rate << std::endl;
}

// ===========================================================================
// C. Stagnation Residual Derivatives
// ===========================================================================

TEST_F(PingTest, ResidualStation_Stagnation) {
    Cppfoil foil_stag("0012", 130);
    foil_stag.SetOper(1.0, -1.0, 1e6, 0.0);
    auto& st = foil_stag.GetState();
    InitThermo(st);
    Param param = BuildParam(st, 1);

    Vector4d U;
    U << 0.00004673616, 0.000104289, 0.0, 0.11977917547;
    double xval = 4.590816441485401e-05;
    Vector2d Aux = Vector2d::Zero();

    rng().seed(42);
    Vector4d dU;
    for (int i = 0; i < 4; ++i)
        dU(i) = std::uniform_real_distribution<double>(0.0, 1.0)(rng());
    double dx_scalar = std::uniform_real_distribution<double>(0.0, 1.0)(rng());

    double ep = 1e-6;

    std::array<Vector3d, 3> v;
    std::array<Vector3d, 3> v_u;

    for (int ie = 0; ie < 3; ++ie) {
        param.simi = true;
        Vector2d xvec;
        xvec << xval, xval;
        Matrix42d Umat;
        Umat.col(0) = U;
        Umat.col(1) = U;
        auto [R, R_U, R_x, R_Ret_unused] =
            ResidualStation(param, xvec, Umat, Aux);
        (void)R_Ret_unused;
        param.simi = false;

        v[ie] = R;
        v_u[ie] = (R_U.leftCols(4) + R_U.rightCols(4)) * dU +
                  (R_x.col(0) + R_x.col(1)) * dx_scalar;

        U += ep * dU;
        xval += ep * dx_scalar;
    }

    double rate = check_ping(ep, v[0], v_u[0], v[1], v_u[1], v[2], v_u[2]);
    EXPECT_GT(rate, MIN_RATE) << "Stagnation residual rate=" << rate;
    std::cout << "ResidualStation_Stagnation: rate=" << rate << std::endl;
}

// ===========================================================================
// D. Global System Jacobians (BuildGlobSys)
// ===========================================================================

TEST_F(PingTest, BuildGlobSys_Jacobians) {
    rng().seed(42);

    // Sub-test 1: explicit R_x test
    {
        FoilState state = CopyState();
        int Nsys = state.glob.Nsys;
        MatrixXd dU = RandMat(4, Nsys);
        VectorXd dx = 0.1 * RandVec(Nsys);
        double ep = 1e-6;

        std::array<VectorXd, 3> v;
        std::array<VectorXd, 3> v_u;

        for (int ie = 0; ie < 3; ++ie) {
            BuildGlobSys(state);
            v[ie] = state.glob.R;

            Eigen::VectorXd dU_flat =
                Eigen::Map<const Eigen::VectorXd>(dU.data(), 4 * Nsys);
            v_u[ie] = state.glob.R_U * dU_flat + state.glob.R_x * dx;

            state.glob.U += ep * dU;
            state.isol.xi += ep * dx;
        }

        double rate = check_ping(ep, v[0], v_u[0], v[1], v_u[1], v[2], v_u[2]);
        EXPECT_GT(rate, MIN_RATE) << "BuildGlobSys explicit rate=" << rate;
        std::cout << "BuildGlobSys (explicit): rate=" << rate << std::endl;
    }

    // Sub-test 2: with StagPointMove (baked-in stag derivatives)
    {
        FoilState state = CopyState();
        int Nsys = state.glob.Nsys;
        MatrixXd dU = RandMat(4, Nsys);
        VectorXd dx = VectorXd::Zero(Nsys);
        double ep = 1e-6;

        StagPointMove(state);

        std::array<VectorXd, 3> v;
        std::array<VectorXd, 3> v_u;

        for (int ie = 0; ie < 3; ++ie) {
            BuildGlobSys(state);
            v[ie] = state.glob.R;

            Eigen::VectorXd dU_flat =
                Eigen::Map<const Eigen::VectorXd>(dU.data(), 4 * Nsys);
            v_u[ie] = state.glob.R_U * dU_flat + state.glob.R_x * dx;

            state.glob.U += ep * dU;
            StagPointMove(state);
        }

        double rate = check_ping(ep, v[0], v_u[0], v[1], v_u[1], v[2], v_u[2]);
        EXPECT_GT(rate, MIN_RATE) << "BuildGlobSys StagPointMove rate=" << rate;
        std::cout << "BuildGlobSys (StagPointMove): rate=" << rate << std::endl;
    }
}

// ===========================================================================
// E. Wake System Derivatives
// ===========================================================================

TEST_F(PingTest, WakeSys_Derivatives) {
    FoilState state = CopyState();
    rng().seed(42);

    int Nsys = state.glob.Nsys;
    MatrixXd dU = RandMat(4, Nsys);
    double ep = 1e-5;

    Param param = BuildParam(state, 1);

    std::array<Vector3d, 3> v;
    std::array<VectorXd, 3> v_u;

    for (int ie = 0; ie < 3; ++ie) {
        auto result = WakeSys(state, param);
        v[ie] = result.R;

        Eigen::Matrix<double, 12, 1> dU_J;
        for (int j = 0; j < 3; ++j) {
            dU_J.segment(4 * j, 4) = dU.col(result.J[j]);
        }
        v_u[ie] = result.R_U * dU_J;

        state.glob.U += ep * dU;
    }

    double rate = check_ping(ep, v[0], v_u[0], v[1], v_u[1], v[2], v_u[2]);
    EXPECT_GT(rate, MIN_RATE) << "WakeSys rate=" << rate;
    std::cout << "WakeSys: rate=" << rate << std::endl;
}

// ===========================================================================
// F. Stagnation State Derivatives
// ===========================================================================

TEST_F(PingTest, StagnationState_Derivatives) {
    rng().seed(42);

    Matrix42d U;
    U.col(0) << 5e-5, 1.1e-4, 0.0, 0.0348;
    U.col(1) << 4.9e-5, 1.09e-4, 0.0, 0.07397;
    Vector2d x;
    x << 5.18e-4, 1.1e-3;

    Eigen::Matrix<double, 4, 2> dU_mat;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < 4; ++i)
            dU_mat(i, j) =
                std::uniform_real_distribution<double>(0.0, 1.0)(rng());
    Vector2d dx;
    dx(0) = std::uniform_real_distribution<double>(0.0, 1.0)(rng());
    dx(1) = std::uniform_real_distribution<double>(0.0, 1.0)(rng());

    double ep = 1e-6;

    std::array<Vector4d, 3> v;
    std::array<Vector4d, 3> v_u;

    for (int ie = 0; ie < 3; ++ie) {
        auto result = StagnationState(U, x);
        v[ie] = result.Ust;

        Eigen::Matrix<double, 8, 1> dU_flat;
        dU_flat << dU_mat.col(0), dU_mat.col(1);
        v_u[ie] = result.Ust_U * dU_flat + result.Ust_x * dx;

        U += ep * dU_mat;
        x += ep * dx;
    }

    Eigen::VectorXd v0(4), v1(4), v2(4), vu0(4), vu1(4), vu2(4);
    v0 = v[0];
    v1 = v[1];
    v2 = v[2];
    vu0 = v_u[0];
    vu1 = v_u[1];
    vu2 = v_u[2];

    double rate = check_ping(ep, v0, vu0, v1, vu1, v2, vu2);
    EXPECT_GT(rate, MIN_RATE) << "StagnationState rate=" << rate;
    std::cout << "StagnationState: rate=" << rate << std::endl;
}

// ===========================================================================
// G. Force Calculation Derivatives
// ===========================================================================

TEST_F(PingTest, CalcForce_Derivatives) {
    rng().seed(42);

    FoilState state = CopyState();
    int Nsys = state.glob.Nsys;
    int N = state.foil.N;

    VectorXd due = RandVec(N);
    MatrixXd dU = MatrixXd::Zero(4, Nsys);
    dU.row(3).head(N) = due.transpose();
    double da = 10.0;
    double ep = 1e-2;

    std::array<double, 3> cl_v;
    std::array<double, 3> cl_vu;
    for (int ie = 0; ie < 3; ++ie) {
        CalcForce(state);
        cl_v[ie] = state.post.cl;
        cl_vu[ie] = state.post.cl_ue.dot(due) + state.post.cl_alpha * da;
        state.glob.U += ep * dU;
        state.oper.alpha += ep * da;
    }

    Eigen::VectorXd v0(1), v1(1), v2(1), vu0(1), vu1(1), vu2(1);
    v0(0) = cl_v[0];
    v1(0) = cl_v[1];
    v2(0) = cl_v[2];
    vu0(0) = cl_vu[0];
    vu1(0) = cl_vu[1];
    vu2(0) = cl_vu[2];
    double rate = check_ping(ep, v0, vu0, v1, vu1, v2, vu2);
    EXPECT_GT(rate, MIN_RATE) << "CalcForce cl rate=" << rate;
    std::cout << "CalcForce (cl): rate=" << rate << std::endl;
}
