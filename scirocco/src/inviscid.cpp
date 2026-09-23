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

    This work has been cleared for public release, distribution unlimited, case
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#include "scirocco/core/inviscid.h"

#include <cmath>
#include <stdexcept>

#include "scirocco/core/geometry.h"
#include "scirocco/core/panel_influence.h"
#include "scirocco/util/math_utils.h"

namespace scirocco {

void InitThermo(FoilState& state) {
    /* Initialize thermodynamic parameters based on operating conditions
       Must match mfoil.py init_thermo exactly */
    double Ma = state.oper.Ma;
    double gam = state.param.gam;
    double gmi = gam - 1.0;
    double rhoinf = state.oper.rho;
    double Vinf = state.oper.Vinf;

    state.param.Minf = Ma;
    state.param.Vinf = Vinf;

    /* Freestream dynamic viscosity */
    if (state.oper.Re > 0) {
        state.param.muinf = rhoinf * Vinf * state.geom.chord / state.oper.Re;
    }

    double finf = 1.0; /* Sutherland's ratio (1 for incompressible) */

    /* Compressibility corrections (only when Ma > 0) */
    if (Ma > 0) {
        state.param.KTb = std::sqrt(1.0 - Ma * Ma); /* Karman-Tsien beta */
        state.param.KTl =
            Ma * Ma /
            ((1.0 + state.param.KTb) * (1.0 + state.param.KTb)); /* KT lambda */

        /* Stagnation enthalpy: H0 = (1 + 0.5*gmi*Ma²) * Vinf² / (gmi * Ma²) */
        state.param.H0 =
            (1.0 + 0.5 * gmi * Ma * Ma) * Vinf * Vinf / (gmi * Ma * Ma);

        /* Temperature ratio and Sutherland's law */
        double Tr = 1.0 - 0.5 * Vinf * Vinf / state.param.H0; /* T_inf / T_0 */
        finf = std::pow(Tr, 1.5) * (1.0 + state.param.Tsrat) /
               (Tr + state.param.Tsrat);

        /* Sonic Cp: cps = 2/(g*Ma²) * (((1+0.5*gmi*Ma²)/(1+0.5*gmi))^(g/gmi) -
           1) */
        double ratio = (1.0 + 0.5 * gmi * Ma * Ma) / (1.0 + 0.5 * gmi);
        state.param.cps =
            2.0 / (gam * Ma * Ma) * (std::pow(ratio, gam / gmi) - 1.0);
    } else {
        /* Incompressible case */
        state.param.KTb = 1.0;
        state.param.KTl = 0.0;
        /* H0 stays at default (0) for incompressible flow */
        /* cps stays at default for incompressible flow */
    }

    /* Stagnation viscosity (using Sutherland's law) */
    state.param.mu0 = state.param.muinf / finf;

    /* Stagnation density: rho0 = rhoinf * (1 + 0.5*gmi*Ma²)^(1/gmi) */
    state.param.rho0 = rhoinf * std::pow(1.0 + 0.5 * gmi * Ma * Ma, 1.0 / gmi);
}

void ClearSolution(FoilState& state) {
    state.isol = Isol();
    state.vsol = Vsol();
    state.glob = Glob();
    state.post = Post();
    state.wake.N = 0;
    state.wake.x.resize(2, 0);
    state.wake.s.resize(0);
    state.wake.t.resize(2, 0);
}

void BuildGamma(FoilState& state, double alpha) {
    int N = state.foil.N;
    if (N <= 0) {
        throw std::runtime_error("No panels defined");
    }

    /* Get trailing edge info */
    TEInfo te = GetTEInfo(state.foil.x);
    bool nogap = (std::abs(te.hTE) < 1e-10 * state.geom.chord);

    /* Build influence matrix A and right-hand sides for alpha=0,90 */
    MatrixXd A = MatrixXd::Zero(N + 1, N + 1);
    MatrixXd rhs = MatrixXd::Zero(N + 1, 2);

    for (int i = 0; i < N; ++i) {
        Vector2d xi = state.foil.x.col(i);

        /* Linear vortex panel contributions */
        for (int j = 0; j < N - 1; ++j) {
            Matrix2d Xj;
            Xj.col(0) = state.foil.x.col(j);
            Xj.col(1) = state.foil.x.col(j + 1);

            double aij, bij;
            PanelLinVortexStream(Xj, xi, aij, bij);
            A(i, j) += aij;
            A(i, j + 1) += bij;
        }

        /* Streamfunction unknown */
        A(i, N) = -1.0;

        /* Right-hand sides: -Vinf * (y*cos(alpha) - x*sin(alpha)) */
        rhs(i, 0) = -xi(1); /* alpha = 0: -y */
        rhs(i, 1) = xi(0);  /* alpha = 90: x */

        /* TE source panel influence */
        Matrix2d XjTE;
        XjTE.col(0) = state.foil.x.col(N - 1);
        XjTE.col(1) = state.foil.x.col(0);
        double a_src = PanelConstSourceStream(XjTE, xi);
        A(i, 0) += -a_src * (0.5 * te.tcp);
        A(i, N - 1) += a_src * (0.5 * te.tcp);

        /* TE vortex panel influence */
        double a_vor, b_vor;
        PanelLinVortexStream(XjTE, xi, a_vor, b_vor);
        A(i, 0) += -(a_vor + b_vor) * (-0.5 * te.tdp);
        A(i, N - 1) += (a_vor + b_vor) * (-0.5 * te.tdp);
    }

    /* Special N-th equation for no-gap case (extrapolation) */
    if (nogap) {
        A.row(N - 1).setZero();
        A(N - 1, 0) = 1.0;
        A(N - 1, 1) = -2.0;
        A(N - 1, 2) = 1.0;
        A(N - 1, N - 3) = -1.0;
        A(N - 1, N - 2) = 2.0;
        A(N - 1, N - 1) = -1.0;
    }

    /* Kutta condition: gamma[0] + gamma[N-1] = 0 */
    A.row(N).setZero();
    A(N, 0) = 1.0;
    A(N, N - 1) = 1.0;

    /* Solve the system */
    state.isol.AIC = A;
    MatrixXd g = A.partialPivLu().solve(rhs);

    /* Store reference solutions (excluding streamfunction row) */
    state.isol.gamref = g.topRows(N);

    /* Compute gamma for current alpha */
    state.isol.gam = state.isol.gamref.col(0) * Cosd(alpha) +
                     state.isol.gamref.col(1) * Sind(alpha);
}

Vector2d InviscidVelocity(const MatrixXd& X, const VectorXd& G, double Vinf,
                          double alpha, const Vector2d& xi, bool compute_lin,
                          MatrixXd* V_G) {
    int N = static_cast<int>(X.cols());
    Vector2d V = Vector2d::Zero();

    if (compute_lin && V_G != nullptr) {
        V_G->resize(2, N);
        V_G->setZero();
    }

    /* Get TE info for source contribution */
    TEInfo te = GetTEInfo(X);

    /* Linear vortex panel contributions */
    for (int j = 0; j < N - 1; ++j) {
        Matrix2d Xj;
        Xj.col(0) = X.col(j);
        Xj.col(1) = X.col(j + 1);

        Vector2d a, b;
        PanelLinVortexVelocity(Xj, xi, nullptr, false, a, b);
        V += a * G(j) + b * G(j + 1);

        if (compute_lin && V_G != nullptr) {
            V_G->col(j) += a;
            V_G->col(j + 1) += b;
        }
    }

    /* TE source panel contribution */
    Matrix2d XjTE;
    XjTE.col(0) = X.col(N - 1);
    XjTE.col(1) = X.col(0);
    Vector2d v_src = PanelConstSourceVelocity(XjTE, xi, nullptr);
    double sigma_te = 0.5 * te.tcp * (G(N - 1) - G(0));
    V += v_src * sigma_te;

    if (compute_lin && V_G != nullptr) {
        V_G->col(0) += v_src * (-0.5 * te.tcp);
        V_G->col(N - 1) += v_src * (0.5 * te.tcp);
    }

    /* TE vortex panel contribution
       Each node has its own factor: node 0 gets +0.5*tdp, node N-1 gets
       -0.5*tdp */
    Vector2d a_vor, b_vor;
    PanelLinVortexVelocity(XjTE, xi, nullptr, false, a_vor, b_vor);
    Vector2d f_vor = a_vor + b_vor;
    V += f_vor * (0.5 * te.tdp) * G(0) + f_vor * (-0.5 * te.tdp) * G(N - 1);

    if (compute_lin && V_G != nullptr) {
        V_G->col(0) += f_vor * (0.5 * te.tdp);
        V_G->col(N - 1) += f_vor * (-0.5 * te.tdp);
    }

    /* Freestream contribution */
    V(0) += Vinf * Cosd(alpha);
    V(1) += Vinf * Sind(alpha);

    return V;
}

VectorXd GetUeInv(const FoilState& state) {
    if (state.isol.gam.size() == 0) {
        throw std::runtime_error("No inviscid solution");
    }

    double alpha = state.oper.alpha;
    Vector2d cs(Cosd(alpha), Sind(alpha));

    /* Airfoil edge velocity from gamma (tangential velocity = gamma) */
    VectorXd uea = state.isol.sgnue.array() * (state.isol.gamref * cs).array();

    /* Wake edge velocity */
    VectorXd uew;
    if (state.oper.viscous && state.wake.N > 0) {
        uew = state.isol.uewiref * cs;
        uew(0) = uea(uea.size() - 1); /* Continuity */
    }

    /* Concatenate */
    VectorXd ueinv(uea.size() + uew.size());
    ueinv << uea, uew;

    return ueinv;
}

MatrixXd GetUeInvRef(const FoilState& state) {
    if (state.isol.gam.size() == 0) {
        throw std::runtime_error("No inviscid solution");
    }

    /* Airfoil reference velocities */
    MatrixXd uearef = state.isol.gamref;
    for (int i = 0; i < uearef.rows(); ++i) {
        uearef.row(i) *= state.isol.sgnue(i);
    }

    /* Wake reference velocities */
    MatrixXd uewref;
    if (state.oper.viscous && state.wake.N > 0) {
        uewref = state.isol.uewiref;
        uewref.row(0) = uearef.row(uearef.rows() - 1); /* Continuity */
    }

    /* Concatenate */
    MatrixXd ueinvref(uearef.rows() + uewref.rows(), 2);
    ueinvref << uearef, uewref;

    return ueinvref;
}

void SolveInviscid(FoilState& state) {
    if (state.foil.N <= 0) {
        throw std::runtime_error("No panels defined");
    }

    state.oper.viscous = false;
    InitThermo(state);

    /* Set sign of ue (all positive for inviscid) */
    state.isol.sgnue = VectorXd::Ones(state.foil.N);

    /* Build and solve for gamma */
    BuildGamma(state, state.oper.alpha);

    state.glob.conv = true;
}

}  // namespace scirocco
