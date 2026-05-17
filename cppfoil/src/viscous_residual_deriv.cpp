/*
    cppfoil: C++ Airfoil Analysis Library

    Copyright 2024-2026 Christopher A. Lupp

    Licensed under the Apache License, Version 2.0 (the "License").
*/
#include "cppfoil/core/viscous_residual_deriv.h"

#include <vector>

#include "cppfoil/core/geometry.h"
#include "cppfoil/core/viscous_geom_deriv.h"
#include "cppfoil/properties/bl_properties.h"

namespace cppfoil {

// ============================================================================
// DRblDFoilX
//
//   dR_BL/d(foil.x)
//     = R_x · DXiDFoilX                       (closed-form xi chain)
//     + (wgap chain at non-initial wake stations from ResidualStation)
//     + (hTE chain at the first wake station from WakeSys)
//
// R_x is the [3*Nsys x Nsys] sparse linearisation already assembled by
// BuildGlobSys.
//
// wgap chain: ResidualStation modifies U via "U.row(1) -= Aux", so
//   ∂R_station/∂Aux(0) = -∂R_station/∂U[1] at the first column (lin_U(:, 1)),
//   ∂R_station/∂Aux(1) = -∂R_station/∂U[1] at the second column (lin_U(:, 5))
// These derivatives live in R_U at rows 3*Is[2][ip]..3*Is[2][ip]+2,
// columns 4*Is[2][ip-1]+1 and 4*Is[2][ip]+1 respectively. Aux at wake-station
// local index `i` equals wgap(i), so the chain through DWgapDFoilX is:
//   contribution to dR_BL row r = -R_U(r, 4*Is[2][ip-1]+1) * dwgap(ip-1)
//                                 -R_U(r, 4*Is[2][ip]+1)   * dwgap(ip)
// for each wake-station index ip = 1..Nw-1 (the first wake station is
// handled by WakeSys, not ResidualStation).
//
// hTE chain at WakeSys: residual is R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE),
// so ∂R(1)/∂hTE = -1, and ∂hTE/∂foil.x[k]_dim comes from GetTEInfoDeriv.
// Only the second component of the first-wake-station residual triple is
// affected.
// ============================================================================
MatrixXd DRblDFoilX(const FoilState& state, const MatrixXd& dfoil_s_dx) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;

    if (state.glob.R_x.rows() != 3 * Nsys || state.glob.R_x.cols() != Nsys) {
        throw std::runtime_error(
            "DRblDFoilX: state.glob.R_x not assembled — call BuildGlobSys "
            "first");
    }

    MatrixXd dxi_dx = DXiDFoilX(state, dfoil_s_dx);    // [Nsys x 2N]
    MatrixXd dwgap_dx = DWgapDFoilX(state, dfoil_s_dx); // [Nw   x 2N]

    // ---- R_x · DXi ----
    MatrixXd dRbl = state.glob.R_x * dxi_dx;  // sparse * dense → dense

    // ---- wgap chain at wake stations ip = 1..Nw-1 ----
    //
    // ResidualStation modifies U via "U.row(1) -= Aux" AND uses Aux directly
    // in GetHw (Hw = Aux/U[0]). So Aux has TWO chain paths into R:
    //   1. U_mod[1] = U[1] - Aux → dU_mod[1]/dAux = -1, contributes -R_U[:, U1_col]
    //   2. Hw = Aux/U[0] → dHw/dAux = 1/U[0], contributes ∂R/∂Hw * 1/U[0]
    //
    // Hw appears in (Rmom, Rshape, Rlag) with coefficients (uelog, -uelog, 0).
    // Total ∂R/∂Aux = -R_U[:, U1_col] + [uelog, -uelog, 0] / U[0].
    //
    // For wake stations with cf=0, these two paths CANCEL for Rmom and Rshape
    // (the U_mod chain through H = (U[1]-Aux)/U[0] gives ±uelog/U[0] that
    // exactly cancels the Hw-direct path). Missing the Hw path produces a
    // spurious 0.5-0.9-magnitude chain term — caught at the Phase 2 test.
    const auto& Iw = state.vsol.Is[2];
    if (static_cast<int>(Iw.size()) != Nw) {
        throw std::runtime_error("DRblDFoilX: Is[wake].size() != Nw");
    }
    const SparseMatrixd& RU = state.glob.R_U;
    for (int ip = 1; ip < Nw; ++ip) {
        int Ig = 3 * Iw[ip];
        int node_prev = Iw[ip - 1];
        int node_curr = Iw[ip];
        int col_prev = 4 * node_prev + 1;
        int col_curr = 4 * node_curr + 1;

        // Compute uelog and U[0] at the two endpoints (Karman-Tsien-corrected
        // edge velocity uk via GetUk to match residual_station internals).
        auto [uk_prev, _uk_U_prev] = GetUk(state.glob.U(3, node_prev), state.param);
        auto [uk_curr, _uk_U_curr] = GetUk(state.glob.U(3, node_curr), state.param);
        double uelog = std::log(uk_curr / uk_prev);
        double U0_prev = state.glob.U(0, node_prev);
        double U0_curr = state.glob.U(0, node_curr);

        // Hw in residual is the AVERAGE: Hw_avg = 0.5*(Hw1 + Hw2). So the
        // per-endpoint Hw chain contributes 0.5*∂R/∂Hw_avg, not the full thing.
        // R coefficients for Hw_avg: +uelog (Rmom), -uelog (Rshape), 0 (Rlag).
        // Per-endpoint Hw_i: factor 0.5 from the average.
        Vector3d Hw_coef(0.5 * uelog, -0.5 * uelog, 0.0);

        for (int r = 0; r < 3; ++r) {
            double cprev = RU.coeff(Ig + r, col_prev);
            double ccurr = RU.coeff(Ig + r, col_curr);
            // Correct ∂R/∂Aux = (-R_U via U_mod chain) + (Hw chain via 1/U[0])
            double dRdAux_prev = -cprev + Hw_coef(r) / U0_prev;
            double dRdAux_curr = -ccurr + Hw_coef(r) / U0_curr;
            // Chain to foil.x: contribute (∂R/∂Aux) * dwgap_dx
            if (dRdAux_prev != 0.0) {
                dRbl.row(Ig + r) += dRdAux_prev * dwgap_dx.row(ip - 1);
            }
            if (dRdAux_curr != 0.0) {
                dRbl.row(Ig + r) += dRdAux_curr * dwgap_dx.row(ip);
            }
        }
    }

    // ---- hTE chain at the first wake station (WakeSys) ----
    // R(1) = Uw[1] - (Ul[1] + Uu[1] + hTE), so dR(1)/dhTE = -1.
    // hTE depends on foil.x[k]_dim only for k ∈ {0, 1, N-2, N-1}.
    int Ig_wake_init = 3 * Iw[0] + 1;  // second residual of first wake station
    for (int k : {0, 1, N - 2, N - 1}) {
        for (int dim = 0; dim < 2; ++dim) {
            double dtcp, dtdp, dhTE, ddtdx;
            GetTEInfoDeriv(state.foil.x, k, dim, dtcp, dtdp, dhTE, ddtdx);
            if (dhTE != 0.0) {
                dRbl(Ig_wake_init, 2 * k + dim) += -dhTE;
            }
        }
    }

    return dRbl;
}

// ============================================================================
// DRueDFoilX — stub (to be filled in next step).
// ============================================================================
MatrixXd DRueDFoilX(const FoilState& state) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int Nsys = N + Nw;
    return MatrixXd::Zero(Nsys, 2 * N);
}

}  // namespace cppfoil
