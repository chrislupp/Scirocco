// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Panel method influence calculations.
// Ported from XFOIL xpanel.f (Mark Drela, 2000).

#ifndef XFOIL_PANEL_H_
#define XFOIL_PANEL_H_

#include "xfoil/config.h"

namespace xfoil {

/// @brief Result of PSILIN panel influence calculation.
struct PsiLinResult {
    double psi;     ///< Stream function contribution
    double psi_ni;  ///< dPsi/dn at field point (normal deriv)
    double qtan;    ///< Tangential velocity contribution
};

/// @brief Result of complete panel influence kernel.
struct PanelInfluence {
    double psi;     ///< Total stream function
    double psi_gam; ///< dPsi/dGam (vorticity sensitivity)
    double psi_sig; ///< dPsi/dSig (source sensitivity)
    double qtan;    ///< Tangential velocity
    double qtan_gam;///< dQtan/dGam
    double qtan_sig;///< dQtan/dSig
};

/// @brief Calculate panel angles from panel node coordinates.
/// @param x X coordinates of panel nodes
/// @param y Y coordinates of panel nodes
/// @param apanel Output: panel angles [rad] (size n-1 for n nodes)
/// @note Panel angle is direction from node i to node i+1
void apcalc(const VectorXd& x, const VectorXd& y, VectorXd& apanel);

/// @brief Calculate atan2 with 0..2*pi range (no branch cuts).
/// @param y Y component
/// @param x X component
/// @return Angle in [0, 2*pi)
double atanc(double y, double x);

/// @brief Single panel influence kernel - vortex and source contributions.
///
/// Calculates the stream function contribution from a constant-strength
/// vortex/source panel to a field point. This is the core influence
/// calculation used to build the AIC and DIC matrices.
///
/// @param xf Field point X
/// @param yf Field point Y
/// @param x1 Panel start X
/// @param y1 Panel start Y
/// @param x2 Panel end X
/// @param y2 Panel end Y
/// @param gam Panel vorticity strength (positive CCW)
/// @param sig Panel source strength
/// @return Stream function and velocity contributions
PanelInfluence psilin(double xf, double yf,
                      double x1, double y1,
                      double x2, double y2,
                      double gam, double sig);

/// @brief Build vortex influence coefficient matrix.
///
/// Computes the AIC matrix where Aij is the stream function at node i
/// due to unit vorticity on panel j. Uses the Kutta condition as the
/// last row (gamma_1 + gamma_n = 0 for sharp TE).
///
/// @param x Panel node X coordinates
/// @param y Panel node Y coordinates
/// @param apanel Panel angles from apcalc()
/// @param aij Output: AIC matrix (n x n)
/// @param sharp_te True if trailing edge is sharp (affects Kutta condition)
void aic_build(const VectorXd& x, const VectorXd& y,
               const VectorXd& apanel, MatrixXd& aij, bool sharp_te);

/// @brief Build source influence coefficient matrix.
///
/// Computes the DIC matrix where Dij is the stream function at node i
/// due to unit source strength on panel j.
///
/// @param x Panel node X coordinates
/// @param y Panel node Y coordinates
/// @param apanel Panel angles from apcalc()
/// @param dij Output: DIC matrix (n x n)
void dic_build(const VectorXd& x, const VectorXd& y,
               const VectorXd& apanel, MatrixXd& dij);

/// @brief Calculate tangential velocities from panel solution.
///
/// Computes surface velocity at each panel node from the vorticity
/// and source distributions using influence matrix multiplication.
///
/// @param x Panel node X coordinates
/// @param y Panel node Y coordinates
/// @param apanel Panel angles
/// @param gam Vorticity distribution
/// @param sig Source distribution (zero for inviscid)
/// @param alpha Angle of attack [rad]
/// @param qinv Output: inviscid surface velocity
void qcalc(const VectorXd& x, const VectorXd& y,
           const VectorXd& apanel, const VectorXd& gam,
           const VectorXd& sig, double alpha, VectorXd& qinv);

/// @brief Calculate surface velocity at a specific arc length.
/// @param si Arc length parameter
/// @param x Panel node X coordinates
/// @param xp X spline derivatives
/// @param y Panel node Y coordinates
/// @param yp Y spline derivatives
/// @param s Arc length array
/// @param gam Vorticity distribution
/// @param sig Source distribution
/// @param alpha Angle of attack [rad]
/// @return Surface velocity magnitude
double qinf(double si, const VectorXd& x, const VectorXd& xp,
            const VectorXd& y, const VectorXd& yp, const VectorXd& s,
            const VectorXd& gam, const VectorXd& sig, double alpha);

/// @brief Calculate freestream velocity components.
/// @param alpha Angle of attack [rad]
/// @param qinf Freestream velocity magnitude
/// @param uinf Output: X component
/// @param vinf Output: Y component
void freestream(double alpha, double qinf, double& uinf, double& vinf);

}  // namespace xfoil

#endif  // XFOIL_PANEL_H_
