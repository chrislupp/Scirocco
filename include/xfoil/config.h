// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef XFOIL_CONFIG_H_
#define XFOIL_CONFIG_H_

#include <Eigen/Dense>

namespace xfoil {

/// @brief Array type aliases using Eigen for linear algebra.
/// Dynamic sizes for arbitrary number of points (no hardcoded limits).

template <typename Scalar, int N>
using Vector = Eigen::Matrix<Scalar, N, 1>;

template <typename Scalar, int Rows, int Cols>
using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

template <typename Scalar>
using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

template <typename Scalar>
using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/// @brief Operating conditions (design variables for optimization).
template <typename Scalar>
struct FlightCondition {
    Scalar alpha;       ///< Angle of attack [rad]
    Scalar re;          ///< Reynolds number
    Scalar mach;        ///< Freestream Mach number
    Scalar xstrip[2];   ///< Forced transition locations (upper/lower)
};

/// @brief Airfoil geometry state (design variables for optimization).
template <typename Scalar>
struct Geometry {
    int n;                          ///< Number of points
    VectorX<Scalar> x, y, s;        ///< Coordinates, arc length
    VectorX<Scalar> xp, yp;         ///< Spline derivatives dx/ds, dy/ds
    VectorX<Scalar> nx, ny;         ///< Normal vectors
    Scalar xle, yle, xte, yte;      ///< LE/TE coordinates
    Scalar chord, sle;              ///< Chord, arc length to LE
    bool sharp_te;                  ///< Sharp trailing edge flag
};

/// @brief Panel solution state.
template <typename Scalar>
struct PanelSolution {
    VectorX<Scalar> gam;    ///< Vorticity strength
    VectorX<Scalar> sig;    ///< Source strength (viscous)
    VectorX<Scalar> qinv;   ///< Inviscid velocity
    VectorX<Scalar> qvis;   ///< Viscous velocity
    VectorX<Scalar> cpi;    ///< Inviscid pressure coefficient
    VectorX<Scalar> cpv;    ///< Viscous pressure coefficient
};

/// @brief Influence coefficient matrices.
template <typename Scalar>
struct InfluenceMatrices {
    MatrixX<Scalar> aij;    ///< Vortex influence (NxN)
    MatrixX<Scalar> dij;    ///< Source influence (NxN)
    bool valid;             ///< True if matrices are up to date
};

/// @brief Boundary layer station state.
template <typename Scalar>
struct BLStation {
    Scalar x;       ///< Arc length position
    Scalar theta;   ///< Momentum thickness
    Scalar dstar;   ///< Displacement thickness
    Scalar uedg;    ///< Edge velocity
    Scalar ctau;    ///< Shear coefficient / amplification
    Scalar cf;      ///< Skin friction coefficient
    Scalar hk;      ///< Shape parameter H = delta*/theta
};

/// @brief BL Newton system (3x5 blocks per station).
template <typename Scalar>
struct BLNewtonBlock {
    Matrix<Scalar, 3, 5> vs1;     ///< d(residual)/d(vars at station 1)
    Matrix<Scalar, 3, 5> vs2;     ///< d(residual)/d(vars at station 2)
    Vector<Scalar, 3> vsrez;      ///< Residual vector
    Vector<Scalar, 3> vsm;        ///< d(residual)/d(Mach)
    Vector<Scalar, 3> vsr;        ///< d(residual)/d(Re)
};

/// @brief Output aerodynamic coefficients.
template <typename Scalar>
struct AeroCoefficients {
    Scalar cl;      ///< Lift coefficient
    Scalar cd;      ///< Total drag coefficient
    Scalar cm;      ///< Moment coefficient
    Scalar cdp;     ///< Pressure drag
    Scalar cdf;     ///< Friction drag
};

}  // namespace xfoil

#endif  // XFOIL_CONFIG_H_
