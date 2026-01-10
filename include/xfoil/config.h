// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef XFOIL_CONFIG_H_
#define XFOIL_CONFIG_H_

#include <Eigen/Dense>

namespace xfoil {

// Vector and matrix type aliases
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using Vector3d = Eigen::Vector3d;
using Matrix35d = Eigen::Matrix<double, 3, 5>;

/// @brief Operating conditions.
struct FlightCondition {
    double alpha;       ///< Angle of attack [rad]
    double re;          ///< Reynolds number
    double mach;        ///< Freestream Mach number
    double xstrip[2];   ///< Forced transition locations (upper/lower)
};

/// @brief Airfoil geometry state.
struct Geometry {
    int n;                  ///< Number of points
    VectorXd x, y, s;       ///< Coordinates, arc length
    VectorXd xp, yp;        ///< Spline derivatives dx/ds, dy/ds
    VectorXd nx, ny;        ///< Normal vectors
    double xle, yle;        ///< LE coordinates
    double xte, yte;        ///< TE coordinates
    double chord, sle;      ///< Chord, arc length to LE
    bool sharp_te;          ///< Sharp trailing edge flag
};

/// @brief Panel solution state.
struct PanelSolution {
    VectorXd gam;    ///< Vorticity strength
    VectorXd sig;    ///< Source strength (viscous)
    VectorXd qinv;   ///< Inviscid velocity
    VectorXd qvis;   ///< Viscous velocity
    VectorXd cpi;    ///< Inviscid pressure coefficient
    VectorXd cpv;    ///< Viscous pressure coefficient
};

/// @brief Influence coefficient matrices.
struct InfluenceMatrices {
    MatrixXd aij;    ///< Vortex influence (NxN)
    MatrixXd dij;    ///< Source influence (NxN)
    bool valid;      ///< True if matrices are up to date
};

/// @brief Boundary layer station state.
struct BLStation {
    double x;       ///< Arc length position
    double theta;   ///< Momentum thickness
    double dstar;   ///< Displacement thickness
    double uedg;    ///< Edge velocity
    double ctau;    ///< Shear coefficient / amplification
    double cf;      ///< Skin friction coefficient
    double hk;      ///< Shape parameter H = delta*/theta
};

/// @brief BL Newton system (3x5 blocks per station).
struct BLNewtonBlock {
    Matrix35d vs1;   ///< d(residual)/d(vars at station 1)
    Matrix35d vs2;   ///< d(residual)/d(vars at station 2)
    Vector3d vsrez;  ///< Residual vector
    Vector3d vsm;    ///< d(residual)/d(Mach)
    Vector3d vsr;    ///< d(residual)/d(Re)
};

/// @brief Output aerodynamic coefficients.
struct AeroCoefficients {
    double cl;      ///< Lift coefficient
    double cd;      ///< Total drag coefficient
    double cm;      ///< Moment coefficient
    double cdp;     ///< Pressure drag
    double cdf;     ///< Friction drag
};

}  // namespace xfoil

#endif  // XFOIL_CONFIG_H_
