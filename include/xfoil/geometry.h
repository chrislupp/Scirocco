// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Airfoil geometry utilities.
// Ported from XFOIL xgeom.f and naca.f (Mark Drela, 2000).

#ifndef XFOIL_GEOMETRY_H_
#define XFOIL_GEOMETRY_H_

#include <string>

#include "xfoil/config.h"

namespace xfoil {

/// @brief Result of NACA airfoil generation.
struct NacaResult {
    VectorXd x;      ///< X coordinates (TE -> upper -> LE -> lower -> TE)
    VectorXd y;      ///< Y coordinates
    int n;           ///< Number of points
    std::string name;  ///< Airfoil name (e.g., "NACA 2412")
};

/// @brief Geometric parameters of an airfoil.
struct GeomParams {
    double sle;      ///< Arc length at leading edge
    double xle;      ///< X coordinate of leading edge
    double yle;      ///< Y coordinate of leading edge
    double xte;      ///< X coordinate of trailing edge
    double yte;      ///< Y coordinate of trailing edge
    double chord;    ///< Chord length
    double radle;    ///< Leading edge radius
    double angte;    ///< Trailing edge angle [rad]
};

/// @brief Trailing edge geometry parameters.
struct TEParams {
    double dxte;     ///< TE gap x component (x(1) - x(n))
    double dyte;     ///< TE gap y component (y(1) - y(n))
    double dste;     ///< Total TE gap (magnitude)
    double ante;     ///< Normal-projected TE gap area
    double aste;     ///< Streamwise-projected TE gap area
    double scs;      ///< cos-like TE direction factor
    double sds;      ///< sin-like TE direction factor
    bool sharp;      ///< True if TE is effectively sharp
};

/// @brief Panel generation parameters.
struct PanelParams {
    int npan = 160;          ///< Number of panel nodes
    double cvpar = 1.0;      ///< Curvature-attraction parameter
    double cterat = 0.15;    ///< TE/LE panel density ratio
    double ctrrat = 0.2;     ///< Refined-area/LE panel density ratio
    double xsref1 = 1.0;     ///< Top side refinement start x/c
    double xsref2 = 1.0;     ///< Top side refinement end x/c
    double xpref1 = 1.0;     ///< Bottom side refinement start x/c
    double xpref2 = 1.0;     ///< Bottom side refinement end x/c
};

/// @brief Result of panel generation.
struct PanelResult {
    VectorXd x;      ///< X coordinates of paneled airfoil
    VectorXd y;      ///< Y coordinates of paneled airfoil
    VectorXd s;      ///< Arc length parameter
    int n;           ///< Number of panel nodes
    double sle;      ///< Arc length at leading edge
    bool converged;  ///< True if Newton iteration converged
};

/// @brief Generate 4-digit NACA airfoil from designation.
/// @param ides 4-digit NACA designation (e.g., 2412)
/// @param nside Number of points per side
NacaResult naca4(int ides, int nside);

/// @brief Generate 4-digit NACA airfoil from parameters.
/// @param m Maximum camber (fraction of chord)
/// @param p Position of max camber (fraction of chord)
/// @param t Maximum thickness (fraction of chord)
/// @param nside Number of points per side
NacaResult naca4_params(double m, double p, double t, int nside);

/// @brief Generate 5-digit NACA airfoil from designation.
/// @param ides 5-digit NACA designation (e.g., 23012)
/// @param nside Number of points per side
NacaResult naca5(int ides, int nside);

/// @brief Locate leading edge spline-parameter value.
/// @param x X coordinates
/// @param xp X spline derivatives
/// @param y Y coordinates
/// @param yp Y spline derivatives
/// @param s Arc length array
/// @return Arc length parameter at leading edge
double lefind(const VectorXd& x, const VectorXd& xp,
              const VectorXd& y, const VectorXd& yp,
              const VectorXd& s);

/// @brief Calculate geometric parameters of an airfoil.
GeomParams geopar(const VectorXd& x, const VectorXd& xp,
                  const VectorXd& y, const VectorXd& yp,
                  const VectorXd& s);

/// @brief Normalize airfoil to unit chord with LE at origin.
void normalize(VectorXd& x, VectorXd& y, VectorXd& s);

/// @brief Calculate normal unit vectors at airfoil panel nodes.
void ncalc(const VectorXd& x, const VectorXd& y, const VectorXd& s,
           VectorXd& nx, VectorXd& ny);

/// @brief Calculate trailing edge geometry parameters.
TEParams tecalc(const VectorXd& x, const VectorXd& xp,
                const VectorXd& y, const VectorXd& yp, double chord);

/// @brief Generate panel distribution from airfoil geometry.
PanelResult pangen(const VectorXd& xb, const VectorXd& yb,
                   const PanelParams& params);

}  // namespace xfoil

#endif  // XFOIL_GEOMETRY_H_
