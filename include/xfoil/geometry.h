// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Airfoil geometry utilities.
// Ported from XFOIL xgeom.f and naca.f (Mark Drela, 2000).

#ifndef XFOIL_GEOMETRY_H_
#define XFOIL_GEOMETRY_H_

#include <cmath>
#include <string>

#include "xfoil/config.h"
#include "xfoil/spline.h"

namespace xfoil {

/// @brief Result of NACA airfoil generation.
template <typename Scalar>
struct NacaResult {
    VectorX<Scalar> x;   ///< X coordinates (TE -> upper -> LE -> lower -> TE)
    VectorX<Scalar> y;   ///< Y coordinates
    int n;               ///< Number of points
    std::string name;    ///< Airfoil name (e.g., "NACA 2412")
};

/// @brief Generate 4-digit NACA airfoil from designation (e.g., 2412).
/// @param ides 4-digit NACA designation (e.g., 2412)
/// @param nside Number of points per side
/// @return NacaResult with airfoil coordinates
template <typename Scalar>
NacaResult<Scalar> naca4(int ides, int nside) {
    NacaResult<Scalar> result;

    // Parse digits
    int n4 = ides / 1000;
    int n3 = (ides - n4 * 1000) / 100;
    int n2 = (ides - n4 * 1000 - n3 * 100) / 10;
    int n1 = ides - n4 * 1000 - n3 * 100 - n2 * 10;

    Scalar m = static_cast<Scalar>(n4) / Scalar(100.0);  // Max camber
    Scalar p = static_cast<Scalar>(n3) / Scalar(10.0);   // Camber position
    Scalar t = static_cast<Scalar>(n2 * 10 + n1) / Scalar(100.0);  // Thickness

    // TE point bunching parameter
    const Scalar an = Scalar(1.5);
    const Scalar anp = an + Scalar(1.0);

    // Generate x-coordinates and thickness/camber distributions
    VectorX<Scalar> xx(nside), yt(nside), yc(nside);

    for (int i = 0; i < nside; ++i) {
        Scalar frac = static_cast<Scalar>(i) / static_cast<Scalar>(nside - 1);
        if (i == nside - 1) {
            xx(i) = Scalar(1.0);
        } else {
            using std::pow;
            xx(i) = Scalar(1.0) - anp * frac * pow(Scalar(1.0) - frac, an) -
                    pow(Scalar(1.0) - frac, anp);
        }

        // Thickness distribution
        using std::sqrt;
        yt(i) = (Scalar(0.29690) * sqrt(xx(i)) - Scalar(0.12600) * xx(i) -
                 Scalar(0.35160) * xx(i) * xx(i) +
                 Scalar(0.28430) * xx(i) * xx(i) * xx(i) -
                 Scalar(0.10150) * xx(i) * xx(i) * xx(i) * xx(i)) *
                t / Scalar(0.20);

        // Camber distribution
        if (p > Scalar(0.0)) {
            if (xx(i) < p) {
                yc(i) = m / (p * p) *
                        (Scalar(2.0) * p * xx(i) - xx(i) * xx(i));
            } else {
                yc(i) = m / ((Scalar(1.0) - p) * (Scalar(1.0) - p)) *
                        ((Scalar(1.0) - Scalar(2.0) * p) +
                         Scalar(2.0) * p * xx(i) - xx(i) * xx(i));
            }
        } else {
            yc(i) = Scalar(0.0);
        }
    }

    // Build airfoil coordinates (upper surface from TE to LE, then lower)
    int nb = 2 * nside - 1;
    result.x.resize(nb);
    result.y.resize(nb);

    // Upper surface (from TE to LE, reversed order)
    int ib = 0;
    for (int i = nside - 1; i >= 0; --i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) + yt(i);
        ++ib;
    }

    // Lower surface (from LE+1 to TE)
    for (int i = 1; i < nside; ++i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) - yt(i);
        ++ib;
    }

    result.n = nb;
    result.name = "NACA " + std::to_string(n4) + std::to_string(n3) +
                  std::to_string(n2) + std::to_string(n1);

    return result;
}

/// @brief Generate 4-digit NACA airfoil from parameters.
/// @param m Maximum camber (fraction of chord)
/// @param p Position of max camber (fraction of chord)
/// @param t Maximum thickness (fraction of chord)
/// @param nside Number of points per side
/// @return NacaResult with airfoil coordinates
template <typename Scalar>
NacaResult<Scalar> naca4_params(Scalar m, Scalar p, Scalar t, int nside) {
    NacaResult<Scalar> result;

    // TE point bunching parameter
    const Scalar an = Scalar(1.5);
    const Scalar anp = an + Scalar(1.0);

    // Generate x-coordinates and thickness/camber distributions
    VectorX<Scalar> xx(nside), yt(nside), yc(nside);

    for (int i = 0; i < nside; ++i) {
        Scalar frac = static_cast<Scalar>(i) / static_cast<Scalar>(nside - 1);
        if (i == nside - 1) {
            xx(i) = Scalar(1.0);
        } else {
            using std::pow;
            xx(i) = Scalar(1.0) - anp * frac * pow(Scalar(1.0) - frac, an) -
                    pow(Scalar(1.0) - frac, anp);
        }

        // Thickness distribution
        using std::sqrt;
        yt(i) = (Scalar(0.29690) * sqrt(xx(i)) - Scalar(0.12600) * xx(i) -
                 Scalar(0.35160) * xx(i) * xx(i) +
                 Scalar(0.28430) * xx(i) * xx(i) * xx(i) -
                 Scalar(0.10150) * xx(i) * xx(i) * xx(i) * xx(i)) *
                t / Scalar(0.20);

        // Camber distribution
        if (p > Scalar(0.0)) {
            if (xx(i) < p) {
                yc(i) = m / (p * p) *
                        (Scalar(2.0) * p * xx(i) - xx(i) * xx(i));
            } else {
                yc(i) = m / ((Scalar(1.0) - p) * (Scalar(1.0) - p)) *
                        ((Scalar(1.0) - Scalar(2.0) * p) +
                         Scalar(2.0) * p * xx(i) - xx(i) * xx(i));
            }
        } else {
            yc(i) = Scalar(0.0);
        }
    }

    // Build airfoil coordinates
    int nb = 2 * nside - 1;
    result.x.resize(nb);
    result.y.resize(nb);

    int ib = 0;
    for (int i = nside - 1; i >= 0; --i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) + yt(i);
        ++ib;
    }
    for (int i = 1; i < nside; ++i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) - yt(i);
        ++ib;
    }

    result.n = nb;
    result.name = "NACA airfoil";

    return result;
}

/// @brief Generate 5-digit NACA airfoil from designation (e.g., 23012).
/// @param ides 5-digit NACA designation
/// @param nside Number of points per side
/// @return NacaResult with airfoil coordinates, empty name if invalid
template <typename Scalar>
NacaResult<Scalar> naca5(int ides, int nside) {
    NacaResult<Scalar> result;

    // Parse digits
    int n5 = ides / 10000;
    int n4 = (ides - n5 * 10000) / 1000;
    int n3 = (ides - n5 * 10000 - n4 * 1000) / 100;
    int n2 = (ides - n5 * 10000 - n4 * 1000 - n3 * 100) / 10;
    int n1 = ides - n5 * 10000 - n4 * 1000 - n3 * 100 - n2 * 10;

    int n543 = 100 * n5 + 10 * n4 + n3;

    // Look up camber line parameters
    Scalar m_cl, c_cl;  // Camber line m and C parameters
    if (n543 == 210) {
        m_cl = Scalar(0.0580);
        c_cl = Scalar(361.4);
    } else if (n543 == 220) {
        m_cl = Scalar(0.1260);
        c_cl = Scalar(51.64);
    } else if (n543 == 230) {
        m_cl = Scalar(0.2025);
        c_cl = Scalar(15.957);
    } else if (n543 == 240) {
        m_cl = Scalar(0.2900);
        c_cl = Scalar(6.643);
    } else if (n543 == 250) {
        m_cl = Scalar(0.3910);
        c_cl = Scalar(3.230);
    } else {
        // Invalid designation
        result.n = 0;
        return result;
    }

    Scalar t = static_cast<Scalar>(n2 * 10 + n1) / Scalar(100.0);

    // TE point bunching parameter
    const Scalar an = Scalar(1.5);
    const Scalar anp = an + Scalar(1.0);

    VectorX<Scalar> xx(nside), yt(nside), yc(nside);

    for (int i = 0; i < nside; ++i) {
        Scalar frac = static_cast<Scalar>(i) / static_cast<Scalar>(nside - 1);
        if (i == nside - 1) {
            xx(i) = Scalar(1.0);
        } else {
            using std::pow;
            xx(i) = Scalar(1.0) - anp * frac * pow(Scalar(1.0) - frac, an) -
                    pow(Scalar(1.0) - frac, anp);
        }

        // Thickness distribution
        using std::sqrt;
        yt(i) = (Scalar(0.29690) * sqrt(xx(i)) - Scalar(0.12600) * xx(i) -
                 Scalar(0.35160) * xx(i) * xx(i) +
                 Scalar(0.28430) * xx(i) * xx(i) * xx(i) -
                 Scalar(0.10150) * xx(i) * xx(i) * xx(i) * xx(i)) *
                t / Scalar(0.20);

        // 5-digit camber distribution
        if (xx(i) < m_cl) {
            yc(i) = (c_cl / Scalar(6.0)) *
                    (xx(i) * xx(i) * xx(i) -
                     Scalar(3.0) * m_cl * xx(i) * xx(i) +
                     m_cl * m_cl * (Scalar(3.0) - m_cl) * xx(i));
        } else {
            yc(i) = (c_cl / Scalar(6.0)) * m_cl * m_cl * m_cl *
                    (Scalar(1.0) - xx(i));
        }
    }

    // Build airfoil coordinates
    int nb = 2 * nside - 1;
    result.x.resize(nb);
    result.y.resize(nb);

    int ib = 0;
    for (int i = nside - 1; i >= 0; --i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) + yt(i);
        ++ib;
    }
    for (int i = 1; i < nside; ++i) {
        result.x(ib) = xx(i);
        result.y(ib) = yc(i) - yt(i);
        ++ib;
    }

    result.n = nb;
    result.name = "NACA " + std::to_string(n5) + std::to_string(n4) +
                  std::to_string(n3) + std::to_string(n2) + std::to_string(n1);

    return result;
}

/// @brief Locates leading edge spline-parameter value.
///
/// The defining condition is that the surface tangent is normal to
/// the chord line connecting the LE point and the TE point.
///
/// @param x X coordinates
/// @param xp X spline derivatives
/// @param y Y coordinates
/// @param yp Y spline derivatives
/// @param s Arc length array
/// @return Arc length parameter at leading edge
template <typename Scalar>
Scalar lefind(const VectorX<Scalar>& x, const VectorX<Scalar>& xp,
              const VectorX<Scalar>& y, const VectorX<Scalar>& yp,
              const VectorX<Scalar>& s) {
    int n = x.size();

    // Convergence tolerance
    Scalar dseps = (s(n - 1) - s(0)) * Scalar(1.0e-5);

    // Trailing edge point coordinates
    Scalar xte = Scalar(0.5) * (x(0) + x(n - 1));
    Scalar yte = Scalar(0.5) * (y(0) + y(n - 1));

    // Get first guess for SLE
    int ile = 2;
    for (int i = 2; i < n - 2; ++i) {
        Scalar dxte = x(i) - xte;
        Scalar dyte = y(i) - yte;
        Scalar dx = x(i + 1) - x(i);
        Scalar dy = y(i + 1) - y(i);
        Scalar dotp = dxte * dx + dyte * dy;
        if (dotp < Scalar(0.0)) {
            ile = i;
            break;
        }
    }

    Scalar sle = s(ile);

    // Check for sharp LE case
    if (s(ile) == s(ile - 1)) {
        return sle;
    }

    // Newton iteration to get exact SLE value
    for (int iter = 0; iter < 50; ++iter) {
        Scalar xle = seval(sle, x, xp, s);
        Scalar yle = seval(sle, y, yp, s);
        Scalar dxds = deval(sle, x, xp, s);
        Scalar dyds = deval(sle, y, yp, s);
        Scalar dxdd = d2val(sle, x, xp, s);
        Scalar dydd = d2val(sle, y, yp, s);

        Scalar xchord = xle - xte;
        Scalar ychord = yle - yte;

        // Drive dot product between chord line and LE tangent to zero
        Scalar res = xchord * dxds + ychord * dyds;
        Scalar ress = dxds * dxds + dyds * dyds + xchord * dxdd + ychord * dydd;

        // Newton delta for SLE
        Scalar dsle = -res / ress;

        using std::abs;
        using std::max;
        using std::min;
        Scalar chord_approx = abs(xchord + ychord);
        dsle = max(dsle, Scalar(-0.02) * chord_approx);
        dsle = min(dsle, Scalar(0.02) * chord_approx);
        sle = sle + dsle;

        if (abs(dsle) < dseps) {
            return sle;
        }
    }

    // Failed to converge, return best guess
    return s(ile);
}

/// @brief Geometric parameters of an airfoil.
template <typename Scalar>
struct GeomParams {
    Scalar sle;      ///< Arc length at leading edge
    Scalar xle;      ///< X coordinate of leading edge
    Scalar yle;      ///< Y coordinate of leading edge
    Scalar xte;      ///< X coordinate of trailing edge
    Scalar yte;      ///< Y coordinate of trailing edge
    Scalar chord;    ///< Chord length
    Scalar radle;    ///< Leading edge radius
    Scalar angte;    ///< Trailing edge angle [rad]
};

/// @brief Calculate geometric parameters of an airfoil.
/// @param x X coordinates
/// @param xp X spline derivatives
/// @param y Y coordinates
/// @param yp Y spline derivatives
/// @param s Arc length array
/// @return Geometric parameters
template <typename Scalar>
GeomParams<Scalar> geopar(const VectorX<Scalar>& x, const VectorX<Scalar>& xp,
                          const VectorX<Scalar>& y, const VectorX<Scalar>& yp,
                          const VectorX<Scalar>& s) {
    GeomParams<Scalar> result;
    int n = x.size();

    result.sle = lefind(x, xp, y, yp, s);

    result.xle = seval(result.sle, x, xp, s);
    result.yle = seval(result.sle, y, yp, s);
    result.xte = Scalar(0.5) * (x(0) + x(n - 1));
    result.yte = Scalar(0.5) * (y(0) + y(n - 1));

    Scalar dx = result.xte - result.xle;
    Scalar dy = result.yte - result.yle;
    result.chord = std::sqrt(dx * dx + dy * dy);

    // Leading edge curvature and radius
    Scalar curvle = curv(result.sle, x, xp, y, yp, s);
    result.radle = Scalar(0.0);
    if (std::abs(curvle) > Scalar(0.001) * (s(n - 1) - s(0))) {
        result.radle = Scalar(1.0) / curvle;
    }

    // Trailing edge angle
    Scalar ang1 = std::atan2(-yp(0), -xp(0));
    Scalar ang2 = std::atan2(yp(n - 1), xp(n - 1));
    result.angte = ang2 - ang1;

    return result;
}

/// @brief Normalize airfoil to unit chord with LE at origin.
/// @param x X coordinates (modified in place)
/// @param y Y coordinates (modified in place)
/// @param s Arc length array (modified in place)
template <typename Scalar>
void normalize(VectorX<Scalar>& x, VectorX<Scalar>& y, VectorX<Scalar>& s) {
    int n = x.size();

    // Calculate arc length and splines
    scalc(x, y, s);
    VectorX<Scalar> xp, yp;
    segspl(x, s, xp);
    segspl(y, s, yp);

    // Find leading edge
    Scalar sle = lefind(x, xp, y, yp, s);

    // Get LE and TE coordinates
    Scalar xmin = seval(sle, x, xp, s);
    Scalar ymin = seval(sle, y, yp, s);
    Scalar xmax = Scalar(0.5) * (x(0) + x(n - 1));

    // Scale factor
    Scalar fudge = Scalar(1.0) / (xmax - xmin);

    // Normalize
    for (int i = 0; i < n; ++i) {
        x(i) = (x(i) - xmin) * fudge;
        y(i) = (y(i) - ymin) * fudge;
        s(i) = s(i) * fudge;
    }
}

/// @brief Calculate normal unit vectors at airfoil panel nodes.
///
/// Normal vectors point outward from the airfoil surface.
/// The normals are computed as perpendicular to the tangent (from splines).
///
/// @param x X coordinates
/// @param y Y coordinates
/// @param s Arc length array
/// @param nx X component of normal (output)
/// @param ny Y component of normal (output)
template <typename Scalar>
void ncalc(const VectorX<Scalar>& x, const VectorX<Scalar>& y,
           const VectorX<Scalar>& s, VectorX<Scalar>& nx, VectorX<Scalar>& ny) {
    int n = x.size();
    if (n <= 1) return;

    nx.resize(n);
    ny.resize(n);

    // Compute spline derivatives (use nx, ny as temp storage first)
    segspl(x, s, nx);
    segspl(y, s, ny);

    // Convert tangent to outward normal: rotate 90 degrees clockwise
    // tangent = (dx/ds, dy/ds)
    // normal = (dy/ds, -dx/ds) normalized
    for (int i = 0; i < n; ++i) {
        Scalar sx = ny(i);   // dy/ds
        Scalar sy = -nx(i);  // -dx/ds
        using std::sqrt;
        Scalar smod = sqrt(sx * sx + sy * sy);
        nx(i) = sx / smod;
        ny(i) = sy / smod;
    }

    // Average normal vectors at corner points (where s values are equal)
    for (int i = 0; i < n - 1; ++i) {
        if (s(i) == s(i + 1)) {
            Scalar sx = Scalar(0.5) * (nx(i) + nx(i + 1));
            Scalar sy = Scalar(0.5) * (ny(i) + ny(i + 1));
            using std::sqrt;
            Scalar smod = sqrt(sx * sx + sy * sy);
            nx(i) = sx / smod;
            ny(i) = sy / smod;
            nx(i + 1) = nx(i);
            ny(i + 1) = ny(i);
        }
    }
}

/// @brief Trailing edge geometry parameters.
template <typename Scalar>
struct TEParams {
    Scalar dxte;     ///< TE gap x component (x(1) - x(n))
    Scalar dyte;     ///< TE gap y component (y(1) - y(n))
    Scalar dste;     ///< Total TE gap (magnitude)
    Scalar ante;     ///< Normal-projected TE gap area
    Scalar aste;     ///< Streamwise-projected TE gap area
    Scalar scs;      ///< cos-like TE direction factor
    Scalar sds;      ///< sin-like TE direction factor
    bool sharp;      ///< True if TE is effectively sharp
};

/// @brief Calculate trailing edge geometry parameters.
///
/// Computes TE gap vector, projected areas, and direction factors
/// used in panel method for TE panel treatment.
///
/// @param x X coordinates
/// @param xp X spline derivatives
/// @param y Y coordinates
/// @param yp Y spline derivatives
/// @param chord Airfoil chord length
/// @return TEParams structure with TE geometry
template <typename Scalar>
TEParams<Scalar> tecalc(const VectorX<Scalar>& x, const VectorX<Scalar>& xp,
                        const VectorX<Scalar>& y, const VectorX<Scalar>& yp,
                        Scalar chord) {
    TEParams<Scalar> te;
    int n = x.size();

    // TE base vector (gap from last point to first point)
    te.dxte = x(0) - x(n - 1);
    te.dyte = y(0) - y(n - 1);

    // TE bisector direction from average of end tangents
    // Note: tangent at first point points into airfoil, so negate
    Scalar dxs = Scalar(0.5) * (-xp(0) + xp(n - 1));
    Scalar dys = Scalar(0.5) * (-yp(0) + yp(n - 1));

    // Normal and streamwise projected TE gap areas
    te.ante = dxs * te.dyte - dys * te.dxte;
    te.aste = dxs * te.dxte + dys * te.dyte;

    // Total TE gap
    using std::sqrt;
    te.dste = sqrt(te.dxte * te.dxte + te.dyte * te.dyte);

    // Sharp TE threshold
    te.sharp = te.dste < Scalar(0.0001) * chord;

    if (te.sharp) {
        te.scs = Scalar(1.0);
        te.sds = Scalar(0.0);
    } else {
        te.scs = te.ante / te.dste;
        te.sds = te.aste / te.dste;
    }

    return te;
}

}  // namespace xfoil

#endif  // XFOIL_GEOMETRY_H_
