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
#include "cppfoil/core/geometry.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "cppfoil/core/inviscid.h"
#include "cppfoil/core/paneling.h"
#include "cppfoil/data/state.h"

namespace cppfoil {

namespace {

// Helper: degrees to radians
constexpr double DEG2RAD = M_PI / 180.0;

// Helper: cos of angle in degrees
inline double cosd(double deg) { return std::cos(deg * DEG2RAD); }

// Helper: sin of angle in degrees
inline double sind(double deg) { return std::sin(deg * DEG2RAD); }

}  // anonymous namespace

void NacaPoints(Geom& geom, const std::string& digits) {
    geom.name = "NACA " + digits;

    const int N = 100;      // Points per side
    const double te = 1.5;  // Trailing-edge bunching factor

    // Generate x distribution (bunched towards TE)
    VectorXd f = VectorXd::LinSpaced(N + 1, 0.0, 1.0);
    VectorXd x(N + 1);
    for (int i = 0; i <= N; ++i) {
        double fi = f(i);
        x(i) = 1.0 - (te + 1.0) * fi * std::pow(1.0 - fi, te) -
               std::pow(1.0 - fi, te + 1.0);
    }

    // Thickness distribution (finite TE gap)
    // t = 0.2969*sqrt(x) - 0.126*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4
    VectorXd t(N + 1);
    for (int i = 0; i <= N; ++i) {
        double xi = x(i);
        t(i) = 0.2969 * std::sqrt(xi) - 0.126 * xi - 0.3516 * xi * xi +
               0.2843 * xi * xi * xi - 0.1015 * xi * xi * xi * xi;
    }

    // Scale by max thickness
    if (digits.length() < 2) {
        throw std::invalid_argument(
            "NACA digits must have at least 2 characters");
    }
    double tmax = std::stod(digits.substr(digits.length() - 2)) * 0.01;
    t = t * (tmax / 0.2);

    // Camber distribution
    VectorXd c = VectorXd::Zero(N + 1);

    if (digits.length() == 4) {
        // 4-digit series
        double m = (digits[0] - '0') * 0.01;
        double p = (digits[1] - '0') * 0.1;

        if (p > 0) {
            for (int i = 0; i <= N; ++i) {
                double xi = x(i);
                if (xi < p) {
                    c(i) = m / (p * p) * (2.0 * p * xi - xi * xi);
                } else {
                    c(i) = m / ((1.0 - p) * (1.0 - p)) *
                           ((1.0 - 2.0 * p) + 2.0 * p * xi - xi * xi);
                }
            }
        }
    } else if (digits.length() == 5) {
        // 5-digit series
        if (digits[0] != '2' || digits[2] != '0') {
            throw std::invalid_argument(
                "5-digit NACA must begin with 2X0, X in 1-5");
        }
        int n = digits[1] - '0';
        if (n < 1 || n > 5) {
            throw std::invalid_argument(
                "5-digit NACA second digit must be 1-5");
        }

        // m and cc values for 5-digit series
        const double mv[] = {0.058, 0.126, 0.2025, 0.29, 0.391};
        const double cv[] = {361.4, 51.64, 15.957, 6.643, 3.23};
        double m_val = mv[n - 1];
        double cc = cv[n - 1];

        for (int i = 0; i <= N; ++i) {
            double xi = x(i);
            if (xi <= m_val) {
                c(i) = (cc / 6.0) * (xi * xi * xi - 3.0 * m_val * xi * xi +
                                     m_val * m_val * (3.0 - m_val) * xi);
            } else {
                c(i) = (cc / 6.0) * m_val * m_val * m_val * (1.0 - xi);
            }
        }
    } else {
        throw std::invalid_argument("Provide 4 or 5 NACA digits");
    }

    // Upper and lower surfaces
    VectorXd zu = c + t;
    VectorXd zl = c - t;

    // Concatenate: from TE lower (reversed) to TE upper
    int npoints = 2 * N + 1;
    geom.npoint = npoints;
    geom.xpoint.resize(2, npoints);

    // Lower surface (reversed, from TE to LE)
    for (int i = 0; i <= N; ++i) {
        geom.xpoint(0, i) = x(N - i);
        geom.xpoint(1, i) = zl(N - i);
    }
    // Upper surface (from LE+1 to TE)
    for (int i = 1; i <= N; ++i) {
        geom.xpoint(0, N + i) = x(i);
        geom.xpoint(1, N + i) = zu(i);
    }

    // Compute chord
    geom.chord = geom.xpoint.row(0).maxCoeff() - geom.xpoint.row(0).minCoeff();
}

void SetCoords(Geom& geom, const MatrixXd& coords) {
    geom.npoint = static_cast<int>(coords.cols());
    geom.xpoint = coords;
    geom.chord = coords.row(0).maxCoeff() - coords.row(0).minCoeff();
    geom.name = "custom";
}

TEInfo GetTEInfo(const MatrixXd& x) {
    TEInfo info;
    int N = static_cast<int>(x.cols());

    // Tangent directions at TE (matching Python: pointing toward TE from
    // interior) t1: lower surface tangent (from node 1 toward node 0, i.e.,
    // toward lower TE) t2: upper surface tangent (from node N-2 toward node
    // N-1, i.e., toward upper TE)
    Vector2d t1 = x.col(0) - x.col(1);          // Lower: toward TE
    Vector2d t2 = x.col(N - 1) - x.col(N - 2);  // Upper: toward TE
    t1.normalize();
    t2.normalize();

    // Average tangent (bisector)
    info.t = 0.5 * (t1 + t2);
    double tnorm = info.t.norm();
    if (tnorm > 1e-10) {
        info.t /= tnorm;
    }

    // Gap vector: from lower TE to upper TE
    Vector2d s = x.col(N - 1) - x.col(0);

    // TE gap (cross product of gap with bisector)
    info.hTE = -s(0) * info.t(1) + s(1) * info.t(0);

    // dtdx: sin of angle between t1 and t2
    info.dtdx = t1(0) * t2(1) - t2(0) * t1(1);

    // TE panel strength factors
    // p: unit vector along TE panel (from lower to upper)
    double s_norm = s.norm();
    if (s_norm > 1e-10) {
        Vector2d p = s / s_norm;
        // tcp: |t cross p|, used for TE source panel strength
        info.tcp = std::abs(info.t(0) * p(1) - info.t(1) * p(0));
        // tdp: t dot p, used for TE vortex panel strength
        info.tdp = info.t.dot(p);
    } else {
        info.tcp = 0.0;
        info.tdp = 0.0;
    }

    return info;
}

void GeomFlap(FoilState& state, const Vector2d& xzhinge, double eta) {
    // Deploy a flap at hinge location with given deflection angle
    // Matches mfoil.mgeom_flap

    MatrixXd X = state.geom.xpoint;  // Make a copy
    int N = static_cast<int>(X.cols());
    double xh = xzhinge(0);  // x hinge location

    // Build rotation matrix (positive eta = clockwise when viewed from above)
    Eigen::Matrix2d R;
    R << cosd(eta), sind(eta), -sind(eta), cosd(eta);

    // First pass: identify which points are on the flap (x > xh) before
    // rotation
    std::vector<bool> isFlap(N, false);
    for (int i = 0; i < N; ++i) {
        isFlap[i] = (X(0, i) > xh);
    }

    // Second pass: rotate flap points around hinge
    for (int i = 0; i < N; ++i) {
        if (isFlap[i]) {
            Vector2d p = X.col(i) - xzhinge;
            X.col(i) = xzhinge + R * p;
        }
    }

    // Third pass: identify which points to keep
    // Keep points that: (a) were not flap points, OR (b) are still right of
    // hinge after rotation
    std::vector<int> keepIndices;
    for (int i = 0; i < N; ++i) {
        if (!isFlap[i] || X(0, i) >= xh) {
            keepIndices.push_back(i);
        }
    }

    // Build new point array with only kept points
    int Nnew = static_cast<int>(keepIndices.size());
    MatrixXd Xnew(2, Nnew);
    for (int i = 0; i < Nnew; ++i) {
        Xnew.col(i) = X.col(keepIndices[i]);
    }

    // Store back
    state.geom.xpoint = Xnew;
    state.geom.npoint = Nnew;

    // Repanel if panels exist
    if (state.foil.N > 0) {
        MakePanels(state.geom, state.foil.N - 1, state.foil);
    }

    // Clear any existing solution
    ClearSolution(state);
}

void GeomAddCamber(FoilState& state, const MatrixXd& xzcamb) {
    // Add camber to airfoil from given camber line increment
    // Matches mfoil.mgeom_addcamber

    // Ensure xzcamb is [2 x Nc] (if given as [Nc x 2], transpose)
    MatrixXd camber = xzcamb;
    if (camber.rows() > camber.cols()) {
        camber.transposeInPlace();
    }

    MatrixXd& X = state.geom.xpoint;
    int N = static_cast<int>(X.cols());
    int Nc = static_cast<int>(camber.cols());

    // For each point, interpolate the camber delta at that x location
    // Using linear interpolation for simplicity
    for (int i = 0; i < N; ++i) {
        double x = X(0, i);

        // Find bracketing indices in camber data
        double dz = 0.0;
        if (x <= camber(0, 0)) {
            dz = camber(1, 0);
        } else if (x >= camber(0, Nc - 1)) {
            dz = camber(1, Nc - 1);
        } else {
            // Linear search for bracket (camber data assumed sorted)
            for (int j = 0; j < Nc - 1; ++j) {
                if (x >= camber(0, j) && x <= camber(0, j + 1)) {
                    // Linear interpolation
                    double t =
                        (x - camber(0, j)) / (camber(0, j + 1) - camber(0, j));
                    dz = (1.0 - t) * camber(1, j) + t * camber(1, j + 1);
                    break;
                }
            }
        }

        // Add camber increment to z coordinate
        X(1, i) += dz;
    }

    // Repanel if panels exist
    if (state.foil.N > 0) {
        MakePanels(state.geom, state.foil.N - 1, state.foil);
    }

    // Clear any existing solution
    ClearSolution(state);
}

void GeomDerotate(FoilState& state) {
    // Derotate airfoil to make chord line horizontal
    // Matches mfoil.mgeom_derotate

    MatrixXd& X = state.geom.xpoint;
    int N = static_cast<int>(X.cols());

    // Find leading edge (minimum x)
    int iLE = 0;
    double xmin = X(0, 0);
    for (int i = 1; i < N; ++i) {
        if (X(0, i) < xmin) {
            xmin = X(0, i);
            iLE = i;
        }
    }
    Vector2d xLE = X.col(iLE);

    // Trailing edge "point" = average of first and last points
    Vector2d xTE = 0.5 * (X.col(0) + X.col(N - 1));

    // Rotation angle to make TE-LE horizontal
    double theta = std::atan2(xTE(1) - xLE(1), xTE(0) - xLE(0));

    // Rotation matrix (negative angle to derotate)
    Eigen::Matrix2d R;
    R << std::cos(theta), std::sin(theta), -std::sin(theta), std::cos(theta);

    // Rotate all points about LE
    for (int i = 0; i < N; ++i) {
        Vector2d p = X.col(i) - xLE;
        X.col(i) = xLE + R * p;
    }

    // Repanel if panels exist
    if (state.foil.N > 0) {
        MakePanels(state.geom, state.foil.N - 1, state.foil);
    }

    // Clear any existing solution
    ClearSolution(state);
}

}  // namespace cppfoil
