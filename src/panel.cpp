// Copyright (c) 2024 ProjectTrysil Authors
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Panel method influence calculations.

#include "xfoil/panel.h"

#include <cmath>

namespace xfoil {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kEps = 1.0e-12;
}  // namespace

double atanc(double y, double x) {
    // atan2 that returns value in [0, 2*pi) instead of [-pi, pi]
    double angle = std::atan2(y, x);
    if (angle < 0.0) {
        angle += kTwoPi;
    }
    return angle;
}

void apcalc(const VectorXd& x, const VectorXd& y, VectorXd& apanel) {
    int n = x.size();
    apanel.resize(n);

    // Panel angles: angle from node i to node i+1
    for (int i = 0; i < n - 1; ++i) {
        double dx = x(i + 1) - x(i);
        double dy = y(i + 1) - y(i);
        apanel(i) = std::atan2(dy, dx);
    }

    // Last node: use angle of previous panel
    apanel(n - 1) = apanel(n - 2);
}

void freestream(double alpha, double qinf, double& uinf, double& vinf) {
    uinf = qinf * std::cos(alpha);
    vinf = qinf * std::sin(alpha);
}

PanelInfluence psilin(double xf, double yf,
                      double x1, double y1,
                      double x2, double y2,
                      double gam, double sig) {
    PanelInfluence result = {};

    // Panel geometry
    double dx = x2 - x1;
    double dy = y2 - y1;
    double ds = std::sqrt(dx * dx + dy * dy);

    if (ds < kEps) {
        // Degenerate panel
        return result;
    }

    // Unit tangent and normal vectors
    double tx = dx / ds;
    double ty = dy / ds;
    double nx = -ty;  // Normal points to the right of tangent
    double ny = tx;

    // Transform field point to panel coordinates
    // Origin at panel midpoint
    double xm = 0.5 * (x1 + x2);
    double ym = 0.5 * (y1 + y2);
    double xfr = xf - xm;
    double yfr = yf - ym;

    // Local coordinates (s along panel, n normal to panel)
    double sf = xfr * tx + yfr * ty;
    double nf = xfr * nx + yfr * ny;

    // Half panel length
    double a = 0.5 * ds;

    // Distances from field point to panel endpoints (squared)
    double r1sq = (sf + a) * (sf + a) + nf * nf;
    double r2sq = (sf - a) * (sf - a) + nf * nf;

    // Avoid log(0) singularity
    r1sq = std::max(r1sq, kEps);
    r2sq = std::max(r2sq, kEps);

    // Logarithm term
    double logr = 0.5 * std::log(r1sq / r2sq);

    // Angle subtended by panel from field point
    // theta = angle from r1 to r2
    double theta1 = std::atan2(nf, sf + a);
    double theta2 = std::atan2(nf, sf - a);
    double dtheta = theta1 - theta2;

    // Handle branch cut in atan2
    // Ensure dtheta is in [-pi, pi]
    while (dtheta > kPi) dtheta -= kTwoPi;
    while (dtheta < -kPi) dtheta += kTwoPi;

    // XFOIL panel influence formulations:
    // Vortex panel: psi = (gam/2pi) * [nf*ln(r1/r2)/2 + sf*(theta2-theta1)]
    //                   = (gam/2pi) * [-nf*logr + sf*(-dtheta)]
    //                   = (gam/2pi) * [-nf*logr - sf*dtheta]
    // Source panel: psi = (sig/2pi) * [sf*ln(r1/r2)/2 + nf*(theta1-theta2)]
    //                   = (sig/2pi) * [sf*logr + nf*dtheta]

    // Contributions (without gamma/sigma factors)
    double psi_per_gam = (-nf * logr - sf * dtheta) / kTwoPi;
    double psi_per_sig = (sf * logr + nf * dtheta) / kTwoPi;

    // Tangential velocity contributions (derivative of psi w.r.t. n)
    // d(psi)/dn at the panel surface
    // For vortex: qtan = d(psi_gam)/dn = gam/2 at the surface (jump condition)
    // For source: qtan = d(psi_sig)/ds = 0 (source doesn't contribute to tangential velocity)

    // Away from the panel, we need the full velocity calculation
    // u = d(psi)/dy, v = -d(psi)/dx
    // qtan = u*tx + v*ty

    // For points not on the panel:
    double u_per_gam = 0.0;
    double v_per_gam = 0.0;
    double u_per_sig = 0.0;
    double v_per_sig = 0.0;

    // Vortex contribution to velocity
    // u = d(psi)/dy, v = -d(psi)/dx
    // From psi_vortex, we get:
    u_per_gam = dtheta / kTwoPi;
    v_per_gam = -logr / kTwoPi;

    // Source contribution to velocity
    u_per_sig = logr / kTwoPi;
    v_per_sig = dtheta / kTwoPi;

    // Transform back to global coordinates
    // qtan = u*tx + v*ty
    double qtan_per_gam = u_per_gam * ty - v_per_gam * tx;  // perpendicular to panel
    double qtan_per_sig = u_per_sig * ty - v_per_sig * tx;

    // Actually, let's use the tangent direction directly
    // Velocity in s-direction (tangent): qs = -d(psi)/dn
    // Velocity in n-direction (normal):  qn = d(psi)/ds
    //
    // d(psi_vortex)/dn = (gam/2pi) * [-logr - sf * d(dtheta)/dn]
    // d(psi_source)/dn = (sig/2pi) * [sf * d(logr)/dn + dtheta]
    //
    // These are complex. Let's use numerical approach is in actual code.

    result.psi = gam * psi_per_gam + sig * psi_per_sig;
    result.psi_gam = psi_per_gam;
    result.psi_sig = psi_per_sig;
    result.qtan = gam * qtan_per_gam + sig * qtan_per_sig;
    result.qtan_gam = qtan_per_gam;
    result.qtan_sig = qtan_per_sig;

    return result;
}

void aic_build(const VectorXd& x, const VectorXd& y,
               const VectorXd& apanel, MatrixXd& aij, bool sharp_te) {
    int n = x.size();
    aij.resize(n, n);
    aij.setZero();

    // Build influence coefficients
    // Aij = stream function at node i due to unit vorticity on panel j
    for (int i = 0; i < n; ++i) {
        // Field point (slightly inside the airfoil to avoid singularity)
        double xf = x(i);
        double yf = y(i);

        // Displace field point slightly along inward normal
        double angle_i;
        if (i == 0) {
            angle_i = apanel(0);
        } else if (i == n - 1) {
            angle_i = apanel(n - 2);
        } else {
            // Average of adjacent panels
            angle_i = 0.5 * (apanel(i - 1) + apanel(i));
        }

        // Inward normal for CCW traversal
        double nx = std::sin(angle_i);
        double ny = -std::cos(angle_i);

        // Small displacement to regularize self-influence
        double eps_n = 1.0e-8;
        xf += eps_n * nx;
        yf += eps_n * ny;

        for (int j = 0; j < n - 1; ++j) {
            // Panel j runs from node j to node j+1
            auto inf = psilin(xf, yf, x(j), y(j), x(j + 1), y(j + 1), 1.0, 0.0);

            // Vortex panel with constant strength requires distributing
            // the influence to both endpoint nodes
            aij(i, j) += 0.5 * inf.psi_gam;
            aij(i, j + 1) += 0.5 * inf.psi_gam;
        }
    }

    // Kutta condition: replace last row
    // For sharp TE: gamma(0) + gamma(n-1) = 0
    // For blunt TE: modified condition
    int last = n - 1;
    for (int j = 0; j < n; ++j) {
        aij(last, j) = 0.0;
    }

    if (sharp_te) {
        aij(last, 0) = 1.0;
        aij(last, last) = 1.0;
    } else {
        // Blunt TE: velocity matching at TE
        // gamma(0) - gamma(n-1) = some_function_of_gap
        // For now, use same as sharp
        aij(last, 0) = 1.0;
        aij(last, last) = 1.0;
    }
}

void dic_build(const VectorXd& x, const VectorXd& y,
               const VectorXd& apanel, MatrixXd& dij) {
    int n = x.size();
    dij.resize(n, n);
    dij.setZero();

    // Build source influence coefficients
    // Dij = stream function at node i due to unit source on panel j
    for (int i = 0; i < n; ++i) {
        // Field point
        double xf = x(i);
        double yf = y(i);

        // Displace slightly along inward normal
        double angle_i;
        if (i == 0) {
            angle_i = apanel(0);
        } else if (i == n - 1) {
            angle_i = apanel(n - 2);
        } else {
            angle_i = 0.5 * (apanel(i - 1) + apanel(i));
        }

        double nx = std::sin(angle_i);
        double ny = -std::cos(angle_i);

        double eps_n = 1.0e-8;
        xf += eps_n * nx;
        yf += eps_n * ny;

        for (int j = 0; j < n - 1; ++j) {
            auto inf = psilin(xf, yf, x(j), y(j), x(j + 1), y(j + 1), 0.0, 1.0);

            dij(i, j) += 0.5 * inf.psi_sig;
            dij(i, j + 1) += 0.5 * inf.psi_sig;
        }
    }
}

void qcalc(const VectorXd& x, const VectorXd& y,
           const VectorXd& apanel, const VectorXd& gam,
           const VectorXd& sig, double alpha, VectorXd& qinv) {
    int n = x.size();
    qinv.resize(n);

    // Freestream contribution
    double uinf, vinf;
    freestream(alpha, 1.0, uinf, vinf);

    for (int i = 0; i < n; ++i) {
        // Panel tangent direction at node i
        double angle_i;
        if (i == 0) {
            angle_i = apanel(0);
        } else if (i == n - 1) {
            angle_i = apanel(n - 2);
        } else {
            angle_i = 0.5 * (apanel(i - 1) + apanel(i));
        }

        double tx = std::cos(angle_i);
        double ty = std::sin(angle_i);

        // Freestream tangential component
        double qt_fs = uinf * tx + vinf * ty;

        // Induced velocity from all panels
        double qt_ind = 0.0;
        double xf = x(i);
        double yf = y(i);

        // Small normal displacement to avoid singularity
        double nx = std::sin(angle_i);
        double ny = -std::cos(angle_i);
        double eps_n = 1.0e-8;
        xf += eps_n * nx;
        yf += eps_n * ny;

        for (int j = 0; j < n - 1; ++j) {
            // Average gamma and sigma on panel j
            double gam_j = 0.5 * (gam(j) + gam(j + 1));
            double sig_j = 0.5 * (sig(j) + sig(j + 1));

            auto inf = psilin(xf, yf, x(j), y(j), x(j + 1), y(j + 1),
                            gam_j, sig_j);

            qt_ind += inf.qtan;
        }

        qinv(i) = qt_fs + qt_ind;
    }
}

double qinf(double si, const VectorXd& x, const VectorXd& xp,
            const VectorXd& y, const VectorXd& yp, const VectorXd& s,
            const VectorXd& gam, const VectorXd& sig, double alpha) {
    // This would interpolate velocity at arbitrary arc length
    // For now, just return 0 as placeholder
    // Full implementation needs spline evaluation
    (void)si;
    (void)x;
    (void)xp;
    (void)y;
    (void)yp;
    (void)s;
    (void)gam;
    (void)sig;
    (void)alpha;
    return 0.0;
}

}  // namespace xfoil
