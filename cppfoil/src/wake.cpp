// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/wake.h"

#include <cassert>
#include <cmath>
#include <algorithm>
#include "cppfoil/core/coupled_solver.h"
#include "cppfoil/core/inviscid.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

void StagPointFind(FoilState& state) {
    int N = state.foil.N;
    const VectorXd& gam = state.isol.gam;

    // Find first index j where gamma becomes positive (matching mfoil)
    // Use a tolerance for near-zero values to handle numerical precision
    // (at the stagnation point, gamma should be exactly zero but numerics give ~1e-13)
    constexpr double gam_tol = 1e-10;

    int j = 0;
    for (; j < N; ++j) {
        // Treat values within tolerance of zero as effectively positive
        // This handles cases where gamma at stagnation is ~1e-13 due to numerics
        if (gam(j) > -gam_tol) break;
    }

    if (j == 0 || j >= N) {
        // Fallback: find minimum |gamma|
        int imin;
        gam.cwiseAbs().minCoeff(&imin);
        j = std::max(1, imin);
    }

    // Stagnation is between j-1 and j
    int i0 = j - 1;
    int i1 = j;

    // Store indices (matching mfoil: I = [j-1, j])
    state.isol.Istag[0] = i0;
    state.isol.Istag[1] = i1;

    // Linear interpolation for stagnation s location
    double g0 = gam(i0);
    double g1 = gam(i1);
    double denom = g1 - g0;
    double w1 = g1 / denom;
    double w2 = -g0 / denom;

    double s0 = state.foil.s(i0);
    double s1 = state.foil.s(i1);
    state.isol.sstag = w1 * s0 + w2 * s1;

    // Stagnation point coordinates
    state.isol.xstag = w1 * state.foil.x.col(i0) + w2 * state.foil.x.col(i1);

    // Linearization of sstag w.r.t. adjacent gammas
    double st_g1 = g1 * (s0 - s1) / (denom * denom);
    state.isol.sstag_g(0) = st_g1;
    state.isol.sstag_g(1) = -st_g1;

    // Set sign of edge velocity: negative on lower (i < j), positive on upper (i >= j)
    state.isol.sgnue.resize(N);
    for (int i = 0; i < N; ++i) {
        if (i < j) {
            state.isol.sgnue(i) = -1.0;  // Lower surface
        } else {
            state.isol.sgnue(i) = 1.0;   // Upper surface
        }
    }

    // Linearization of sstag w.r.t. edge velocity (ue)
    // ga_ue = dgamma/due = sign(gam) for inviscid
    double ga_ue0 = (std::abs(g0) > 1e-12) ? g0 / std::abs(g0) : 0.0;
    double ga_ue1 = (std::abs(g1) > 1e-12) ? g1 / std::abs(g1) : 0.0;
    state.isol.sstag_ue(0) = state.isol.sstag_g(0) * ga_ue0;
    state.isol.sstag_ue(1) = state.isol.sstag_g(1) * ga_ue1;

    // Compute xi (distance from stagnation point) for foil points only
    // Wake points will be added after BuildWake is called
    state.isol.xi.resize(N);
    for (int i = 0; i < N; ++i) {
        state.isol.xi(i) = std::abs(state.foil.s(i) - state.isol.sstag);
    }
}

void ExtendXiToWake(FoilState& state) {
    // Extend xi to include wake points
    // mfoil formula: xi = wake.s - sstag (wake.s continues from foil.s[N-1])
    int N = state.foil.N;
    int Nw = state.wake.N;
    if (Nw == 0) return;

    VectorXd xi_new(N + Nw);
    xi_new.head(N) = state.isol.xi;

    // Wake xi = wake.s - sstag (matching mfoil exactly)
    for (int i = 0; i < Nw; ++i) {
        xi_new(N + i) = state.wake.s(i) - state.isol.sstag;
    }
    state.isol.xi = xi_new;
}

namespace {
// Geometric spacing function matching mfoil's space_geom
// Returns Np points from [0, L] with first interval dx0
VectorXd SpaceGeom(double dx0, double L, int Np) {
    if (Np <= 1) {
        VectorXd x(1);
        x(0) = 0.0;
        return x;
    }

    int N = Np - 1;  // Number of intervals
    double d = L / dx0;

    // Solve for geometric ratio r using Newton iteration
    // L = dx0 * (r^N - 1) / (r - 1)
    // Initial guess from polynomial approximation
    double a = N * (N - 1.0) * (N - 2.0) / 6.0;
    double b = N * (N - 1.0) / 2.0;
    double c = N - d;
    double disc = std::max(b * b - 4.0 * a * c, 0.0);
    double r = 1.0 + (-b + std::sqrt(disc)) / (2.0 * a);

    // Newton iterations
    for (int k = 0; k < 10; ++k) {
        double rN = std::pow(r, N);
        double R = rN - 1.0 - d * (r - 1.0);
        double R_r = N * std::pow(r, N - 1) - d;
        double dr = -R / R_r;
        if (std::abs(dr) < 1e-6) break;
        r += dr;
    }

    // Generate points
    VectorXd x(Np);
    x(0) = 0.0;
    double dx = dx0;
    for (int i = 1; i < Np; ++i) {
        x(i) = x(i - 1) + dx;
        dx *= r;
    }

    return x;
}
}  // namespace

void BuildWake(FoilState& state) {
    int N = state.foil.N;
    double wakelen = state.geom.wakelen;
    double chord = state.geom.chord;

    // Number of wake points (matching mfoil exactly)
    int Nw = static_cast<int>(std::ceil(N / 10.0 + 10.0 * wakelen));

    // TE point (average of upper and lower TE)
    Vector2d xy1 = state.foil.x.col(0);      // Lower TE
    Vector2d xyN = state.foil.x.col(N - 1);  // Upper TE
    Vector2d xte = 0.5 * (xy1 + xyN);

    // TE direction: tangent pointing downstream
    Vector2d n = xyN - xy1;  // Normal to TE gap
    Vector2d t_te(n(1), -n(0));  // Rotate 90 degrees to get tangent
    // Ensure tangent points downstream (positive x)
    if (t_te(0) < 0) t_te = -t_te;

    // First nominal wake panel size (average of TE panels)
    double ds1 = 0.5 * ((state.foil.s(1) - state.foil.s(0)) +
                        (state.foil.s(N - 1) - state.foil.s(N - 2)));

    // Geometric spacing for wake
    VectorXd sv = SpaceGeom(ds1, wakelen * chord, Nw);

    // Wake s values continue from airfoil (matching mfoil)
    VectorXd sw(Nw);
    for (int i = 0; i < Nw; ++i) {
        sw(i) = state.foil.s(N - 1) + sv(i);
    }

    // Initialize wake arrays
    state.wake.N = Nw;
    state.wake.x.resize(2, Nw);
    state.wake.s = sw;
    state.wake.t.resize(2, Nw);

    // First wake point slightly behind TE (like mfoil)
    state.wake.x.col(0) = xte + 1e-5 * t_te * chord;

    double Vinf = state.oper.Vinf;
    double alpha = state.oper.alpha;

    // Generate wake points via predictor-corrector streamline tracing
    for (int i = 1; i < Nw; ++i) {
        Vector2d x0 = state.wake.x.col(i - 1);
        double ds = sv(i) - sv(i - 1);  // Distance to next point

        // Predictor: velocity at current point
        Vector2d V0 = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, x0);
        double V0mag = V0.norm();
        if (V0mag < 1e-10) V0mag = 1e-10;
        Vector2d t0 = V0 / V0mag;

        // Predictor step
        Vector2d x1_pred = x0 + ds * t0;

        // Corrector: velocity at predicted point
        Vector2d V1 = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, x1_pred);
        double V1mag = V1.norm();
        if (V1mag < 1e-10) V1mag = 1e-10;
        Vector2d t1 = V1 / V1mag;

        // Corrector step (average tangent)
        Vector2d t_avg = 0.5 * (t0 + t1);
        double tmag = t_avg.norm();
        if (tmag < 1e-10) tmag = 1e-10;
        t_avg /= tmag;

        state.wake.x.col(i) = x0 + ds * t_avg;
        state.wake.t.col(i - 1) = t_avg;
    }

    // First tangent (same as first computed tangent)
    {
        Vector2d V0 = InviscidVelocity(state.foil.x, state.isol.gam,
                                        Vinf, alpha, state.wake.x.col(0));
        double V0mag = V0.norm();
        if (V0mag < 1e-10) V0mag = 1e-10;
        state.wake.t.col(0) = V0 / V0mag;
    }
    // Last tangent
    state.wake.t.col(Nw - 1) = state.wake.t.col(Nw - 2);

    // Compute wake inviscid edge velocities
    state.isol.uewi.resize(Nw);
    state.isol.uewiref.resize(Nw, 2);

    Vector2d cs0(1.0, 0.0);  // cos(0), sin(0)
    Vector2d cs90(0.0, 1.0); // cos(90), sin(90)

    for (int i = 0; i < Nw; ++i) {
        Vector2d xi = state.wake.x.col(i);
        Vector2d ti = state.wake.t.col(i);

        // Reference velocities (alpha=0 and alpha=90)
        Vector2d V0 = InviscidVelocity(state.foil.x, state.isol.gamref.col(0),
                                        Vinf, 0.0, xi);
        Vector2d V90 = InviscidVelocity(state.foil.x, state.isol.gamref.col(1),
                                         Vinf, 90.0, xi);

        state.isol.uewiref(i, 0) = V0.dot(ti);
        state.isol.uewiref(i, 1) = V90.dot(ti);

        // Current alpha
        state.isol.uewi(i) = state.isol.uewiref(i, 0) * Cosd(alpha) +
                             state.isol.uewiref(i, 1) * Sind(alpha);
    }
}

void IdentifySurfaces(FoilState& state) {
    int N = state.foil.N;
    int Nw = state.wake.N;
    int istag = state.isol.Istag[0];

    // Clear surface indices
    for (int si = 0; si < 3; ++si) {
        state.vsol.Is[si].clear();
    }

    // Lower surface: from stagnation to TE (decreasing index)
    for (int i = istag; i >= 0; --i) {
        state.vsol.Is[0].push_back(i);
    }

    // Upper surface: from stagnation to TE (increasing index)
    for (int i = istag + 1; i < N; ++i) {
        state.vsol.Is[1].push_back(i);
    }

    // Wake: all wake points
    for (int i = 0; i < Nw; ++i) {
        state.vsol.Is[2].push_back(N + i);
    }
}

void SetWakeGap(FoilState& state) {
    int Nw = state.wake.N;
    if (Nw <= 0) return;

    // TE gap
    int N = state.foil.N;
    double hTE = (state.foil.x.col(N - 1) - state.foil.x.col(0)).norm();
    state.geom.hTE = hTE;  // Store in geometry

    // Exponential decay of gap along wake
    state.vsol.wgap.resize(Nw);
    double decay_rate = 5.0;  // Controls how fast gap closes

    for (int i = 0; i < Nw; ++i) {
        double s = state.wake.s(i);
        double smax = state.wake.s(Nw - 1);
        double xi = s / std::max(smax, 1e-10);
        state.vsol.wgap(i) = hTE * std::exp(-decay_rate * xi);
    }
}

void StagPointMove(FoilState& state) {
    // Moves the LE stagnation point on the airfoil using the global solution ue
    // This is called during coupled iterations to move the stagnation point
    // based on where the viscous edge velocity goes to zero

    int N = state.foil.N;
    int Nw = state.wake.N;
    (void)state.isol.sstag;  // Original stagnation point (unused but kept for debugging)

    // Get edge velocity from global solution
    VectorXd ue = state.glob.U.row(3).head(N).transpose();

    // Current adjacent node indices
    int I0 = state.isol.Istag[0];
    int I1 = state.isol.Istag[1];

    bool newpanel = true;

    if (ue(I1) < 0) {
        // Move stagnation point up (larger s, new panel)
        int j = I1;
        for (; j < N; ++j) {
            if (ue(j) > 0) break;
        }
        assert(j < N && "no stagnation point");
        int I1_new = j;
        // Flip sign of ue on points between old and new
        for (int k = I1; k < I1_new; ++k) {
            ue(k) *= -1.0;
        }
        I0 = I1_new - 1;
        I1 = I1_new;
    } else if (ue(I0) < 0) {
        // Move stagnation point down (smaller s, new panel)
        int j = I0;
        for (; j >= 0; --j) {
            if (ue(j) > 0) break;
        }
        assert(j > 0 && "no stagnation point");
        int I0_new = j;
        // Flip sign of ue on points between old and new
        for (int k = I0_new + 1; k <= I0; ++k) {
            ue(k) *= -1.0;
        }
        I0 = I0_new;
        I1 = I0_new + 1;
    } else {
        newpanel = false;  // Staying on the current panel
    }

    // Move point along panel by interpolating where ue = 0
    double ue0 = ue(I0);
    double ue1 = ue(I1);
    double S0 = state.foil.s(I0);
    double S1 = state.foil.s(I1);

    assert(ue0 > 0 && ue1 > 0 && "stagpoint_move: velocity error");

    double den = ue0 + ue1;
    double w1 = ue1 / den;
    double w2 = ue0 / den;
    state.isol.sstag = w1 * S0 + w2 * S1;
    state.isol.xstag = state.foil.x.col(I0) * w1 + state.foil.x.col(I1) * w2;
    state.isol.sstag_ue.resize(2);
    state.isol.sstag_ue(0) = ue1 * (S1 - S0) / (den * den);
    state.isol.sstag_ue(1) = -ue0 * (S1 - S0) / (den * den);

    // Set new xi coordinates for every point
    state.isol.xi.resize(N + Nw);
    for (int i = 0; i < N; ++i) {
        state.isol.xi(i) = std::abs(state.foil.s(i) - state.isol.sstag);
    }
    for (int i = 0; i < Nw; ++i) {
        state.isol.xi(N + i) = state.wake.s(i) - state.isol.sstag;
    }

    // If moving to a new panel, update indices and re-identify surfaces
    if (newpanel) {
        state.isol.Istag[0] = I0;
        state.isol.Istag[1] = I1;
        for (int i = 0; i <= I0; ++i) {
            state.isol.sgnue(i) = -1;
        }
        for (int i = I0 + 1; i < N; ++i) {
            state.isol.sgnue(i) = 1;
        }
        IdentifySurfaces(state);
        // Update ue in glob.U with potentially flipped signs
        state.glob.U.row(3).head(N) = ue.transpose();
        state.glob.realloc = true;
        RebuildUeM(state);
    }
}

}  // namespace cppfoil
