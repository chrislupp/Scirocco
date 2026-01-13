// cppfoil: C++ airfoil analysis library
// Copyright (c) 2024
// SPDX-License-Identifier: MIT

#include "cppfoil/core/wake.h"

#include <cmath>
#include <algorithm>
#include "cppfoil/core/inviscid.h"
#include "cppfoil/util/math_utils.h"

namespace cppfoil {

void StagPointFind(FoilState& state) {
    int N = state.foil.N;
    const VectorXd& gam = state.isol.gam;

    // Find where gamma changes sign (stagnation point)
    int istag = -1;
    for (int i = 0; i < N - 1; ++i) {
        if (gam(i) * gam(i + 1) <= 0) {
            istag = i;
            break;
        }
    }

    if (istag < 0) {
        // No sign change found; use minimum |gamma|
        int imin;
        gam.cwiseAbs().minCoeff(&imin);
        istag = imin;
    }

    // Linear interpolation for stagnation s location
    double g1 = gam(istag);
    double g2 = gam(istag + 1);
    double denom = g2 - g1;
    double frac = (std::abs(denom) > 1e-12) ? -g1 / denom : 0.5;
    frac = Clamp(frac, 0.0, 1.0);

    double s1 = state.foil.s(istag);
    double s2 = state.foil.s(istag + 1);
    state.isol.sstag = s1 + frac * (s2 - s1);

    // Store indices
    state.isol.Istag[0] = istag;
    state.isol.Istag[1] = istag + 1;

    // Linearization of sstag w.r.t. adjacent gammas
    if (std::abs(denom) > 1e-12) {
        double ds = s2 - s1;
        state.isol.sstag_g(0) = ds * g2 / (denom * denom);
        state.isol.sstag_g(1) = -ds * g1 / (denom * denom);
    } else {
        state.isol.sstag_g.setZero();
    }

    // Linearization of sstag w.r.t. edge velocity (ue)
    // ga_ue = dgamma/due = sign(gam) for inviscid
    double ga_ue0 = (std::abs(g1) > 1e-12) ? g1 / std::abs(g1) : 0.0;
    double ga_ue1 = (std::abs(g2) > 1e-12) ? g2 / std::abs(g2) : 0.0;
    state.isol.sstag_ue(0) = state.isol.sstag_g(0) * ga_ue0;
    state.isol.sstag_ue(1) = state.isol.sstag_g(1) * ga_ue1;

    // Stagnation point coordinates
    Vector2d x1 = state.foil.x.col(istag);
    Vector2d x2 = state.foil.x.col(istag + 1);
    state.isol.xstag = x1 + frac * (x2 - x1);

    // Set sign of edge velocity: positive on upper, negative on lower
    state.isol.sgnue.resize(N);
    for (int i = 0; i < N; ++i) {
        if (i <= istag) {
            state.isol.sgnue(i) = -1.0;  // Lower surface
        } else {
            state.isol.sgnue(i) = 1.0;   // Upper surface
        }
    }

    // Compute xi (distance from stagnation point) for foil points only
    // Wake points will be added after BuildWake is called
    state.isol.xi.resize(N);
    for (int i = 0; i < N; ++i) {
        state.isol.xi(i) = std::abs(state.foil.s(i) - state.isol.sstag);
    }
}

void ExtendXiToWake(FoilState& state) {
    // Extend xi to include wake points: xi_wake = wake.s + (foil perimeter - sstag)
    int N = state.foil.N;
    int Nw = state.wake.N;
    if (Nw == 0) return;

    // xi for wake = wake.s + distance from stag to TE (along upper surface)
    // Actually mfoil uses: wake.s - sstag, but wake.s starts at 0 at TE
    // So effectively xi_wake = wake.s + (foil.s[N-1] - sstag)
    // But simpler: xi_wake = wake.s + xi_at_TE
    double xi_TE = state.isol.xi(N - 1);  // xi at upper TE

    VectorXd xi_new(N + Nw);
    xi_new.head(N) = state.isol.xi;
    for (int i = 0; i < Nw; ++i) {
        xi_new(N + i) = xi_TE + state.wake.s(i);
    }
    state.isol.xi = xi_new;
}

void BuildWake(FoilState& state) {
    // Number of wake points
    int Nw = 32;  // Default number of wake points
    // double wakelen = state.geom.wakelen * state.geom.chord;  // TODO: use for spacing

    // TE point (average of upper and lower TE)
    int N = state.foil.N;
    Vector2d xy1 = state.foil.x.col(0);      // Lower TE
    Vector2d xyN = state.foil.x.col(N - 1);  // Upper TE
    Vector2d xte = 0.5 * (xy1 + xyN);

    // TE direction: tangent pointing downstream
    Vector2d n = xyN - xy1;  // Normal to TE gap
    Vector2d t_te(n(1), -n(0));  // Rotate 90 degrees to get tangent
    double t_mag = t_te.norm();
    if (t_mag < 1e-12) {
        t_te = Vector2d(1.0, 0.0);  // Default to +x if TE is sharp
    } else {
        t_te /= t_mag;
    }
    // Ensure tangent points downstream (positive x)
    if (t_te(0) < 0) t_te = -t_te;

    // Wake stretching ratio
    double ds0 = (state.foil.x.col(N - 1) - state.foil.x.col(N - 2)).norm();
    double stretch = 1.2;  // Geometric stretching factor

    // Generate wake points via predictor-corrector streamline tracing
    state.wake.N = Nw;
    state.wake.x.resize(2, Nw);
    state.wake.s.resize(Nw);
    state.wake.t.resize(2, Nw);

    // First wake point slightly behind TE (like mfoil)
    state.wake.x.col(0) = xte + 1e-5 * t_te * state.geom.chord;
    state.wake.s(0) = 0.0;
    state.wake.t.col(0) = t_te;  // Initial tangent

    double ds = ds0;
    double Vinf = state.oper.Vinf;
    double alpha = state.oper.alpha;

    for (int i = 1; i < Nw; ++i) {
        Vector2d x0 = state.wake.x.col(i - 1);

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
        state.wake.s(i) = state.wake.s(i - 1) + ds;

        // Stretch for next step
        ds *= stretch;
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
    // Update stagnation point based on viscous solution
    // This is called during coupled iterations to move the stagnation point
    // based on where the viscous edge velocity goes to zero

    int N = state.foil.N;
    int Nw = state.wake.N;

    // Re-find stagnation point based on inviscid gamma
    // (simplified - full version would use viscous ue)
    StagPointFind(state);

    // Update xi values for foil points (already done in StagPointFind)
    // Now update wake xi values too
    if (Nw > 0) {
        double xi_TE = state.isol.xi(N - 1);  // xi at upper TE
        // Resize if needed (StagPointFind only allocates N)
        if (state.isol.xi.size() < N + Nw) {
            VectorXd xi_new(N + Nw);
            xi_new.head(N) = state.isol.xi;
            state.isol.xi = xi_new;
        }
        for (int i = 0; i < Nw; ++i) {
            state.isol.xi(N + i) = xi_TE + state.wake.s(i);
        }
    }
}

}  // namespace cppfoil
