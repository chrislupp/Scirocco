/*
    Scirocco: C++ Airfoil Analysis Library

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
    number: AFRL-2026-0405.

    The views expressed are those of the authors and do not reflect the
    official guidance or position of the United States Government, the
    Department of Defense or of the United States Air Force.

    Statement from DoD: The Appearance of external hyperlinks does not
    constitute endorsement by the United States Department of Defense (DoD) of
    the linked websites, of the information, products, or services contained
    therein. The DoD does not exercise any editorial, security, or other
    control over the information you may find at these locations.
*/
#pragma once

namespace scirocco {

/** Shape of the onset of Tollmien-Schlichting amplification in the e^N
 *  transition model: how the amplification rate rises from zero as
 *  log10(Re_theta) passes the critical value lrc(Hk). */
enum class OnsetModel {
    /** XFOIL/mfoil cubic ramp over lrc +/- 0.1 (default, matches mfoil). */
    Cubic,
    /** Ramp fitted to Orr-Sommerfeld envelope N-factor curves of the
     *  Falkner-Skan profiles: starts at the critical Re_theta and approaches
     *  the same asymptote, so downstream amplification is unchanged. */
    Envelope
};

/** Laminar-turbulent transition treatment in the viscous solve. */
enum class TransitionModel {
    /** XFOIL/mfoil cut cell: one interval holds the transition point and
     *  stations switch from laminar to turbulent equations (default). */
    CutCell,
    /** Captured: amplification and shear stress at every node, closures
     *  blended smoothly by the amplification factor; no station switches.
     *  Prescribed angle of attack only. See captured_transition.h. */
    Captured
};

/** Under-relaxation rule for coupled Newton updates. */
enum class StepLimiter {
    /** mfoil: absolute caps of 2 on the amplification change and 0.05 on the
     *  shear-stress change at every node, 0.2 Vinf on ue, and at most a 50%
     *  decrease in theta and delta*. */
    Absolute,
    /** XFOIL: every change is normalized (amplification by 10, shear stress,
     *  theta and delta* by their values, ue by 0.25 Vinf) and the step is
     *  scaled so that no normalized change exceeds +1.5 or -0.5. */
    Relative
};

/* Parameters */
struct Param {
    /* Solver parameters */
    int verb = 0; /**< Printing verbosity level (higher = more verbose) */
    double rtol = 1e-11;      /**< Residual tolerance for Newton */
    double rtol_floor = 1e-8; /**< Also converged once the residual norm is
                                   below this and a Newton step no longer
                                   halves it (the round-off floor) */
    bool scale_coupled_system =
        true; /**< Equilibrate equations and variables */
    bool use_trust_region =
        false; /**< Experimental scaled dogleg step control */
    bool trust_update_topology = true; /**< Evaluate the complete trial state */
    double trust_radius_initial = 10.0;
    int trust_max_trials = 16;
    bool use_continuation =
        true; /**< Retry failed cold starts along an angle/Cl path */
    int continuation_max_attempts = 24;
    int continuation_max_iterations =
        800;          /**< Total extra coupled iterations */
    int niglob = 200; /**< Maximum number of global iterations (raised from
                           50 to give PTC headroom on rough intermediate
                           iterates that need extra Newton steps to ramp
                           tau back up to full Newton) */

    /* Pseudo-transient continuation (PTC) for the coupled BL Newton solve.
       Adds (1/tau) * I to the Jacobian diagonal so that early Newton steps
       are regularized: (R_U + (1/tau) I) dU = -R. tau is grown by SER
       (switched evolution-relaxation) as the residual decreases, so the
       iteration recovers full Newton near the solution and the converged
       R = 0 is unchanged. */
    bool use_ptc = true; /**< PTC on/off (use_ptc=false -> pure Newton) */
    /* Start with effectively no damping (1/tau ~ 1e-6); SER shrinks tau only
       when the residual increases. This preserves full-Newton behavior on
       well-behaved cases (the test reference) while still catching the
       divergence pattern (R: 0.97 -> 98 -> 1e15) we observed on rough
       intermediate iterates: when R jumps by a large factor, SER shrinks
       tau hard (clamped by tau_shrink_clamp) and the next step is damped. */
    double tau_init = 1e6;         /**< initial pseudo-time step */
    double tau_min = 1e-2;         /**< floor on tau (max regularization) */
    double tau_max = 1e8;          /**< cap on tau (effectively "PTC off") */
    double tau_grow_clamp = 4.0;   /**< cap per-iter tau growth ratio */
    double tau_shrink_clamp = 0.1; /**< cap per-iter tau shrink ratio */

    /* Newton globalization for the cut-cell transition model */
    StepLimiter step_limiter = StepLimiter::Absolute; /**< Under-relaxation */
    bool bl_march = true;    /**< Re-solve the boundary layer by a mixed-mode
                                  march before each Newton step (see
                                  bl_march.h); false gives mfoil's iteration */
    bool warm_start = false; /**< Scirocco::Solve starts a viscous solve
                                  from the last converged boundary layer of
                                  the same object (cold start if that
                                  fails). Off by default: near stall the
                                  result can then depend on earlier solves */

    /* Viscous parameters */
    double ncrit = 9.0;                   /**< Critical amplification factor */
    OnsetModel onset = OnsetModel::Cubic; /**< Amplification onset ramp */
    TransitionModel transition = TransitionModel::CutCell; /**< Transition */
    /** Forced transition (trip) on the upper surface at
     *  x = xtrip_upper * chord: the boundary layer turns turbulent there
     *  unless e^N transition comes first. 1 or more: free transition only
     *  (see FindTrip). With TransitionModel::Captured the closures blend
     *  over capt_L momentum thicknesses behind the trip. */
    double xtrip_upper = 1.0;
    /** Forced transition on the lower surface, as xtrip_upper. */
    double xtrip_lower = 1.0;
    /** Turbulent skin friction never below the laminar value at the same
     *  Hk and Re_theta, as in XFOIL. It acts at small Re_theta (a trip close
     *  behind the stagnation point) and just behind a separation bubble.
     *  false gives mfoil's closure. Both transition models use it. */
    bool cf_laminar_floor = true;
    double capt_L = 60.0; /**< Captured: 10-90% transition length [theta],
                               calibrated on E387 drag (UIUC Vol 5) */
    double capt_ws = 1.0; /**< Captured: onset switch width in n */
    double capt_wc = 1.0; /**< Captured: closure blend width in n */
    double capt_b = 0.0;  /**< Captured: blend start after onset [ws] */
    double Cuq = 1.0;     /**< Scales the uq term in the shear lag equation */
    double Dlr = 0.9;     /**< Wall/wake dissipation length ratio */
    double SlagK = 5.6;   /**< Shear lag constant */

    /* Initial Ctau after transition */
    double CtauC = 1.8; /**< Ctau constant */
    double CtauE = 3.3; /**< Ctau exponent */

    /* G-beta locus: G = GA*sqrt(1+GB*beta) + GC/(H*Rt*sqrt(cf/2)) */
    double GA = 6.7;  /**< G-beta A constant */
    double GB = 0.75; /**< G-beta B constant */
    double GC = 18.0; /**< G-beta C constant */

    /* Operating conditions and thermodynamics */
    double Minf = 0.0;   /**< Freestream Mach number */
    double Vinf = 0.0;   /**< Freestream speed */
    double muinf = 0.0;  /**< Freestream dynamic viscosity */
    double mu0 = 0.0;    /**< Stagnation dynamic viscosity */
    double rho0 = 1.0;   /**< Stagnation density */
    double H0 = 0.0;     /**< Stagnation enthalpy */
    double Tsrat = 0.35; /**< Sutherland Ts/Tref */
    double gam = 1.4;    /**< Ratio of specific heats */
    double KTb = 1.0;    /**< Karman-Tsien beta = sqrt(1-Minf^2) */
    double KTl = 0.0;    /**< Karman-Tsien lambda = Minf^2/(1+KTb)^2 */
    double cps = 0.0;    /**< Sonic cp */

    /* Station information */
    bool simi = false; /**< True at a similarity station */
    bool turb = false; /**< True at a turbulent station */
    bool wake = false; /**< True at a wake station */
};

}  // namespace scirocco
