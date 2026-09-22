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

/* Parameters */
struct Param {
    /* Solver parameters */
    int verb = 0; /**< Printing verbosity level (higher = more verbose) */
    double rtol = 1e-11; /**< Residual tolerance for Newton */
    bool scale_coupled_system = true; /**< Equilibrate equations and variables */
    bool use_trust_region = false; /**< Experimental scaled dogleg step control */
    bool trust_update_topology = true; /**< Evaluate the complete trial state */
    double trust_radius_initial = 10.0;
    int trust_max_trials = 16;
    bool use_continuation = true; /**< Retry failed cold starts along an angle/Cl path */
    int continuation_max_attempts = 24;
    int continuation_max_iterations = 800; /**< Total extra coupled iterations */
    int niglob = 200;    /**< Maximum number of global iterations (raised from
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

    /* Viscous parameters */
    double ncrit = 9.0; /**< Critical amplification factor */
    double Cuq = 1.0;   /**< Scales the uq term in the shear lag equation */
    double Dlr = 0.9;   /**< Wall/wake dissipation length ratio */
    double SlagK = 5.6; /**< Shear lag constant */

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
