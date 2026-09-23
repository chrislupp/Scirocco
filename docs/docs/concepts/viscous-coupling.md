---
draft: true
sidebar_position: 2
---

# Viscous coupling

When `Re > 0`, Scirocco extends the inviscid solve with an **integral
boundary-layer (BL) system** marched along each side of the airfoil and
along the wake, and couples it back to the inviscid edge velocity through
displacement effects. The whole thing is solved as a single nonlinear
system with Newton's method. Implementation in
[`core/boundary_layer.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/boundary_layer.h)
and
[`core/coupled_solver.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/coupled_solver.h).

## State variables

At every BL station the unknowns are

| Symbol | Meaning |
| --- | --- |
| $\theta$ | Momentum thickness |
| $\delta^*$ | Displacement thickness |
| $\tilde n$ (laminar) / $C_\tau$ (turbulent) | Amplification factor / shear-stress coefficient |
| $u_e$ | Edge velocity (also an inviscid unknown — this is the coupling slot) |

These are stacked into a global state vector $\mathbf{U}$ together with the
nodal vortex strengths $\gamma$ from the inviscid system.

## Governing equations

The integral BL equations are the **momentum** and **shape-parameter**
equations:

$$
\frac{d\theta}{ds} + (2 + H - M_e^2)\,\frac{\theta}{u_e}\,\frac{du_e}{ds} = \frac{C_f}{2}
$$

$$
\theta \frac{dH^*}{ds} + (2H^{**} + H^*(1-H))\,\frac{\theta}{u_e}\,\frac{du_e}{ds} = 2 C_D - H^* C_f / 2
$$

In laminar regions an amplification-factor transport equation tracks
$\tilde n(s)$ for transition prediction (see
[Transition](./transition)). In turbulent regions a shear-lag equation for
$C_\tau$ replaces it, using the G-$\beta$ closure
(constants `GA`, `GB`, `GC`, `SlagK`, `Cuq`, `Dlr` in
[`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h)).

Closure relations $H^* = H^*(H, Re_\theta)$, $C_f = C_f(H, Re_\theta)$,
$C_D = C_D(H, Re_\theta)$ come from Drela/Giles fits (mfoil's
correlations, ported in `bl_properties.h`).

## Coupling slot

Inviscid panel methods solve for an edge velocity $u_e^{\text{inv}}$
assuming the wall is the body surface. A boundary layer thickens the
effective body by $\delta^*(s)$, which perturbs $u_e$. Scirocco does this
not by physically moving the wall, but with the **mass-source equivalent**:
the boundary-layer displacement effect is represented by a distribution of
sources $m(s) = u_e\,\delta^*$ on the original surface, and the inviscid
system gains a sensitivity matrix

$$
u_e = u_e^{\text{inv}}(\gamma, \alpha) + \mathbf{U}_{e,m}\, \mathbf{m}.
$$

`CalcUeM` and `RebuildUeM` in
[`core/coupled_solver.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/coupled_solver.h)
build and maintain $\mathbf{U}_{e,m}$ each iteration. The result is that
$u_e$ appears in both the inviscid Kutta system and the BL closure
equations — that's the coupling.

## The Newton iteration

`SolveCoupled` builds the residual $\mathbf{R}(\mathbf{U})$ for the joint
system and steps

$$
\mathbf{U}^{k+1} = \mathbf{U}^k - \alpha_k\, \mathbf{R}_U^{-1}(\mathbf{U}^k)\, \mathbf{R}(\mathbf{U}^k)
$$

`alpha_k` is an under-relaxation factor chosen by the line search in
`UpdateState`. The Jacobian $\mathbf{R}_U$ is sparse — block-banded along
each surface march, plus the dense $\mathbf{U}_{e,m}$ coupling — and is
factored with sparse LU. The factorization is kept around after
convergence for the adjoint solve (see
[Sensitivities tutorial](../tutorials/sensitivities) and
`StoreAdjointFactorization`).

If the residual stops decreasing, **pseudo-transient continuation** (PTC)
adds a $1/\tau$ ridge to the diagonal of $\mathbf{R}_U$ — same converged
solution, much steadier path to it. See
[PTC](./ptc).

## Stagnation, transition, and wake march

Each `SolveCoupled` iteration:

1. Finds the stagnation point (the panel node where inviscid $u_e$ changes
   sign), splits the airfoil into upper and lower surface arrays
2. Initializes laminar BL with Thwaites from the stagnation point
3. Marches $\tilde n$ along each laminar segment; triggers transition where
   $\tilde n$ crosses `param.ncrit` (default 9; see
   [Transition](./transition))
4. Continues with the turbulent + shear-lag closure to the trailing edge
5. Marches both surfaces through the wake, using a wake-specific closure
   (no $C_f$, modified $H^*$)
6. Recomputes $u_e$ via the coupling sensitivity and updates the residual

## Outputs

After convergence, force coefficients are integrated:

| Quantity | Source |
| --- | --- |
| $C_l$, $C_m$ | Pressure integration (as in the inviscid case) |
| $C_{d,p}$ | Pressure integration, **with** displacement effect on $u_e$ |
| $C_{d,f}$ | $\int C_f\, \hat{\mathbf{t}} \cdot \hat{\mathbf{x}}\, ds$ |
| $C_d$ | Far-wake momentum deficit (Squire–Young extrapolation) |

The far-wake estimate is the textbook "best" 2D drag number: it includes
parasitic, pressure, and the wake-momentum closing terms in one figure.
