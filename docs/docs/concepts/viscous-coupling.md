---
title: Viscous coupling
description: How the boundary layer, wake, and potential flow are solved together.
---

# Viscous coupling

With `Re > 0`, Scirocco adds a boundary layer on both surfaces and a wake
behind the trailing edge. The boundary layer displaces the outer flow, and
the outer flow drives the boundary layer, so both are solved together as one
nonlinear system with Newton's method.

<CpPlot />

## The state

Each station on the airfoil and in the wake carries four unknowns:

| Unknown | Meaning |
| --- | --- |
| $\theta$ | Momentum thickness |
| $\delta^*$ | Displacement thickness |
| $\tilde n$ or $C_\tau^{1/2}$ | Amplification factor where the flow is laminar, shear-stress coefficient where it is turbulent |
| $u_e$ | Edge velocity |

The residual combines the integral momentum and shape-parameter equations,
an amplification or shear-lag equation, and the coupling to the outer flow.
The outer flow sees the boundary layer through the mass defect
$m = u_e \delta^*$: a precomputed influence matrix maps changes in $m$ to
changes in $u_e$.

## The solve

1. Solve potential flow and lay out the wake along a streamline.
2. Find the stagnation point and march a first guess for the boundary layer
   down each surface.
3. Run Newton iterations on the full coupled system, updating the
   stagnation point and the transition location as they move.

The Jacobian is exact, so a well-posed case converges in a handful of steps:

<ConvergencePlot variant="typical" />

Newton steps are damped with
[pseudo-transient continuation](./ptc) when the residual grows. If a cold
start still fails, Scirocco retries along a continuation path in angle of
attack. `Solve()` returns `true` only when the residual meets the
tolerance.

## Outputs

| Quantity | Where it comes from |
| --- | --- |
| `cl`, `cm` | Surface pressure integration |
| `cd` | Momentum deficit at the end of the wake, extrapolated to infinity (Squire–Young) |
| `cdf` | Integrated skin friction |
| `cdp` | `cd − cdf` |

## Limits

Integral boundary-layer models are reliable for attached flow and small
separation bubbles. Near and past stall, the solve may still converge, but
the answer is less trustworthy. Converging is not the same as being right;
compare against data for your regime.

Source: [`boundary_layer.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/boundary_layer.cpp),
[`coupled_solver.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/coupled_solver.cpp),
[`forces.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/forces.cpp).
