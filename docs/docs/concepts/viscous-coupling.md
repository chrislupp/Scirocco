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
   down each surface. With `warm_start` (the default), a second `Solve()` on
   the same object starts from the previous solution instead, shifted along
   each surface to the new stagnation point.
3. Run Newton iterations on the full coupled system, updating the
   stagnation point and the transition location as they move.

The Jacobian is exact, so a well-posed case converges in a handful of steps:

<ConvergencePlot variant="typical" />

Each Newton step solves a linear system with $4N$ unknowns. The coupling
rows are dense: every node's $u_e$ depends on the mass defect of every other
node. Substituting $m$ for $\delta^*$ turns that block into a diagonal, and
eliminating the boundary-layer unknowns, which couple only to the station
upstream, leaves an $N \times N$ dense system. The step costs about 2 ms at
200 panels and is exact; the adjoint uses the same factorization.

Far from the solution, a Newton step can move the boundary layer somewhere
its own transition criterion disagrees with, and the next update then throws
transition far upstream or downstream. While the residual is above
`march_residual`, each update is therefore followed by a coupled march, as
in XFOIL: every station is re-solved together with its own edge-velocity
equation, holding the other nodes' mass defect fixed, and transition is
placed where the amplification reaches $n_{crit}$. The marched state is kept
only if it lowers the residual. It doesn't change the converged answer;
transition is checked once more before the solve is accepted.

Newton steps are damped with
[pseudo-transient continuation](./ptc) when the residual grows. If a cold
start still fails, Scirocco retries along a continuation path in angle of
attack. `Solve()` returns `true` only when the residual meets the
tolerance, or stops falling below `rtol_floor`.

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
[`coupled_linear.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/coupled_linear.cpp),
[`coupled_march.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/coupled_march.cpp),
[`forces.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/forces.cpp).
