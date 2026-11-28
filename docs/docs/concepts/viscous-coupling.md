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
   down each surface. With [`warm_start`](#warm-starts) on, a second
   `Solve()` on the same object starts from the previous solution instead.
   It keeps that solution's stagnation panel, and is shifted along each
   surface only if the predicted edge velocity moves the stagnation point
   to another panel.
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

Between Newton steps the boundary layer is re-solved by a march down each
surface, as in XFOIL. Each station is solved with the station upstream fixed
and its edge velocity held close to the current iterate: a combination of
$u_e$ and $H_k$ chosen so that the station is solved in direct mode in
attached flow and in inverse mode near separation, where $u_e$ alone would
make the problem singular. On laminar stretches the march places transition
where the amplification reaches $n_{crit}$, with the same test used to accept
the converged solution. Newton then only has to settle the coupling between
the boundary layer and the outer flow, and typical cases converge in 7-11
steps. The march doesn't change the converged answer: its equations are the
global ones, so a converged state is a fixed point.

Newton steps are damped with
[pseudo-transient continuation](./ptc) when the residual grows. If a cold
start still fails, Scirocco retries along a continuation path in angle of
attack. `Solve()` returns `true` only when the residual meets the
tolerance, or stops falling below `rtol_floor`.

## Warm starts

By default every `Solve()` starts cold, so the same inputs give the same
result whatever the object solved before. Set `param.warm_start = true` to
start each viscous solve from the last converged boundary layer of the same
object, as XFOIL does. A polar step of 0.5° then mostly takes 5 or 6 Newton
iterations instead of about 10. If the warm attempt hasn't converged within
30 iterations, the solve starts again cold.

Near stall the coupled equations can have more than one converged solution,
and which one a warm start reaches depends on the earlier points. On NACA
0012 at Re = 2×10⁵ with 160 panels, sweeping in 0.25° steps:

| α | cold | warm, sweeping up | warm, sweeping down |
| --- | --- | --- | --- |
| 14.00° | 0.746 | 0.983 | 0.746 |
| 14.25° | 0.710 | 0.887 | 0.710 |

Sweeping up, the warm start stays on the attached solution until it ends
between 14.25° and 14.5°. Sweeping down, it stays on the stalled one until
13.75°. Each of these values is a converged solution, and XFOIL shows the
same hysteresis. In this band the solution a warm start reaches is sensitive
enough that it can change with the compiler or the BLAS library. Away from stall a warm and a cold start agree to round-off. Use
warm starts for polars where speed matters, and keep cold starts where a
result must not depend on the order of evaluation, such as in optimization.

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
[`bl_march.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/bl_march.cpp),
[`forces.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/forces.cpp).
