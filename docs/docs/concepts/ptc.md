---
title: Pseudo-transient continuation
description: The damping Scirocco adds to Newton steps when the residual grows.
---

# Pseudo-transient continuation

Newton's method converges fast once it is close to the answer, but a bad
step far from it can throw the iteration off course. Pseudo-transient
continuation (PTC) adds damping when the residual grows and removes it as
the residual falls. It's on by default.

## The step

A plain Newton step solves

$$
R_U\, \Delta U = -R .
$$

A PTC step adds a diagonal term,

$$
\left(R_U + \frac{1}{\tau} I\right) \Delta U = -R .
$$

where $I$ puts $1/\tau$ on each equation's own unknown: boundary-layer
equation $k$ of a station on unknown $k$ of that station, and the coupling
equation on $u_e$. Large $\tau$ gives back plain Newton. Small $\tau$
shortens the step towards $-\tau R$, a cautious move down the residual. The converged answer,
$R = 0$, doesn't depend on $\tau$, so PTC changes the path, not the
destination.

## How τ changes

After each step, Scirocco compares the residual with the previous one,
$r = \|R^k\| / \|R^{k-1}\|$:

- If the residual fell ($r < 1$), $\tau$ grows by $1/r$, at most by a
  factor of `tau_grow_clamp`.
- If it didn't, $\tau$ shrinks by $1/r$, at most by a factor of
  `tau_shrink_clamp`.
- $\tau$ always stays between `tau_min` and `tau_max`.

It starts at `tau_init = 1e6`, which is effectively plain Newton. On an
easy case PTC never engages.

## What it does in practice

PTC is a safeguard, not a speed-up. On this difficult start both variants
spend about ten steps near $\|R\| \approx 4$ before Newton takes over, and
finish within a step of each other:

<ConvergencePlot variant="hard" />

Its real value is on starts where an undamped step would diverge, such as
the odd shapes an optimizer tries along the way.

## Tuning

| Symptom | Try |
| --- | --- |
| Diverges in the first one or two steps | Start damped: `param.tau_init = 1e2` |
| Residual oscillates between two levels | Slow the return to Newton: `param.tau_grow_clamp = 2.0` |
| Still falling when it hits the iteration limit | Raise `param.niglob` |
| You want a plain-Newton baseline | `param.use_ptc = false` |

All fields are listed under [Solver parameters](../api/solver-params#pseudo-transient-continuation).

## Gradients aren't affected

The adjoint uses the plain Jacobian $R_U$ at the converged state, without
the $1/\tau$ term. If PTC and plain Newton reach the same solution, they
give the same gradients.
