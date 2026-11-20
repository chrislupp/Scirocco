---
title: The adjoint
description: How Scirocco differentiates its discrete equations, and how the result is tested.
---

# The adjoint

Scirocco's gradients come from a discrete adjoint: it differentiates the
equations the solver actually solves, so the derivatives are exact for the
discretized model, to the solver's tolerance.

## The idea

Let $U$ be the converged state, $p$ the inputs (coordinates, $\alpha$, $Re$,
$M_\infty$), and $R(U, p) = 0$ the discrete equations. For an output
$J(U, p)$ such as $C_d$, the chain rule gives

$$
\frac{dJ}{dp} = J_p + J_U \frac{dU}{dp},
\qquad
R_U \frac{dU}{dp} = -R_p .
$$

Computing $dU/dp$ directly costs one linear solve per input. Instead, solve
one transposed system for the adjoint vector $\lambda$:

$$
R_U^T \lambda = J_U^T ,
\qquad
\frac{dJ}{dp} = J_p - \lambda^T R_p .
$$

That's one solve per **output**, however many inputs there are. Scirocco
reuses the factorization of $R_U$ from the last Newton step, so each
adjoint solve is a back-substitution. Most of the cost is in the product
$\lambda^T R_p$: differentiating the residual, and above all the panel
influence coefficients, with respect to every node coordinate. See
[Performance](../validation/performance).

## The chain

Moving a node changes more than the panel it sits on. The derivative chain
runs through:

1. The panel influence coefficients and the inviscid solution.
2. Arc lengths along the surface and the wake, which follow the nodes.
3. The boundary-layer residuals at every station.
4. The force integration.

With a lift target, the chain also includes the angle-of-attack change that
keeps $C_l$ on target. The adjoint always uses the plain converged
Jacobian, without the [PTC](./ptc) damping term.

## How it's tested

The reference is mfoil, the Python code Scirocco was ported from, run with
complex-step differentiation. Complex step gives derivatives accurate to
machine precision, with no step-size trade-off. The test suites compare
against it piece by piece (influence coefficients, residual partials, the
geometric chain, force partials) and end to end, for coordinates, $\alpha$,
$Re$, and $M_\infty$.

## Limits

The gradient describes the discretized model at the converged state. Where
transition jumps from one panel to the next, the discrete model has a kink,
and a gradient taken on one side doesn't describe the other. The gradient
says nothing about whether the model itself is physically accurate.

Source: [`adjoint.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/adjoint.cpp).
