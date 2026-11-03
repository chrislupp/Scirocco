---
sidebar_position: 1
---

# Discrete adjoint

Scirocco differentiates the discrete solver equations. This page gives the
adjoint relationship; [Sensitivities](../api/sensitivity) documents the API,
geometry modes, and lift-target derivatives.

## Residual and objective

Let $U$ be the solved state, $p$ the inputs, and $R(U,p)=0$ the discrete
residual at convergence. For one scalar output $J(U,p)$, differentiating gives

$$
R_U\frac{dU}{dp} + R_p = 0,
\qquad
\frac{dJ}{dp} = J_p + J_U\frac{dU}{dp}.
$$

Rather than solving for a state derivative for every input, solve the transpose
system for the adjoint $\lambda$:

$$
R_U^T\lambda = J_U^T.
$$

Then the total derivative is

$$
\frac{dJ}{dp} = J_p - \lambda^T R_p.
$$

One scalar output requires one adjoint right-hand side. The number of these
solves does not grow with the number of design variables, although assembly,
factorization, and derivative contractions still depend on problem size.

## Applying the relationship

The viscous adjoint uses the unregularized converged Jacobian. Damping used to
reach the solution is not included as an extra physical term in that residual.
Derivative chains account for geometry, wake, boundary-layer, and force terms.
The [sensitivity reference](../api/sensitivity#geometry-and-validity) explains
which reference quantities are held fixed and how arc lengths are differentiated.

For a geometry parameterization $x=x(A)$, OpenMDAO combines the panel-coordinate
derivative with the upstream geometry Jacobian:

$$
\frac{dJ}{dA} = \frac{dJ}{dx}\frac{dx}{dA}.
$$

The remote discipline supplies prescribed-angle derivatives. Native lift-target
analysis additionally differentiates the angle adjustment that enforces the
lift constraint. A singular lift/angle slope makes that trim derivative fail.

## Validation and limits

The repository compares derivative components and complete gradients against
mfoil complex-step reference data. The Python tutorial additionally checks the
remote geometry-to-analysis chain with finite differences.

These derivatives describe the implemented discretization at the converged
state. They do not establish the physical accuracy of the underlying flow model,
and branch changes such as transition moving between panels can complicate
numerical checks. See [Troubleshooting](../getting-started/troubleshooting).

Source: [`adjoint.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/adjoint.cpp).
