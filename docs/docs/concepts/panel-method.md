---
draft: true
sidebar_position: 1
---

# Panel method

The inviscid backbone of Scirocco is a 2D linear-vortex panel method. This
page describes what `Solve()` does on an inviscid call ($Re = 0$): paneling,
the AIC matrix, the Kutta condition, and how forces are integrated. The
implementation lives in
[`core/inviscid.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/inviscid.h),
[`core/panel_influence.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/panel_influence.h),
and
[`core/paneling.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/paneling.h).

## Discretization

The airfoil contour $\mathbf{x}(s)$ is approximated by $N$ panel nodes,
ordered TE-lower → LE → TE-upper. Between consecutive nodes lies a flat
panel; on each panel, a vortex sheet of strength $\gamma$ varies linearly
between the two end-node values.

Panel placement is **curvature-graded** by `MakePanels`: the cubic spline of
the input geometry is reparameterized by an arc-length-plus-curvature metric,
so the leading edge and other tight curves get extra nodes automatically.
You ask for a node count; the paneler decides where to put them.

## Influence coefficients

The disturbance velocity at any point $\mathbf{x}$ produced by a unit
linear-vortex panel between nodes $j$ and $j+1$ has a closed-form
expression. Summing across all $N$ panels, the tangential velocity at panel
node $i$ is a linear combination of the $N$ unknown nodal vortex strengths
$\gamma_j$, plus a contribution from the freestream:

$$
u_t(\mathbf{x}_i) = \sum_{j=1}^{N} A_{ij}\, \gamma_j
                  + (U_\infty \cos\alpha,\, U_\infty \sin\alpha) \cdot \hat{\mathbf{t}}_i
$$

The matrix $\mathbf{A}$ is the **aerodynamic influence-coefficient (AIC)**
matrix. It depends only on geometry, so for a fixed shape it's built once
per call.

## Boundary condition and Kutta closure

Inviscid flow tangency at every node gives $N$ equations. Linear-vortex
panels share end-node strengths, so the unknowns are also $N$. That's
exactly square — and singular: any constant $\gamma_j \equiv c$ produces
zero net tangential velocity, so the system has a one-dimensional null
space.

The **Kutta condition** closes it: at the trailing edge, the upper and
lower-surface vortex strengths must be equal and opposite,

$$
\gamma_1 + \gamma_N = 0,
$$

which is added as one extra equation and one extra (Lagrange-multiplier-like)
unknown. The resulting square non-singular system is solved with a dense LU
factorization.

In code, this is the work of `BuildGamma` in
[`core/inviscid.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/inviscid.h);
`SolveInviscid` is the top-level wrapper that initializes thermodynamics,
builds the system, solves, and computes edge velocities $u_e$ at every
panel node and wake node.

## Force integration

With $u_e$ known at every node, the local pressure coefficient comes from
Bernoulli with the Kármán–Tsien correction at non-zero Mach (see
[Compressibility](./compressibility)):

$$
C_p = 1 - (u_e / U_\infty)^2 \quad \text{(incompressible form)}
$$

Pressure is integrated around the contour to give force coefficients. The
moment is computed about the quarter-chord. The implementation lives in
[`core/forces.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/forces.h).

## What's *not* in the inviscid solve

A pure inviscid run does **not** march any boundary-layer state, never
solves a wake system, and produces $C_d \equiv 0$ — that's D'Alembert's
paradox showing up in the discrete formulation. To get real drag you need
the viscous coupling described in [Viscous coupling](./viscous-coupling).

## Wake

Even for inviscid runs Scirocco builds a discrete wake centerline
downstream of the trailing edge — it's needed as soon as you turn on
viscosity to march the BL through the wake. For inviscid use the wake is
present but its only role is to support sensible far-field velocity
distributions. See
[`core/wake.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/wake.h).
