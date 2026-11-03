---
sidebar_position: 1
---

# Panel method

Scirocco represents two-dimensional potential flow with a linear vortex sheet
on straight panels. It provides the inviscid solution and the outer-flow part
of the viscous analysis.

## Geometry and unknowns

An airfoil with `npanel` panels has `N = npanel + 1` nodes, ordered lower
trailing edge → leading edge → upper trailing edge. Vortex strength varies
linearly between nodes. The C++ constructors use a curvature-based paneler;
the remote `x` input replaces panel nodes directly. See
[Input and output conventions](../api/conventions).

`BuildGamma` assembles the streamfunction influence matrix. The unknowns are
`N` nodal vortex strengths and one common surface streamfunction. Requiring a
constant streamfunction along the airfoil enforces flow tangency. The trailing
edge adds the Kutta condition:

$$
\gamma_0 + \gamma_{N-1} = 0.
$$

The implementation includes trailing-edge source/vortex contributions and a
special extrapolation equation for a closed trailing edge. It solves the dense
system for two freestream directions and combines them for the requested angle.
See [`inviscid.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/inviscid.cpp)
for the assembly and treatment of the trailing edge.

## Pressure and forces

In incompressible flow, surface edge velocity gives

$$
C_p = 1 - (u_e / V_\infty)^2.
$$

At nonzero Mach number, Scirocco applies the
[Kármán–Tsien correction](./compressibility). Pressure integration provides
lift and moment; the default moment reference is `(0.25, 0)` for unit chord.

Inviscid total `Cd` and its derivatives are zero by convention. The separate
near-field pressure-drag quantity `cdpi` can contain discretization effects.
An inviscid run has no boundary-layer distributions or wake coordinates.
For profile drag and wake results, use [viscous analysis](./viscous-coupling).

## Practical limits

The inviscid model cannot predict viscous separation or stall. A nearly linear
lift curve at large angles is not evidence that the real airfoil remains
attached. Check panel-count sensitivity and use a model appropriate to the
flow regime before interpreting results.
