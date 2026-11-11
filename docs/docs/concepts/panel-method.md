---
title: Panel method
description: The linear-vortex panel method behind every Scirocco solve.
---

# Panel method

Every solve starts with potential flow around the airfoil. Scirocco covers
the surface with straight panels and places a vortex sheet on them whose
strength varies linearly from node to node.

<PanelPlot />

## The equations

An airfoil with `npanel` panels has $N = \text{npanel} + 1$ nodes, numbered
from the lower trailing edge, around the leading edge, to the upper trailing
edge. The unknowns are the $N$ nodal vortex strengths $\gamma_i$ and one
constant, the streamfunction value $\Psi_0$ on the surface.

Flow can't cross the surface, so the streamfunction must equal $\Psi_0$ at
every node. That gives $N$ equations. The **Kutta condition** closes the
system by making the flow leave the trailing edge smoothly:

$$
\gamma_0 + \gamma_{N-1} = 0
$$

A finite trailing-edge gap is modelled with an extra source and vortex panel
across the gap. Scirocco factors the dense $(N+1) \times (N+1)$ system once
and solves it for freestream at 0° and 90°. Any angle of attack is then a
weighted sum of the two, which is why an inviscid solve takes about a
millisecond once the airfoil is paneled.

## From vortex strength to forces

On the surface, the edge velocity equals the vortex strength, so the
pressure coefficient is

$$
C_p = 1 - \left(\frac{u_e}{V_\infty}\right)^2 .
$$

At nonzero Mach number, the [Kármán–Tsien correction](./compressibility) is
applied to this. Integrating pressure around the surface gives $C_l$ and
$C_m$; the moment is taken about the quarter chord, $(0.25, 0)$.

In potential flow, the drag of a 2D body is zero (d'Alembert's paradox), so
inviscid `Cd` is reported as exactly zero. A separate near-field pressure
drag, `cdpi`, is available in `GetState().post` and contains only
discretization error.

## Limits

Potential flow knows nothing about the boundary layer, so it can't predict
drag, separation, or stall. An inviscid lift curve stays straight long after
the real airfoil has stalled. Use a [viscous solve](./viscous-coupling) for
anything that depends on the boundary layer.

Source: [`inviscid.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/inviscid.cpp),
[`panel_influence.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/panel_influence.cpp).
