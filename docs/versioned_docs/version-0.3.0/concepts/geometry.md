---
title: Geometry and CST
description: How Scirocco represents shapes, panels them, and parameterizes them with CST.
---

# Geometry and CST

Scirocco keeps two versions of the shape. The **source points** are what you
gave it: a NACA section or your coordinates (`geom.xpoint`). The **panel
nodes** are what it solves on (`foil.x`), laid out by a spline-based paneler
that clusters nodes where curvature is high.

<PanelPlot />

Flaps, camber changes, and derotation edit the source points and re-panel.
The Philote server works differently: its `x` input replaces the panel
nodes directly, without re-paneling, so a geometry component upstream
controls the layout exactly.

## NACA sections

`"2412"` and `"23012"` generate 4- and 5-digit sections from the classic
formulas. The generator produces 201 source points; the paneler then
redistributes them to the panel count you ask for. Five-digit sections must
be the non-reflexed `2X0YY` family.

## CST (Kulfan) parameterization

For optimization, the CST parameterization describes each surface with a
handful of coefficients. With $\psi = x/c$:

$$
z(\psi) = C(\psi)\,S(\psi) + \psi\, z_\text{TE},
\qquad C(\psi) = \sqrt{\psi}\,(1 - \psi),
\qquad S(\psi) = \sum_{i=0}^{n} A_i \binom{n}{i} \psi^i (1-\psi)^{n-i}
$$

The class function $C$ gives a round leading edge and a pointed trailing
edge. The shape function $S$ is a Bernstein polynomial whose coefficients
$A_i$ are the design variables. Scirocco subtracts the lower surface, so
**non-negative coefficients on both surfaces keep the upper surface above
the lower one**. They don't stop the section getting arbitrarily thin, so
optimizations still need a thickness or area constraint.

CST runs as its own Philote discipline,
[`KulfanDiscipline`](../api/kulfan-discipline). It takes `cst_upper` and
`cst_lower` and returns panel coordinates, with exact partials, ready for
the analysis server.

### Smoothness outputs

The discipline also returns the second differences of each coefficient
vector,

$$
D_i = A_{i+2} - 2A_{i+1} + A_i ,
$$

as `smoothness_upper` and `smoothness_lower`. Since
$S''(\psi) = n(n-1)\sum_i D_i\, B^{n-2}_i(\psi)$, bounding $|D_i|$ bounds the
curvature of the shape function and stops the optimizer from building
wiggles. It doesn't bound the surface curvature directly, because the class
function contributes too. See [CST design](../optimization/cst-design).

Source: [`geometry.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/geometry.cpp),
[`paneling.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/paneling.cpp).
