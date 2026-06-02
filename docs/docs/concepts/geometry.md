---
sidebar_position: 5
---

# Geometry

Scirocco accepts geometry in three forms — NACA digit strings, raw point
clouds, and a parametric CST/Kulfan basis — and exposes a small set of
modifier operations (flap, camber, derotation) that act on whichever form
you started with.

This page surveys the **conceptual** picture. For the hands-on usage see
the tutorials [Custom airfoil](../tutorials/custom-airfoil) and
[Geometry modifications](../tutorials/geometry-mods); for signatures see
[API → Geometry](../api/geometry).

## NACA 4- and 5-digit

`NacaPoints(geom, digits)` in
[`core/geometry.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/geometry.h)
generates analytical NACA points and writes them into the
`Geom` substruct. The two families:

- **4-digit** (e.g. `"2412"`): camber, position of max camber, thickness —
  the classic Jacobs formulas.
- **5-digit** (e.g. `"23012"`): adds a design lift coefficient and a
  reflex bit, otherwise the same thickness distribution.

Coordinates are evaluated at cosine-clustered $x$-stations so the leading
edge is naturally well resolved.

## Raw point clouds

`SetCoords(geom, coords)` ingests a `2 × N` matrix, TE-lower → LE →
TE-upper. Scirocco does not require regularly-spaced points: the cubic-
spline paneler `MakePanels` rebuilds the spacing on its own. Bring an
XFOIL `.dat` file or a CAD export and you're set.

The internal `Geom` struct keeps both the **input** point cloud
(`geom.xpoint`) and the **paneled** representation (`foil.x`) separate.
After any modification — `GeomFlap`, `GeomAddCamber`, `GeomDerotate` — the
input cloud is mutated and the paneled representation is rebuilt.

## CST / Kulfan parameterization

Class–shape transformation (Kulfan, 2008) writes each surface as

$$
y(\psi) = C(\psi)\, S(\psi) + \psi \cdot y_\text{TE},
\qquad C(\psi) = \psi^{0.5}(1-\psi)^{1},
\qquad S(\psi) = \sum_{i=0}^{n} A_i \binom{n}{i} \psi^i (1-\psi)^{n-i}
$$

with $\psi = x/c \in [0,1]$. The class function $C(\psi)$ gives a round
leading edge and a sharp trailing edge; the shape function $S(\psi)$ is a
Bernstein polynomial whose coefficients $A_i$ become the optimization
design variables.

A useful property of Bernstein bases: the coefficients are non-negative
multipliers of non-negative basis functions, so requiring $A_i \ge 0$ for
all $i$ on the upper surface and $A_i \ge 0$ on the lower surface (with
opposite sign convention) is **sufficient** for the airfoil to be
non-self-intersecting at every iterate — a cheap and powerful
constraint.

Scirocco implements CST as a Philote-Cpp discipline,
[`KulfanDiscipline`](../api/kulfan-discipline), which takes the per-surface
$A_i$ vectors as inputs and emits `x = (x, z)` panel-node coordinates as
output. You feed those coordinates straight into the analysis discipline.

### Smoothness constraint

The discipline also emits second-difference quantities
$D_i = A_{i+2} - 2 A_{i+1} + A_i$ as the `smoothness_*` outputs. Because

$$
S''(\psi) = n(n-1) \sum_i D_i\, B^{n-2}_i(\psi),
$$

bounding $|D_i|$ caps the curvature of $S$ (and hence the airfoil) with a
linear constraint. This is the right knob to stabilize a CST-driven
viscous optimizer that wants to fluctuate the basis coefficients — see
[Optimization → CST design](../optimization/cst-design).

## Modifier operations

All three modifiers in
[`core/geometry.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/geometry.h)
operate on the current `geom.xpoint` and re-panel afterwards. They compose
in call order:

| Function | Effect |
| --- | --- |
| `GeomFlap(state, hinge, eta)` | Rotates everything aft of `hinge` by `eta` degrees |
| `GeomAddCamber(state, xzcamb)` | Linearly interpolates a $(x, z)$ camber increment onto the chord; keeps thickness distribution |
| `GeomDerotate(state)` | Aligns the LE–TE line with the $x$-axis |

For everyday use, the wrappers on the `Scirocco` class (`foil.GeomFlap(...)`,
etc.) are equivalent and more convenient.
