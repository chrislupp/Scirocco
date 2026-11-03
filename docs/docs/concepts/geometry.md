---
sidebar_position: 5
---

# Geometry

Scirocco accepts geometry in three forms — NACA digit strings, ordered source
coordinates, and a parametric CST/Kulfan basis — and exposes a small set of
modifier operations (flap, camber, derotation) that act on whichever form
you started with. Native modifiers act on source geometry; they are not
inputs on the remote discipline. For a remote CST model, change the coefficients
or supply coordinates from a Python geometry component.

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
- **5-digit** (e.g. `"23012"`): the implemented non-reflex family begins
  with `2X0`, where `X` is 1–5. Other five-digit families are rejected.

The generator produces 201 source points using a nonuniform spacing law.
The constructor then re-panels those points to the requested panel count.

## Raw point clouds

`SetCoords(geom, coords)` ingests a `2 × N` matrix, TE-lower → LE →
TE-upper. Scirocco does not require regularly-spaced points: the cubic-
spline paneler `MakePanels` rebuilds the spacing on its own. Parse imported
files into the required matrix, check contour order, and normalize the chord
as described in [Conventions](../api/conventions).

The internal `FoilState` keeps **input** coordinates in `geom.xpoint`
and **panel-node** coordinates in `foil.x`.
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

Bernstein basis functions are non-negative on the chord. With Scirocco's
opposite upper/lower sign convention, non-negative coefficients and
`te_upper >= te_lower` keep the upper surface at or above the lower surface.
This still allows zero thickness; practical designs need thickness or area
constraints as well.

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

bounding $|D_i|$ bounds the second derivative of the shape polynomial and
discourages coefficient oscillations. This is not a direct surface-curvature
bound, because the class function also contributes derivatives. See
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
