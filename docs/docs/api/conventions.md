---
title: Conventions
description: Node order, units, valid input ranges, reference quantities, and array sizes.
---

# Conventions

## Node order and layout

Nodes run from the **lower trailing edge**, around the leading edge, to the
**upper trailing edge**. In C++, coordinates are a `2 × N` matrix: row 0 is
x, row 1 is z. On the Philote server, `x` is flat: `[x_0 … x_N, z_0 … z_N]`.

| Interface | Geometry input | What happens to it |
| --- | --- | --- |
| C++ coordinate constructor | `2 × M` points, at least 4 | Re-paneled to `npanel + 1` nodes |
| Philote `coordinates` option | `{"x": [...], "z": [...]}` | Re-paneled; only initializes the server |
| Philote `x` input | Exactly `2 * (npanel + 1)` values | Used as the panel nodes, not re-paneled |
| CST server `x` output | `npanel + 1` cosine-spaced nodes | Connects straight to the analysis `x` input |

Many `.dat` files list points from the upper trailing edge first; reverse
them.

## Valid geometry

- Panel count from 4 to 4096. The CST server also needs an even count and
  3–256 coefficients per surface.
- Finite values, positive chord, no repeated neighbouring points, no
  self-intersection, and the leading edge listed once.
- A finite trailing-edge gap is fine, and recommended: viscous solves on an
  exactly closed trailing edge can fail (a known issue).
- For the Philote `x` input, cluster nodes at both edges yourself, since the
  server doesn't re-panel.

## Units and references

| Quantity | Convention |
| --- | --- |
| Coordinates | Fractions of chord. Leading edge near `x = 0`, trailing edge near `x = 1`. |
| `alpha` | Degrees, as are all `d*/dalpha` derivatives |
| `Re` | Based on chord |
| Moment reference | `(0.25, 0)`, in `geom.xref` |
| Reference chord | `geom.chord` |
| Wake length | One chord, in `geom.wakelen` |

The Philote `x` variable is labelled with length units, but values are
still read as fractions of chord. Nothing rescales them.

`Cl` and `cdpi` are normalized by the reference chord and `Cm` by its
square. `Cd` keeps mfoil's convention of a unit reference length, so to
compare coefficients across shapes use unit-chord geometry rather than
scaling a section.

## Operating point

| Input | Valid values |
| --- | --- |
| `alpha` | Any finite angle |
| `cl` (C++ only) | Negative: fixed angle. Zero or positive: lift target. |
| `Re` | Finite. Positive: viscous. Zero or negative: inviscid. |
| `Ma` | `0 <= Ma < 1`. The model is valid only while the surface stays subsonic. |

The Philote server always runs at a fixed angle. For a lift target, make
`alpha` a design variable and constrain `cl`.

## Array sizes

`N` is `npanel + 1` airfoil nodes; `Nw` is the number of wake nodes.

| Result | Inviscid | Viscous |
| --- | --- | --- |
| `GetCp()` | `N` | `N` |
| `GetCf()`, `GetTheta()`, `GetDeltaStar()` | empty | `N + Nw`, airfoil first |
| `GetWakeCoords()` | empty | `2 × Nw` |
| `GetCd()` | `0` | Wake momentum deficit |
| `GetCdp()` | `0` | `Cd − Cdf` |
| `d*_dx` gradients | `2 × N` | `2 × N` |
