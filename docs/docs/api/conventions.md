---
sidebar_position: 2
---

# Input and output conventions

Use this reference when preparing geometry, interpreting results, or connecting
an upstream geometry component. The [C++ API](./scirocco) and
[Philote server](../optimization/server) share solver behavior but accept
geometry at different stages.

## Geometry layout and panel counts

| Interface | Geometry input | What happens next |
| --- | --- | --- |
| C++ coordinate constructor | Finite `2 × M` source coordinates, at least four points | Re-panels into `npanel + 1` nodes |
| Philote `coordinates` option | Dictionary with equal-length `x` and `z` lists | Initializes and re-panels the seed geometry |
| Philote `x` input | Exactly `2 * (npanel + 1)` values | Replaces panel nodes directly; does not re-panel |
| Kulfan output | `npanel + 1` cosine-spaced nodes, flattened as below | Connects directly to the analysis `x` input |

Coordinates run from the **lower trailing edge, around the leading edge, to
the upper trailing edge**. Matrix row 0 contains x and row 1 contains z.
The remote array is `[x_0, ..., x_N, z_0, ..., z_N]`, not interleaved pairs.
Check the order when importing a `.dat` file; file conventions vary.

Panel counts must be integers from 4 through 4096. Kulfan additionally requires
an even count and 3–256 coefficients per surface. Source coordinates must have
positive chord and no repeated adjacent points. Do not duplicate the leading
edge, reverse the contour, or supply a self-intersecting shape. Input validation
does not replace checking the geometry visually. Upper and lower trailing-edge
points may have a finite gap.

For the remote interface, cluster nodes near the leading and trailing edges
before sending them. Use the [custom-coordinate tutorial](../tutorials/custom-airfoil)
for a complete example. Set exactly one of `naca` or `coordinates` to initialize
the analysis discipline; its connected `x` input supplies the evaluated shape.

## Normalization and references

Use **unit-chord geometry**, with the leading edge near x = 0 and the trailing
edge near x = 1, for the tutorials and conventional coefficient comparisons.
Philote labels `x` with length units, but it does not normalize an arbitrary
physical airfoil automatically. Reynolds number is based on chord; angle of
attack and its derivatives use degrees.

The native state holds reference chord in `geom.chord`, moment reference in
`geom.xref` (default `(0.25, 0)`), and wake extent in `geom.wakelen` (default one
chord). Pressure-integrated lift and near-field pressure drag use reference
chord; moment uses its square. Total drag retains the mfoil/XFOIL unit-reference
convention. Scaling a shape alone is therefore not a general way to obtain
consistently renormalized coefficients.

Reference quantities are held fixed when differentiating coordinates. Moving
an endpoint does not redefine the reference chord or moment reference. See
[Sensitivity validity](./sensitivity#geometry-and-validity) for arc-length modes
and the distinction between panel-node and source-geometry derivatives.

## Operating conditions

| Input | Contract |
| --- | --- |
| `alpha` | Finite angle in degrees |
| `cl` (C++ only) | Finite negative value selects prescribed angle; non-negative value requests a lift-target solve |
| `Re` | Finite; positive selects viscous analysis, nonpositive selects inviscid analysis |
| `Ma` | Finite, `0 <= Ma < 1` |

The remote discipline exposes prescribed-angle analysis. To impose a target
lift in OpenMDAO, make `alpha` a design variable and constrain `cl`.

The accepted Mach range is an input constraint, not a validity guarantee.
Scirocco does not resolve shocks or model three-dimensional flow. Review
[Compressibility](../concepts/compressibility) before interpreting high-Mach
results; difficult separated flows can also fail to converge.

## Results and derivatives

`N` denotes airfoil nodes and `N_w` wake nodes:

| Result | Inviscid | Viscous |
| --- | --- | --- |
| `GetCp()` | `N` entries | `N` entries |
| `GetCf()`, `GetTheta()`, `GetDeltaStar()` | Empty | `N + N_w`, airfoil first |
| `GetWakeCoords()` | Empty | `2 × N_w` |
| `GetCd()` | Zero | Far-wake drag estimate |
| `GetCdp()` | Zero | Total drag minus friction drag |

Near-field pressure drag, `state.post.cdpi`, is a separate quantity from total
`cd` and from the viscous remainder `cdp`.

Check `Solve()` before consuming outputs. Compute sensitivities only after a
successful solve; cached gradients are invalidated by a new solve, changed
operating conditions or geometry, and mutable `GetState()` access. See
[Troubleshooting](../getting-started/troubleshooting) for failure handling.
