---
title: Sensitivities
description: Reference for the Sensitivity struct, gradient getters, and what the derivatives hold fixed.
---

# Sensitivities

```cpp
#include <scirocco/scirocco.h>   // Sensitivity is in scirocco/data/sensitivity.h
```

`Scirocco::ComputeSensitivities()` returns a `Sensitivity` with the
derivatives of $C_l$, $C_d$, $C_m$, and near-field pressure drag. For a
walkthrough, see [Compute gradients](../guides/gradients).

## `Sensitivity`

| Field | Shape | Derivative |
| --- | --- | --- |
| `dcl_dx`, `dcd_dx`, `dcm_dx`, `dcdpi_dx` | `2 × N` | With respect to node coordinates: row 0 is x, row 1 is z |
| `dcl_dalpha`, `dcd_dalpha`, `dcm_dalpha`, `dcdpi_dalpha` | scalar | Per degree |
| `dcl_dRe`, `dcd_dRe`, `dcm_dRe`, `dcdpi_dRe` | scalar | Per unit Reynolds number |
| `dcl_dMa`, `dcd_dMa`, `dcm_dMa`, `dcdpi_dMa` | scalar | Per unit Mach number |

In lift-target mode there are also:

| Field | Derivative |
| --- | --- |
| `dcl_dcltarget`, `dcd_dcltarget`, `dcm_dcltarget`, `dcdpi_dcltarget` | Of each output with respect to the target |
| `dalpha_dx`, `dalpha_dRe`, `dalpha_dMa`, `dalpha_dcltarget`, `dalpha_dalpha` | Of the solved angle of attack |

`sens.dcd_dx(1, j)` is $\partial C_d / \partial z_j$. Nodes are in the same
order as `GetFoilCoords()`.

## Requirements

- Call after a successful `Solve()`. Otherwise it throws.
- If the adjoint linear solve fails, it throws. You never get stale, zero,
  or regularized gradients.
- The cached getters (`GetDcdDx()` and so on) throw after anything that
  invalidates the gradient: a new solve, `SetOper`, a geometry change, or a
  call to the non-const `GetState()`.

## What's held fixed

The derivatives describe moving nodes, not redefining the problem:

- Reference chord, moment reference point, and wake length stay fixed. Moving
  the trailing-edge node doesn't change the chord used for the coefficients.
- They are derivatives with respect to **panel nodes**, after paneling. The
  paneler that turned your source points into nodes isn't differentiated.

Surface arc lengths follow the nodes in one of three ways
(`ArcLengthMode`):

| Mode | Set by | Arc lengths |
| --- | --- | --- |
| `Fixed` | The constructors | Held at the values from the original paneling, as in mfoil |
| `Spline` | `RefreshPanelGeometry()` | Re-integrated along a spline through the current nodes. The Philote server uses this. |
| `Cumulative` | `RecomputeArcLength()` | Sums of straight panel lengths |

Each mode is differentiated consistently with how it computes arc length.

```cpp
#include <scirocco/core/paneling.h>

// Nodes are your design variables, as in the Philote server:
scirocco::RefreshPanelGeometry(foil.GetState().foil);
foil.SetOper(alpha, -1.0, Re, Ma);
if (foil.Solve()) {
    auto sens = foil.ComputeSensitivities();
}
```

## Lift-target mode

With a lift target, the derivatives include the change in $\alpha$ that
keeps $C_l$ on target, so `dcl_dx`, `dcl_dRe`, and `dcl_dMa` are zero. The
derivative fails, with an exception, if $C_l$ doesn't change with $\alpha$.

In a viscous lift-target solve, the wake is laid out at the starting angle
and not moved (`oper.redowake = false`, as in mfoil). The solution then
depends on the starting angle, and `d*_dalpha` is the derivative with
respect to that starting angle. With `redowake = true`, the wake follows the
solved angle and the starting angle drops out.

## Near-field pressure drag

`cdpi` integrates surface pressure in the drag direction, in both inviscid
and viscous mode. It is separate from `cd`, which comes from the wake
momentum deficit and is zero in inviscid mode. The Philote server exposes
`cl`, `cd`, and `cm` only.

## Accuracy

Every derivative is tested against complex-step reference data from mfoil,
which is exact to machine precision. End-to-end comparisons use a tolerance
of `2e-11 + 2e-8 × |reference|`. Influence coefficients between distant
panels use a quadrature that avoids cancellation on short panels, and their
derivatives are checked against 60-digit references.
