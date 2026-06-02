---
sidebar_position: 3
---

# Custom airfoil coordinates

Skip the NACA generator and hand Scirocco an arbitrary point cloud. Useful
for analyzing measured geometry, CST-fit airfoils, or shapes coming out of an
external CAD tool.

## The constructor

```cpp
scirocco::Scirocco(const Eigen::MatrixXd& coords, int npanel = 199);
```

`coords` is a `2 × N` matrix (row 0 is $x$, row 1 is $z$) ordered **TE-lower
→ LE → TE-upper** (the standard counter-clockwise loop). Internally
Scirocco calls `SetCoords` (see [API → Geometry](../api/geometry)) to load
the points, then runs the same cubic-spline curvature-graded paneler that
the NACA constructor uses, producing `npanel` evenly-resolved panel nodes.

## Example: round-trip a NACA 0012

To check that the constructor is doing the right thing, generate NACA 0012
points yourself and feed them back in. The two paths should produce the same
$C_l$ on the same operating point.

```cpp
#include <Eigen/Dense>

#include "scirocco/scirocco.h"

Eigen::MatrixXd MakeCircle(int n_per_side) {
    // Trivial example: a 12% thick parabolic-arc airfoil.
    int n = 2 * n_per_side - 1;
    Eigen::MatrixXd coords(2, n);

    for (int i = 0; i < n_per_side; ++i) {
        double s = static_cast<double>(i) / (n_per_side - 1);
        double x = 1.0 - s;
        double z = -0.06 * 4.0 * x * (1.0 - x);   // lower surface
        coords(0, i) = x;
        coords(1, i) = z;
    }
    for (int i = 1; i < n_per_side; ++i) {
        double s = static_cast<double>(i) / (n_per_side - 1);
        double x = s;
        double z = +0.06 * 4.0 * x * (1.0 - x);   // upper surface
        coords(0, n_per_side - 1 + i) = x;
        coords(1, n_per_side - 1 + i) = z;
    }
    return coords;
}

int main() {
    auto coords = MakeCircle(80);
    scirocco::Scirocco foil(coords, 199);

    foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
    foil.Solve();

    std::cout << "Cl = " << foil.GetCl()
              << "  Cd = " << foil.GetCd() << "\n";
}
```

## Point ordering — the gotcha

The single most common mistake is the **ordering convention**. Scirocco
expects:

```
index 0   : trailing edge, lower surface
index ... : along the lower surface forward to the leading edge
index k   : leading edge (or somewhere near it — the exact LE is found
            automatically; you do not need to put a node there)
index ... : along the upper surface aft to the trailing edge
index N-1 : trailing edge, upper surface
```

This is the same convention used by XFOIL `.dat` files and Selig's UIUC
database. If your CAD tool exports clockwise, reverse the matrix before
passing it in:

```cpp
Eigen::MatrixXd reordered = coords.rowwise().reverse();
```

If you get an obviously wrong $C_l$ (negative at $\alpha = 0$ on a symmetric
shape, for instance), reversed coordinates are the first thing to check.

## Closed vs. open trailing edges

Scirocco does not require the upper and lower TE points to coincide. A small
finite TE gap is handled by `GetTEInfo` automatically and contributes the
expected base-pressure terms. If your shape has a deliberately blunt TE,
just keep both points and Scirocco will treat the gap as a real geometric
feature.

## Recovering the paneled geometry

After construction you may want to see what Scirocco actually paneled (the
input points may not appear as nodes — the paneler chooses its own
distribution):

```cpp
const Eigen::MatrixXd& paneled = foil.GetFoilCoords();   // 2 × npanel
```

Plotting the input cloud and `paneled` together on the same axes is a quick
sanity check.

## Full source

There is no canonical example file for this workflow; the snippet above is
self-contained. The constructor lives in
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h),
and the underlying coordinate intake is `SetCoords` in
[`scirocco/core/geometry.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/geometry.h).

## Next: [Geometry modifications](./geometry-mods) — flaps, camber, derotation.
