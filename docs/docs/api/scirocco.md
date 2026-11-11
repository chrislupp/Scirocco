---
title: Scirocco
description: Reference for the scirocco::Scirocco class.
---

# `Scirocco`

```cpp
#include <scirocco/scirocco.h>
```

The one class you need. It holds the geometry, the operating point, and the
solution. Header:
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h).

```cpp
scirocco::Scirocco foil("2412", 160);
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
if (foil.Solve()) {
    std::cout << foil.GetCl() << ' ' << foil.GetCd() << '\n';
}
```

## Constructors

```cpp
explicit Scirocco(const std::string& naca = "0012", int npanel = 199);
explicit Scirocco(const Eigen::MatrixXd& coords, int npanel = 199);
```

| Argument | |
| --- | --- |
| `naca` | A 4-digit or `2X0YY` 5-digit NACA designation |
| `coords` | `2 × M` points, lower trailing edge → leading edge → upper trailing edge |
| `npanel` | Panel count, 4 to 4096. The section gets `npanel + 1` nodes. |

Both constructors panel the section straight away. See
[Define the geometry](../guides/geometry).

## Operating point

```cpp
void SetOper(double alpha = 0.0, double cl = -1.0, double Re = 0.0, double Ma = 0.0);
```

| Argument | |
| --- | --- |
| `alpha` | Angle of attack in degrees. With a lift target, the starting guess. |
| `cl` | Negative: solve at `alpha`. Zero or positive: find the `alpha` that gives this lift coefficient. |
| `Re` | Chord Reynolds number. Positive: viscous. Zero or negative: inviscid. |
| `Ma` | Freestream Mach number, `0 <= Ma < 1` |

`SetOper` checks the inputs and clears any previous solution and gradients.

## Solve

```cpp
bool Solve();
```

Runs the inviscid solve and, when `Re > 0`, the coupled viscous Newton
iteration. Returns `false` if it doesn't converge or a lift target can't be
reached; then the results are not valid. Throws on invalid inputs or
non-finite results.

## Results

| Getter | Returns |
| --- | --- |
| `GetCl()`, `GetCd()`, `GetCm()` | Force and moment coefficients; moment about (0.25, 0) |
| `GetCdf()`, `GetCdp()` | Friction and pressure parts of `Cd` |
| `GetAlpha()` | Angle of attack of the solution, useful in lift-target mode |
| `IsConverged()` | Same as the last `Solve()` return |
| `GetCp()` | `VectorXd`, `N`: pressure coefficient at each node |
| `GetCf()` | `VectorXd`, `N + Nw`: skin friction, airfoil then wake |
| `GetTheta()` | `VectorXd`, `N + Nw`: momentum thickness |
| `GetDeltaStar()` | `VectorXd`, `N + Nw`: displacement thickness |
| `GetFoilCoords()` | `MatrixXd`, `2 × N`: panel nodes |
| `GetWakeCoords()` | `MatrixXd`, `2 × Nw`: wake nodes |

In an inviscid solve, `GetCd()` is zero and the boundary-layer and wake
arrays are empty.

## Geometry changes

| Method | Effect |
| --- | --- |
| `GeomFlap(const Vector2d& hinge, double eta)` | Deflect everything aft of `hinge` by `eta` degrees, positive down |
| `GeomAddCamber(const MatrixXd& xzcamb)` | Add a camber-line increment, keeping thickness |
| `GeomDerotate()` | Rotate so the chord line is horizontal |

Each one re-panels the section. Call `SetOper` and `Solve` again afterwards.

## Gradients

```cpp
Sensitivity ComputeSensitivities();
```

Call after a successful `Solve()`. Returns the
[`Sensitivity`](./sensitivity) struct and caches it for these getters:

| Getter | Returns |
| --- | --- |
| `GetDclDx()`, `GetDcdDx()`, `GetDcmDx()` | `2 × N` derivatives with respect to node x and z |
| `GetDclDalpha()`, `GetDcdDalpha()`, `GetDcmDalpha()` | Per degree |
| `GetDclDRe()`, `GetDcdDRe()`, `GetDcmDRe()` | Per unit Reynolds number |
| `GetDclDMa()`, `GetDcdDMa()`, `GetDcmDMa()` | Per unit Mach number |

The getters throw if there's no valid cached gradient. See
[Compute gradients](../guides/gradients).

## Solver state

```cpp
FoilState&       GetState();
const FoilState& GetState() const;
```

Use it to set [solver parameters](./solver-params) before solving, or to read
quantities without a getter:

```cpp
foil.GetState().param.verb  = 1;      // print the residual each step
foil.GetState().param.ncrit = 11.0;   // later transition

const scirocco::Scirocco& solved = foil;
const Eigen::VectorXd& hk = solved.GetState().post.Hk;   // shape factor
```

The non-const overload clears cached gradients, since you might change the
state through it. Read through a `const` reference to keep them. See
[`FoilState`](./foil-state) for the layout.
