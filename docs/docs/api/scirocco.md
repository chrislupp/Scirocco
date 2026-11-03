---
sidebar_position: 2
---

# `Scirocco` class

`scirocco::Scirocco` is the single user-facing entry point. It wraps an
internal `FoilState` and exposes the most common operations as members.
Header:
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h).

```cpp
#include <scirocco/scirocco.h>
```

## Construction

| Constructor | Use it when |
| --- | --- |
| `Scirocco(const std::string& naca = "0012", int npanel = 199)` | Generating a NACA 4- or 5-digit airfoil |
| `Scirocco(const Eigen::MatrixXd& coords, int npanel = 199)` | Bringing your own `2 × N` point cloud (TE-lower → LE → TE-upper) |

Both constructors run the curvature-graded paneler immediately, so after
construction `foil.GetFoilCoords()` is populated.

## Operating point

```cpp
void SetOper(double alpha = 0.0,
             double cl    = -1.0,
             double Re    = 0.0,
             double Ma    = 0.0);
```

| Argument | Notes |
| --- | --- |
| `alpha` | Angle of attack, degrees |
| `cl`    | Target lift coefficient. Any finite negative value (conventionally `-1.0`) selects prescribed-angle mode; a non-negative value requests a solve at the target lift. |
| `Re`    | Reynolds number based on chord. Any finite nonpositive value selects inviscid analysis. |
| `Ma`    | Freestream Mach. Activates the Kármán–Tsien correction. |

`SetOper` validates and stores the operating conditions and invalidates previous
convergence and sensitivity state. Thermodynamic quantities are initialized
during the subsequent solve. See [Input and output conventions](./conventions)
for valid ranges and normalization.

## Solve

```cpp
bool Solve();
```

Performs the inviscid linear solve, plus (if `Re > 0`) the coupled BL
Newton iteration. Returns `true` on success and `false` if convergence fails,
including an inviscid lift-target solve that cannot reach its target. Invalid
inputs and nonfinite results can throw exceptions. Check the return value
before using results; see [Troubleshooting](../getting-started/troubleshooting).

## Result accessors — scalars

| Method | Returns |
| --- | --- |
| `GetCl()`, `GetCd()`, `GetCm()` | Force coefficients |
| `GetCdf()`, `GetCdp()` | Friction- and pressure-drag components |
| `GetAlpha()` | The $\alpha$ that produced these results (relevant in Cl-mode) |
| `IsConverged()` | Mirrors the `Solve()` return value |

## Result accessors — distributions (`const Eigen::VectorXd&`)

| Method | Quantity |
| --- | --- |
| `GetCp()` | Pressure coefficient (Kármán–Tsien-corrected at $M_\infty > 0$) |
| `GetCf()` | Skin friction coefficient |
| `GetTheta()` | Momentum thickness $\theta$ |
| `GetDeltaStar()` | Displacement thickness $\delta^*$ |

Here `N = npanel + 1` and `N_w` is the number of wake nodes. After a
successful solve:

| Method | Inviscid size | Viscous size |
| --- | --- | --- |
| `GetCp()` | `N` | `N` |
| `GetCf()`, `GetTheta()`, `GetDeltaStar()` | Empty | `N + N_w` |

For boundary-layer arrays, airfoil nodes come first, followed by wake nodes.
Wake coordinates are empty in inviscid mode.

## Geometry accessors

| Method | Returns |
| --- | --- |
| `GetFoilCoords()` | `2 × N` airfoil panel-node coordinates |
| `GetWakeCoords()` | `2 × N_w` wake-node coordinates |

## Geometry modification

| Method | Effect |
| --- | --- |
| `GeomFlap(const Vector2d& hinge, double eta)` | Deflect everything aft of `hinge` by `eta` degrees |
| `GeomAddCamber(const MatrixXd& xzcamb)` | Add a $(x, z)$ camber increment, preserving thickness |
| `GeomDerotate()` | Re-align LE–TE with the $x$-axis |

All three re-panel automatically. The
[parameterized geometry tutorial](../tutorials/geometry-mods) demonstrates
the corresponding Python/Philote workflow; the methods above are the direct
C++ path.

## Sensitivities

```cpp
Sensitivity ComputeSensitivities();
```

Call after `Solve()`. Returns the [`Sensitivity`](./sensitivity) struct
and caches it internally so the per-output getters work:

| Method | Returns |
| --- | --- |
| `GetDclDx()`, `GetDcdDx()`, `GetDcmDx()` | `2 × N` matrices — derivative w.r.t. each panel-node coordinate |
| `GetDclDalpha()`, `GetDcdDalpha()`, `GetDcmDalpha()` | per-degree derivatives |
| `GetDclDMa()`, `GetDcdDMa()`, `GetDcmDMa()` | derivatives w.r.t. Mach |
| `GetDclDRe()`, `GetDcdDRe()`, `GetDcmDRe()` | derivatives w.r.t. Reynolds |

See the [sensitivities tutorial](../tutorials/sensitivities) for the
Python/Philote derivative interface and a finite-difference verification
recipe.

## Reaching the internal state

```cpp
FoilState&        GetState();
const FoilState&  GetState() const;
```

Use this to set fields like `param.verb`, `param.ncrit`, `param.use_ptc`,
or to read deeper distributions (`post.Ret`, `post.Hk`, etc.) that aren't
exposed as direct accessors.

```cpp
foil.GetState().param.verb   = 1;     // print solver progress
foil.GetState().param.ncrit  = 11.0;  // cleaner tunnel
```

See [`FoilState`](./foil-state) for the full struct layout.

## Minimal example

```cpp
#include <iostream>
#include <scirocco/scirocco.h>

int main() {
    scirocco::Scirocco foil("0012", 199);
    foil.SetOper(5.0, -1.0, 1.0e6, 0.0);
    if (foil.Solve()) {
        std::cout << "Cl = " << foil.GetCl()
                  << "  Cd = " << foil.GetCd() << "\n";
    }
}
```
