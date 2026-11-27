---
title: Run an analysis
description: Angle-of-attack sweeps, lift targets, Mach number, and surface distributions in C++.
---

# Run an analysis

Everything goes through one class, `scirocco::Scirocco`. Create it with a
geometry, call `SetOper`, call `Solve`, and read the results.

```cpp
#include <scirocco/scirocco.h>

scirocco::Scirocco foil("2412", 160);   // NACA 2412, 160 panels
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);    // alpha, Cl target, Re, Mach
bool ok = foil.Solve();
```

Always check what `Solve()` returns. A failed solve leaves the results
undefined.

## Sweep angle of attack

Reuse one object and call `SetOper` for each point:

```cpp
scirocco::Scirocco foil("2412", 160);

std::printf(" alpha      Cl       Cd       Cm\n");
for (double alpha = -2.0; alpha <= 10.0; alpha += 2.0) {
    foil.SetOper(alpha, -1.0, 1.0e6, 0.0);
    if (!foil.Solve()) {
        std::printf("%6.1f  did not converge\n", alpha);
        continue;
    }
    std::printf("%6.1f  %7.4f  %7.5f  %7.4f\n", alpha, foil.GetCl(),
                foil.GetCd(), foil.GetCm());
}
```

```output
 alpha      Cl       Cd       Cm
  -2.0   0.0231  0.00657  -0.0542
   0.0   0.2354  0.00563  -0.0516
   2.0   0.4491  0.00577  -0.0480
   4.0   0.7202  0.00693  -0.0585
   6.0   0.9028  0.00902  -0.0507
   8.0   1.0844  0.01224  -0.0438
  10.0   1.2628  0.01547  -0.0370
```

Each point starts cold, so its result doesn't depend on the points before
it. For long polars, `foil.GetState().param.warm_start = true` starts each
point from the last one, which usually saves a few Newton iterations. Near stall
that can land on a different converged solution than a cold start; see
[warm starts](../concepts/viscous-coupling#warm-starts).

The same sweep from −4° to 14° in 1° steps, with the inviscid lift for
comparison:

<PolarPlot />

## Inviscid or viscous

`Re` picks the model:

- `Re > 0` couples the panel method to the boundary layer and wake. You get
  drag, skin friction, and boundary-layer thicknesses.
- `Re = 0` solves potential flow only. It's much faster, drag is zero, and
  there is no boundary layer or wake.

## Solve for a lift coefficient

Pass a non-negative `cl` and Scirocco finds the angle of attack:

```cpp
scirocco::Scirocco foil("2412", 160);
foil.SetOper(0.0, 0.5, 1.0e6, 0.0);   // find alpha for Cl = 0.5
if (foil.Solve()) {
    std::printf("alpha = %.4f deg, Cl = %.4f, Cd = %.5f\n",
                foil.GetAlpha(), foil.GetCl(), foil.GetCd());
}
```

```output
alpha = 2.3754 deg, Cl = 0.5000, Cd = 0.00597
```

The `alpha` you pass is the starting guess. Pass `cl = -1.0` to go back to
fixed-angle mode.

## Add Mach number

The fourth argument turns on the Kármán–Tsien compressibility correction:

```cpp
for (double mach : {0.0, 0.3, 0.5, 0.6}) {
    foil.SetOper(2.0, -1.0, 3.0e6, mach);
    if (foil.Solve()) {
        std::printf("Ma = %.1f  Cl = %.4f  Cd = %.5f\n",
                    mach, foil.GetCl(), foil.GetCd());
    }
}
```

```output
Ma = 0.0  Cl = 0.4617  Cd = 0.00506
Ma = 0.3  Cl = 0.4872  Cd = 0.00522
Ma = 0.5  Cl = 0.5437  Cd = 0.00558
Ma = 0.6  Cl = 0.5945  Cd = 0.00591
```

The correction doesn't model shocks. Keep to subsonic flow without strong
supersonic pockets; see [Compressibility](../concepts/compressibility).

## Read the distributions

After a solve, surface and boundary-layer data are available as Eigen
vectors:

```cpp
const Eigen::MatrixXd& xz = foil.GetFoilCoords();  // 2 x N panel nodes
const Eigen::VectorXd& cp = foil.GetCp();          // N
const Eigen::VectorXd& cf = foil.GetCf();          // N + wake nodes
```

```output
161 airfoil nodes, 27 wake nodes
cp: 161 values, cf: 188 values
```

| Accessor | Size | Contents |
| --- | --- | --- |
| `GetCp()` | `N` | Pressure coefficient at each panel node |
| `GetCf()` | `N + Nw` | Skin friction, airfoil nodes then wake nodes |
| `GetTheta()` | `N + Nw` | Momentum thickness |
| `GetDeltaStar()` | `N + Nw` | Displacement thickness |
| `GetFoilCoords()` | `2 × N` | Node coordinates, lower trailing edge → leading edge → upper trailing edge |
| `GetWakeCoords()` | `2 × Nw` | Wake node coordinates |

`N` is the number of panels plus one. In an inviscid solve the
boundary-layer arrays and wake are empty. More quantities, such as the shape
factor `Hk` and edge velocity `ue`, are in `GetState().post`; see
[`FoilState`](../api/foil-state).

## See the solver working

Set the verbosity before solving to print the residual each Newton step:

```cpp
foil.GetState().param.verb = 1;
```

A healthy viscous solve converges in under ten steps. If it doesn't
converge, see [Troubleshooting](../getting-started/troubleshooting#the-solve-doesnt-converge).

<Next items={[
  { to: "/docs/guides/geometry", title: "Define the geometry", body: "NACA sections, your own coordinates, flaps, and camber." },
  { to: "/docs/guides/gradients", title: "Compute gradients", body: "Derivatives of the results for optimization." },
]} />
