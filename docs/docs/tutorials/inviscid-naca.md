---
sidebar_position: 1
---

# Inviscid NACA analysis

A guided walk through the simplest possible Scirocco program: a NACA 0012
alpha sweep with no viscous coupling. This mirrors
[`scirocco/examples/naca_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/naca_example.cpp).

By the end you will:

- Instantiate a NACA airfoil and pick a panel count
- Sweep angle of attack and read out $C_l$, $C_d$, $C_m$
- Understand why inviscid drag is always zero

## The program

```cpp
#include <iomanip>
#include <iostream>

#include "scirocco/scirocco.h"

int main() {
    scirocco::Scirocco foil("0012", 199);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Alpha      CL        CD        CM\n";

    for (double alpha = -4.0; alpha <= 12.0; alpha += 2.0) {
        foil.SetOper(alpha, -1.0, 0.0, 0.0);   // alpha, cl_target, Re, Ma
        foil.Solve();

        std::cout << std::setw(7) << alpha << "  "
                  << std::setw(8) << foil.GetCl() << "  "
                  << std::setw(8) << foil.GetCd() << "  "
                  << std::setw(8) << foil.GetCm() << "\n";
    }
}
```

## Walkthrough

### Construction

```cpp
scirocco::Scirocco foil("0012", 199);
```

The first argument is a NACA digit string. Both 4- and 5-digit families work
(`"23012"` is a valid 5-digit). Internally this calls `NacaPoints` to generate
the analytical airfoil, then runs curvature-graded paneling with cubic-spline
arc-length parameterization to produce 199 panel nodes.

Why 199? An odd panel count gives a single node at the leading edge, which the
stagnation tracker prefers. You can pick any positive integer; 99 is typical
for cheap inviscid sweeps, 199 is comfortable for viscous work, and going much
above ≈ 400 starts hitting diminishing returns versus runtime.

### Operating point

```cpp
foil.SetOper(alpha, -1.0, 0.0, 0.0);
```

`SetOper(α, Cl_target, Re, Ma)`:

| Slot | Value | Meaning |
| --- | --- | --- |
| `alpha` | `-4 … 12` | Angle of attack in degrees |
| `cl_target` | `-1.0` | Sentinel: "drive at the supplied α; do not search for a target lift". A non-negative value here switches to **Cl-mode**, where Scirocco runs an outer secant loop to find the α that yields that $C_l$ — useful for fixed-lift drag polars. |
| `Re` | `0.0` | Zero disables viscous coupling. Any positive value turns it on. |
| `Ma` | `0.0` | Zero is incompressible. Up to ≈ 0.7 the Kármán–Tsien correction applies. |

### Solve

```cpp
foil.Solve();
```

For an inviscid run this builds the aerodynamic influence-coefficient (AIC)
matrix, applies the Kutta condition at the trailing edge, solves the linear
system for the vortex strengths, and integrates pressure to get force
coefficients. `Solve()` returns `true` on convergence; for pure inviscid it is
always `true` because the system is linear.

### Read results

```cpp
foil.GetCl();   // lift coefficient
foil.GetCd();   // total drag (always zero for inviscid — see below)
foil.GetCm();   // pitching moment about c/4
```

For richer output, also grab the distributions:

```cpp
const auto& cp = foil.GetCp();           // 2 × N — pressure coefficient
const auto& x  = foil.GetFoilCoords();   // 2 × N — node x,z
```

You can pipe `(x, cp)` straight into Matplotlib / Tecplot / PGFPlots to draw
the classic chord-wise $C_p$ plot.

## Sample output

```text
=== scirocco Example: NACA Airfoil Analysis ===

  Alpha      CL        CD        CM
  -----   -------   -------   -------
-4.0000   -0.4823    0.0000   -0.0000
-2.0000   -0.2407    0.0000    0.0000
 0.0000    0.0000    0.0000    0.0000
 2.0000    0.2407    0.0000   -0.0000
 4.0000    0.4823    0.0000   -0.0000
 ...
=== Analysis Complete ===
```

Two things to notice:

1. **$C_l$ is linear in $\alpha$ with slope $\approx 0.121$/deg.** This is the
   thin-airfoil prediction $C_{l,\alpha} = 2\pi$ /rad
   $\approx 0.110$ /deg, plus a thickness correction the panel method picks
   up.
2. **$C_d$ is identically zero.** D'Alembert's paradox: a 2D body in inviscid,
   incompressible, irrotational flow produces no drag. To see real drag,
   continue to the [Viscous analysis tutorial](./viscous-analysis).

## Full source

Canonical copy: [`scirocco/examples/naca_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/naca_example.cpp)
in the repo, built as the `naca_example` CMake target.

## Next: [Viscous analysis](./viscous-analysis) — turning on the boundary-layer solver.
