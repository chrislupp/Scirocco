---
sidebar_position: 2
---

# Quickstart

A 5-minute end-to-end run: build the library, build the bundled NACA example,
and inspect its alpha-sweep output.

## Prerequisite: a configured build

If you haven't already, follow the [Installation](./installation) page through
`cmake --build`. You should have an executable at
`scirocco/build-default/naca_example`.

## Run the example

```bash
./scirocco/build-default/naca_example
```

You should see a NACA 0012 alpha sweep from $-4°$ to $12°$ in $2°$ steps:

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

$C_d$ is zero because the example runs **inviscid** ($Re = 0$) — D'Alembert's
paradox in action. To see non-zero drag, you turn on viscosity, which is the
subject of the next tutorial.

## What the example is doing

The whole source is 35 lines; here is the load-bearing part:

```cpp
#include <scirocco/scirocco.h>

int main() {
    // NACA 0012 with 199 panels (curvature-graded by default)
    scirocco::Scirocco foil("0012", 199);

    for (double alpha = -4.0; alpha <= 12.0; alpha += 2.0) {
        // alpha [deg], target Cl (-1 means "use alpha"), Re, Mach
        foil.SetOper(alpha, -1.0, 0.0, 0.0);
        foil.Solve();

        std::cout << alpha << "  " << foil.GetCl() << "  "
                  << foil.GetCd() << "  " << foil.GetCm() << "\n";
    }
}
```

Three calls do all the work:

| Call | What happens |
| --- | --- |
| `Scirocco("0012", 199)` | Builds NACA 0012 coordinates, runs cubic-spline curvature-graded paneling to 199 nodes. |
| `SetOper(α, Cl_target, Re, Ma)` | Stores operating conditions. `Cl_target = -1.0` is the sentinel for "use the supplied α directly"; passing a positive value instead triggers an α-search loop to hit that lift coefficient. |
| `Solve()` | Builds the AIC matrix, applies the Kutta condition, solves for the vortex strengths. With `Re > 0` it also Newton-couples the viscous boundary-layer system to the inviscid solve. Returns `true` on convergence. |

Result accessors (`GetCl`, `GetCd`, `GetCm`, `GetCdf`, `GetCdp`, plus the
distribution getters `GetCp`, `GetCf`, `GetTheta`, `GetDeltaStar`) read out of
the post-processing struct that `Solve()` populates.

## Try a viscous run

Edit `naca_example.cpp` (or write your own program against the library), and
swap the `SetOper` call:

```cpp
foil.SetOper(alpha, -1.0, 1.0e6, 0.0);   // Re = 1e6, incompressible
```

Rebuild and re-run. Drag is no longer zero:

```text
  Alpha      CL        CD        CM
  -----   -------   -------   -------
 0.0000    0.0000    0.0055    0.0000
 2.0000    0.2410    0.0058   -0.0001
 4.0000    0.4839    0.0066   -0.0005
 ...
```

That's the friction + pressure drag from the coupled boundary-layer solve.

## Where to next

- [Inviscid NACA tutorial](../tutorials/inviscid-naca) — same example, more
  annotation, with `GetCp` and `GetCf` distributions plotted.
- [Viscous analysis tutorial](../tutorials/viscous-analysis) — the viscous
  swap above, with a $Re$ sweep and drag breakdown ($C_{d,f}$ vs $C_{d,p}$).
- [Concepts → Panel method](../concepts/panel-method) — what `Solve()` is
  actually doing.
