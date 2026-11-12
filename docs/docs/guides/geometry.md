---
title: Define the geometry
description: NACA sections, coordinate files, flaps, camber, and panel counts in C++.
---

# Define the geometry

Pass either a NACA designation or a coordinate matrix to the constructor.
Scirocco then re-panels the shape with nodes clustered where curvature is
high.

<PanelPlot />

## NACA sections

```cpp
scirocco::Scirocco a("0012");        // 4-digit, 199 panels (the default)
scirocco::Scirocco b("2412", 160);   // 4-digit, 160 panels
scirocco::Scirocco c("23012", 160);  // 5-digit
```

Five-digit sections must be the standard non-reflexed family, `2X0YY` with
`X` from 1 to 5. The generator gives a small trailing-edge gap, 0.25% of
chord for a 12% section.

## Your own coordinates

Pass a `2 × M` matrix: row 0 is x, row 1 is z. Points run from the **lower**
trailing edge, around the leading edge, to the **upper** trailing edge.
Use unit chord, with the leading edge near `x = 0`.

Most `.dat` files use the Selig format, which runs the other way: upper
trailing edge first. This loader reads one and reverses it:

```cpp
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <scirocco/scirocco.h>

// Selig .dat: a name line, then "x z" pairs from the upper trailing edge
// around to the lower one. Scirocco wants the reverse.
Eigen::MatrixXd ReadSelig(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::string name;
    std::getline(in, name);
    std::vector<double> x, z;
    for (double xi, zi; in >> xi >> zi;) {
        x.push_back(xi);
        z.push_back(zi);
    }
    const long n = static_cast<long>(x.size());
    Eigen::MatrixXd xz(2, n);
    for (long i = 0; i < n; ++i) {
        xz(0, i) = x[n - 1 - i];
        xz(1, i) = z[n - 1 - i];
    }
    return xz;
}

int main() {
    scirocco::Scirocco foil(ReadSelig("naca0012.dat"), 160);
    foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
    if (foil.Solve()) {
        std::printf("Cl = %.4f  Cd = %.5f\n", foil.GetCl(), foil.GetCd());
    }
}
```

```output
Cl = 0.4277  Cd = 0.00728
```

The input points don't need any particular spacing; the paneler
redistributes them. They do need to be a clean loop: no repeated
neighbouring points, no self-intersections, positive chord.

:::warning Leave a gap at the trailing edge
Viscous solves on coordinates with an **exactly closed** trailing edge often
fail to converge, and can converge to a wrong answer. Give the trailing edge
a small gap, as the built-in NACA generator does. This is a known issue
([#18](https://github.com/chrislupp/Scirocco/issues/18)).
:::

### Why a closed trailing edge fails

The limitation comes from mfoil, which Scirocco ports. The wake leaves the
airfoil along the normal to the trailing-edge gap. When the gap is zero,
that direction is undefined:

- the first wake point lands on the trailing edge itself;
- in potential flow that is a stagnation point, so the wake streamline has
  no direction to follow.

mfoil stops on the same coordinates with "Wrong wake direction". Scirocco
has no such check, so it carries on with a bad wake. Usually the solve then
fails, but sometimes it reports convergence at the wrong Cl.

Starting the wake along the trailing-edge bisector instead doesn't fix it.
mfoil gets past the wake construction, but the coupled viscous solve still
fails on NACA 0012 at 120–199 panels. Supporting a sharp trailing edge
properly would take a new formulation of the trailing-edge and wake
junction, not a port fix.

## Flaps, camber, and alignment

These modify the current shape and re-panel it. They apply in the order you
call them.

```cpp
scirocco::Scirocco foil("0012", 160);
foil.SetOper(2.0, -1.0, 1.0e6, 0.0);
foil.Solve();
std::printf("clean:       Cl = %.4f\n", foil.GetCl());

foil.GeomFlap(Eigen::Vector2d(0.75, 0.0), 10.0);  // hinge at 75% chord, 10° down
foil.SetOper(2.0, -1.0, 1.0e6, 0.0);
foil.Solve();
std::printf("flap 10 deg: Cl = %.4f\n", foil.GetCl());
```

```output
clean:       Cl = 0.2119
flap 10 deg: Cl = 0.8287
```

| Method | Effect |
| --- | --- |
| `GeomFlap(hinge, eta)` | Rotates everything aft of `hinge` by `eta` degrees; positive is down |
| `GeomAddCamber(xzcamb)` | Adds a camber-line increment (a `2 × K` matrix of x, Δz) and keeps the thickness |
| `GeomDerotate()` | Rotates the section so the chord line is horizontal |

After a flap, `GeomDerotate()` puts the chord line back on the x-axis, so
`alpha` is measured from the new chord.

## How many panels

The default is 199. More panels cost more time, mostly in the viscous solve
and the gradient; see [Performance](../concepts/performance). Whatever count
you pick, check that your answer doesn't move when you change it.

Panel counts must be between 4 and 4096.

<Next items={[
  { to: "/docs/guides/gradients", title: "Compute gradients", body: "Derivatives with respect to every node coordinate." },
  { to: "/docs/concepts/geometry", title: "Geometry and CST", body: "The paneler and the CST parameterization." },
]} />
