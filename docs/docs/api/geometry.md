---
sidebar_position: 4
---

# Geometry

Free functions for generating NACA airfoils, loading user coordinates,
and modifying geometry in place. Header:
[`scirocco/core/geometry.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/geometry.h).

```cpp
#include <scirocco/core/geometry.h>
```

All functions live in the `scirocco::` namespace and operate on `Geom`
or `FoilState` references. The same operations are wrapped as members on
[`Scirocco`](./scirocco) — use the wrappers for typical work, drop down
to the free functions when you want fine control without going through
the `Scirocco` constructor.

## NACA generation

```cpp
void NacaPoints(Geom& geom, const std::string& digits);
```

Generates a NACA 4- or 5-digit airfoil and writes the result into
`geom.xpoint`. Both families are supported:

```cpp
scirocco::Geom geom;
scirocco::NacaPoints(geom, "2412");   // NACA 4-digit, m=2%, p=40%, t=12%
scirocco::NacaPoints(geom, "23012");  // NACA 5-digit
```

The output contains 201 source points in a `2 × N` matrix, using a nonuniform
spacing law. The constructor subsequently re-panels these to the requested
count. Supported five-digit designations begin with `2X0`, with `X` from 1 to 5.

## User-supplied coordinates

```cpp
void SetCoords(Geom& geom, const Eigen::MatrixXd& coords);
```

Loads an arbitrary `2 × N` point cloud, ordered TE-lower → LE → TE-upper
(Scirocco ordering; check and reverse imported files if necessary). Internally this just copies the
matrix into `geom.xpoint`; the paneler that runs later does the actual
discretization.

```cpp
Eigen::MatrixXd coords = LoadDatFile("rg15.dat");   // your loader
scirocco::Geom geom;
scirocco::SetCoords(geom, coords);
```

## Trailing-edge information

```cpp
struct TEInfo {
    Vector2d t;      // unit tangent at TE (average of upper and lower)
    double   hTE;    // TE gap height
    double   dtdx;   // d(theta)/dx at TE — used in BL closure
    double   tcp;    // TE constant source-panel strength factor
    double   tdp;    // TE vortex-panel strength factor
};

TEInfo GetTEInfo(const Eigen::MatrixXd& x);
```

Used inside the inviscid solver to set up the Kutta closure and TE
constant source. Users rarely call it directly, but it's the right way to
inspect the TE gap and tangent if you need to verify geometry intake.

## Geometry modification

These three functions mutate `state.geom.xpoint` and re-panel
`state.foil` if panels already exist. All three are also exposed as
[`Scirocco`](./scirocco) member functions; the wrapper versions are
identical.

### Flap

```cpp
void GeomFlap(FoilState& state, const Vector2d& xzhinge, double eta);
```

Rotates everything aft of `xzhinge` by `eta` degrees. Positive `eta` is
flap-down. The hinge does not need to coincide with a panel node;
interpolation is automatic.

### Camber addition

```cpp
void GeomAddCamber(FoilState& state, const MatrixXd& xzcamb);
```

`xzcamb` is a `2 × Nc` camber-line **increment** sampled along the chord.
The increment is linearly interpolated and added to the existing mean
line; the thickness distribution is preserved.

### Derotation

```cpp
void GeomDerotate(FoilState& state);
```

Rotates the airfoil so the LE–TE line is parallel to the $x$-axis. Useful
for canonicalizing measured coordinates.

## Composition

The modifier operations compose in call order — each one re-panels
afterwards so the next call sees the updated discretization:

```cpp
scirocco::Scirocco foil("0012", 199);

Vector2d hinge(0.75, 0.0);
foil.GeomFlap(hinge, 10.0);           // (1) deflect
foil.GeomDerotate();                  // (2) re-canonicalize
foil.SetOper(0.0, -1.0, 1.0e6, 0.0);
foil.Solve();
```

## Panel geometry refresh

```cpp
void RefreshPanelGeometry(Panel& panel);
```

Header:
[`scirocco/core/paneling.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/paneling.h).

Rebuilds `panel.s` (arc-length array) and `panel.t` (unit tangent vectors)
from the current `panel.x` node coordinates without re-spacing the
panels. Use this after directly overwriting `foil.x` from an upstream
geometry discipline — for example, after `SiroccoDiscipline::Compute`
receives new coordinates from a `KulfanDiscipline` — so the boundary-layer
solver sees consistent geometric state without breaking upstream gradient
chains.

```cpp
// After overwriting panel.x with new CST coordinates:
scirocco::RefreshPanelGeometry(state.foil);
// Now panel.s and panel.t match panel.x — safe to call the BL solver.
```

`RefreshPanelGeometry` is called automatically by `SiroccoDiscipline::Compute`
when chained downstream of a geometry discipline; you only need to call
it directly when manipulating `panel.x` outside the discipline interface.

## See also

- Tutorial: [Supply coordinates from Python](../tutorials/custom-airfoil)
- Tutorial: [Modify parameterized geometry](../tutorials/geometry-mods)
