---
sidebar_position: 4
---

# Geometry modifications

Once an airfoil exists you can deflect a flap, layer additional camber on
top of it, or remove an offset rotation. All three operations mutate the
geometry in place and re-panel automatically, so subsequent `Solve()` calls
see the modified shape.

## Flap deployment

```cpp
void GeomFlap(const Eigen::Vector2d& xzhinge, double eta);
```

Hinges on a single point and rotates everything aft of it by `eta` degrees.
Positive `eta` is **flap-down** (trailing edge moves toward $-z$, increasing
camber and lift).

```cpp
scirocco::Scirocco foil("0012", 199);

Eigen::Vector2d hinge(0.75, 0.0);     // 75% chord
foil.GeomFlap(hinge, 10.0);           // 10° down

foil.SetOper(0.0, -1.0, 1.0e6, 0.0);
foil.Solve();

std::cout << "Flapped Cl at alpha=0: " << foil.GetCl() << "\n";
```

Compared to the clean 0012 (where $C_l(0°) = 0$), you should see a positive
$C_l$ — the canonical effect of a positive flap deflection.

Implementation notes:

- The hinge does not have to coincide with a panel node; Scirocco interpolates.
- Multiple `GeomFlap` calls compose, but each one re-panels the modified
  geometry. For very large deflections that introduce kinks, refine panels.
- A negative `eta` is a flap-up (spoiler-like) deflection.

## Adding camber

```cpp
void GeomAddCamber(const Eigen::MatrixXd& xzcamb);
```

`xzcamb` is a `2 × Nc` matrix sampling a camber-line **increment** along
the chord. The increment is added to the existing mean line; the thickness
distribution is preserved. Useful for retrofitting reflex or extra rear
loading onto a base section.

```cpp
// A small parabolic reflex camber centered at x = 0.7
int Nc = 51;
Eigen::MatrixXd camb(2, Nc);
for (int i = 0; i < Nc; ++i) {
    double x = static_cast<double>(i) / (Nc - 1);
    camb(0, i) = x;
    camb(1, i) = -0.01 * std::pow(x - 0.7, 2) * (x > 0.4 ? 1.0 : 0.0);
}

scirocco::Scirocco foil("0012", 199);
foil.GeomAddCamber(camb);
```

The function discards points outside $[0, 1]$ and treats the camber-line as
piecewise linear in $x$.

## Derotation

```cpp
void GeomDerotate();
```

Removes any pitch offset by re-aligning the chord (LE–TE line) with the
$x$-axis. This is mostly useful when ingesting measured coordinates that
were captured at an unknown installation angle, or when chaining geometry
modifications that introduce a net rotation.

```cpp
scirocco::Scirocco foil(measured_coords, 199);
foil.GeomDerotate();    // canonical orientation
```

After derotation, $\alpha$ in `SetOper` again means **angle relative to the
chord line**, which is the usual convention.

## Order matters

Geometry modifications run against the **current** shape:

```cpp
foil.GeomDerotate();              // (1) canonicalize
foil.GeomAddCamber(reflex);       // (2) layer camber on the clean shape
foil.GeomFlap(hinge, 10.0);       // (3) deflect the cambered shape
```

If you want to swap stages — e.g. flap before camber — call them in the
desired order. Each call re-panels at the end, so the next call sees a fully
discretized airfoil.

## Full source

The class wrappers live in
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h);
the underlying free functions are in
[`scirocco/core/geometry.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/geometry.h).

## Next: [HDF5 output](./hdf5-output) — persist iteration history to disk.
