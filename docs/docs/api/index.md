---
title: C++ reference
sidebar_label: Overview
description: The public C++ API of Scirocco.
---

# C++ reference

Include one header and use one class:

```cpp
#include <scirocco/scirocco.h>
```

[`scirocco::Scirocco`](./scirocco) covers geometry, solving, results, and
gradients. Everything is in the `scirocco` namespace, and arrays are Eigen
types (`Eigen::VectorXd`, `Eigen::MatrixXd`, `Eigen::Vector2d`).

| Page | What's on it |
| --- | --- |
| [`Scirocco`](./scirocco) | The class you use for everything |
| [Sensitivities](./sensitivity) | The `Sensitivity` struct and gradient getters |
| [Geometry functions](./geometry) | NACA generation, coordinates, flaps, camber |
| [`FoilState`](./foil-state) | The full solver state, for reading beyond the getters |
| [Solver parameters](./solver-params) | Tolerances, iteration limits, damping, closures |
| [`Hdf5Writer`](./hdf5-writer) | Optional HDF5 output |
| [Internals](./internals) | Solver kernel headers, for contributors |
| [Conventions](./conventions) | Node order, units, valid ranges, array sizes |
| [Build options](./build-options) | CMake options |
| [`KulfanDiscipline`](./kulfan-discipline) | The CST geometry discipline behind the Python server |

Under the hood, `Scirocco` holds a [`FoilState`](./foil-state), and the
solver is a set of free functions that operate on it. Reach for those only
when the class doesn't expose what you need.
