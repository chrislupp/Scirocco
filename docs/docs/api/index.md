---
sidebar_position: 1
---

# API Reference

Hand-written reference for the public surface of `scirocco::`. Internal
solver kernels are summarized in a single index page with links to source.

## Pages

| Page | Contents |
| --- | --- |
| [`Scirocco`](./scirocco) | The user-facing entry-point class |
| [`FoilState`](./foil-state) | The aggregate state struct + sub-structs (`Panel`, `Post`, `Oper`, `Param`) |
| [Geometry](./geometry) | `NacaPoints`, `SetCoords`, `GeomFlap`, `GeomAddCamber`, `GeomDerotate` |
| [Sensitivities](./sensitivity) | `ComputeSensitivities`, `Sensitivity` struct, and the gradient getters |
| [`Hdf5Writer`](./hdf5-writer) | Optional HDF5 iteration output |
| [`KulfanDiscipline`](./kulfan-discipline) | CST parameterization + Philote gRPC interface |
| [Internals](./internals) | Index of solver kernel headers (boundary-layer, transition, coupled solver, adjoint, ...) with source links |

## Conventions

- All public symbols live in the `scirocco::` namespace.
- The single user-facing header is `<scirocco/scirocco.h>`. Pulling that in
  transitively brings the data structures you need for typical usage.
- Eigen types are used directly (`Eigen::VectorXd`, `Eigen::MatrixXd`,
  `Eigen::Vector2d`).
- Free functions in `core/` operate on `FoilState` by reference; the
  `Scirocco` class is a thin wrapper that holds a `FoilState` and exposes the
  most common operations as member functions.

## Where to start

If you've just installed Scirocco, start at [`Scirocco`](./scirocco) — that's
the entire surface most users need. The other pages are there when you want to
reach below the wrapper for specific control.
