---
title: Build options
description: CMake options for the Scirocco build.
---

# Build options

Pass these to CMake at configure time, for example
`-DSCIROCCO_HDF5_OUTPUT=ON`.

| Option | Default | What it does |
| --- | --- | --- |
| `SCIROCCO_BUILD_TESTS` | `ON` | Solver, geometry, Jacobian, and adjoint tests. The reference suites also need HDF5. |
| `SCIROCCO_BUILD_EXAMPLES` | `ON` | `quickstart`, `naca_example`, `viscous_example`, `bench`, and `docs_figure_data` |
| `SCIROCCO_BUILD_CLI` | `ON` | `fauXFOIL`, the [XFOIL-style command line](../guides/fauxfoil.md) |
| `SCIROCCO_CLI_GRAPHICS` | `ON` | Give `fauXFOIL` a [plot window](../guides/fauxfoil.md#plots) when Qt 6 Widgets is found; text-only otherwise |
| `SCIROCCO_HDF5_OUTPUT` | `OFF` | The `scirocco_hdf5` library and `hdf5_output_example` |
| `SCIROCCO_BUILD_PHILOTE` | `OFF` | The Philote gRPC servers. Needs Protobuf, gRPC, and Philote-Cpp. |
| `SCIROCCO_NATIVE_ARCH` | `ON` | Compile Release builds with `-march=native` |
| `SCIROCCO_USE_BLAS` | `ON` | Use BLAS/LAPACK (the Apple Accelerate framework on macOS, OpenBLAS or MKL elsewhere) for the dense part of the coupled solve; falls back to Eigen if none is found |
| `SCIROCCO_SANITIZERS` | `OFF` | AddressSanitizer and UndefinedBehaviorSanitizer (GCC or Clang) |
| `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS` | `OFF` | Extra tests that isolate each piece of the adjoint chain |
| `SCIROCCO_BUILD_PROBES` | `OFF` | 14 `probe_*` diagnostic programs for adjoint debugging |

`SCIROCCO_USE_BLAS` speeds up viscous solves by 10-15% at 200 panels and
more at higher panel counts. CMake reports which library it found. With a
static library, the installed package configuration finds the same LAPACK
for consumers.

`SCIROCCO_NATIVE_ARCH` makes binaries specific to the build machine's CPU.
Turn it off when you build binaries to run elsewhere, such as Docker images
or release artifacts.

Build directories are disposable. If CMake's cache or test discovery gets
confused, delete the directory and configure again.
