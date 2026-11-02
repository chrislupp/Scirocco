# Scirocco

[![Build](https://github.com/chrislupp/Scirocco/actions/workflows/build.yml/badge.svg)](https://github.com/chrislupp/Scirocco/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/chrislupp/Scirocco/branch/main/graph/badge.svg)](https://codecov.io/gh/chrislupp/Scirocco/branch/main)

A modern C++ library for inviscid and viscous aerodynamic analysis of airfoils.

## Overview

Scirocco provides computational tools for analyzing airfoil performance across different operating conditions. It predicts aerodynamic coefficients (lift, drag, moment) and boundary layer distributions using panel methods coupled with integral boundary layer equations.

Key capabilities:
- Inviscid (potential flow) analysis using linear vortex panel methods
- Viscous analysis with coupled boundary layer solver
- Laminar-to-turbulent transition prediction (e^N method)
- Compressibility corrections for subsonic flow
- NACA 4-digit and 5-digit airfoil generation
- Geometry modification (flaps, camber, derotation)
- Philote-Cpp gRPC disciplines with OpenMDAO/Python bindings

## Requirements

- CMake 3.18+
- C++17 compatible compiler
- Eigen 3.4+ (linear algebra library)

Optional:
- Google Test (for unit tests)
- HDF5 (for reference validation tests)
- Protobuf, gRPC, and Philote-Cpp (for the Python/MDO server interface)

## Building

From the repository root:

```bash
cmake -S scirocco -B scirocco/build-default -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build-default -j
```

Local build directories are ignored by git and are disposable; delete and
regenerate them whenever CMake or test discovery state gets stale.

Build options:

| Option | Default | Purpose |
|---|---|---|
| `SCIROCCO_BUILD_PHILOTE` | `OFF` | Build gRPC servers; requires Protobuf, gRPC, and Philote-Cpp |
| `SCIROCCO_BUILD_TESTS` | `ON` | Build the default test suite (4 gtest executables — primal + Jacobian unit + end-to-end adjoint vs CS) |
| `SCIROCCO_BUILD_EXAMPLES` | `ON` | Build example programs |
| `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS` | `OFF` | Build chain-piece sensitivity tests (`sensitivity_ping`, `viscous_chain`, `viscous_residual_deriv`, `viscous_output_ping`) — useful for localizing adjoint regressions; otherwise redundant with the end-to-end test |
| `SCIROCCO_BUILD_PROBES` | `OFF` | Build the 14 standalone `probe_*` diagnostic binaries |
| `SCIROCCO_COVERAGE` | `OFF` | Instrument library + tests with `-fprofile-instr-generate -fcoverage-mapping` for `llvm-cov` |
| `SCIROCCO_NATIVE_ARCH` | `ON` | Pass `-march=native` in Release builds |

Example — lean default build:
```bash
cmake -S scirocco -B scirocco/build-default -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build-default -j
```

Example — full debug build for sensitivity work:
```bash
cmake -S scirocco -B scirocco/build-debug -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON \
      -DSCIROCCO_BUILD_PROBES=ON
```

To build the servers, follow the [Philote dependency instructions](docs/docs/getting-started/installation.md)
and configure with `-DSCIROCCO_BUILD_PHILOTE=ON -DPhiloteCpp_DIR=/path/to/Philote-Cpp/build`.

Install the core library with `cmake --install scirocco/build-default --prefix /your/prefix`.
Consumers can use `find_package(scirocco CONFIG REQUIRED)` and link
`scirocco::scirocco`, with the prefix on `CMAKE_PREFIX_PATH` and Eigen installed.

Operating conditions must be finite, with `0 <= Ma < 1`. Finite negative `cl`
selects alpha mode, and finite nonpositive `Re` selects inviscid mode. Panel
counts must be integers from 4 through 4096 (even for Kulfan). Custom coordinates
must be finite 2-by-N arrays with at least four points, positive chord, and no
repeated adjacent points. Kulfan accepts 3–256 coefficients per surface.
Boundary-layer distributions and wake coordinates are empty in inviscid mode.
Invalid inputs throw exceptions; failed or nonfinite solutions cannot supply gradients.

## Usage

### Python through Philote

The primary Python interface uses
`philote_mdo.openmdao.RemoteExplicitComponent` to connect OpenMDAO to the
compiled Scirocco servers. One server creates parameterized Kulfan/CST
coordinates and the other analyzes them:

```python
import grpc
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent

problem = om.Problem()
problem.model.add_subsystem(
    "geometry",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50052"),
        n_cst_upper=6, n_cst_lower=6, npanel=100,
    ),
    promotes_inputs=["cst_upper", "cst_lower"],
    promotes_outputs=["x"],
)
problem.model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50051"),
        naca="0012", npanel=100,
    ),
    promotes_inputs=["x"],
)
```

Scirocco does not currently ship an in-process `import scirocco` module. See
the [Python quickstart](https://chrislupp.github.io/Scirocco/docs/getting-started/quickstart)
for installation, server launch commands, inputs, and a complete first run.

The C++ library remains available directly for native applications:

### Basic Inviscid Analysis

```cpp
#include <scirocco/scirocco.h>

scirocco::Scirocco foil("0012", 199);  // NACA 0012, 199 panels
foil.SetOper(5.0, -1.0, 0.0, 0.0);  // alpha=5°, no Re (inviscid)
foil.Solve();

std::cout << "Cl = " << foil.GetCl() << std::endl;
std::cout << "Cm = " << foil.GetCm() << std::endl;
```

### Viscous Analysis

```cpp
#include <scirocco/scirocco.h>

scirocco::Scirocco foil("2412", 129);  // NACA 2412, 129 panels
double Re = 1.0e6;
double Ma = 0.0;

foil.SetOper(4.0, -1.0, Re, Ma);  // alpha=4°, Re=1e6
bool converged = foil.Solve();

std::cout << "Cl = " << foil.GetCl() << std::endl;
std::cout << "Cd = " << foil.GetCd() << std::endl;
std::cout << "Cdf = " << foil.GetCdf() << std::endl;  // friction drag
std::cout << "Cdp = " << foil.GetCdp() << std::endl;  // pressure drag
```

### Custom Airfoil Coordinates

```cpp
Eigen::MatrixXd coords(2, n);  // row 0 is x; row 1 is z
// ... fill coordinates ...

scirocco::Scirocco foil(coords, 199);
```

### Geometry Modification

```cpp
scirocco::Scirocco foil("0012", 199);

// Deploy a flap
Eigen::Vector2d hinge(0.75, 0.0);  // hinge at 75% chord
double deflection = 10.0;          // 10 degrees
foil.GeomFlap(hinge, deflection);

// Add camber
Eigen::MatrixXd camber_line(n, 2);
// ... define camber ...
foil.GeomAddCamber(camber_line);
```

### Accessing Results

```cpp
// Scalar results
double cl = foil.GetCl();
double cd = foil.GetCd();
double cm = foil.GetCm();
double alpha = foil.GetAlpha();

// Distributions (along airfoil surface)
const Eigen::VectorXd& cp = foil.GetCp();          // pressure coefficient
const Eigen::VectorXd& cf = foil.GetCf();          // skin friction
const Eigen::VectorXd& theta = foil.GetTheta();    // momentum thickness
const Eigen::VectorXd& dstar = foil.GetDeltaStar(); // displacement thickness

// Geometry
const Eigen::MatrixXd& foil_coords = foil.GetFoilCoords();
const Eigen::MatrixXd& wake_coords = foil.GetWakeCoords();
```

Solver progress output is quiet by default. Advanced/debug callers can enable
the previous verbose behavior through the internal state:

```cpp
foil.GetState().param.verb = 1;
```

## Optimization examples

- [Single-point viscous CST optimization](examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py): minimize drag at a prescribed lift with an area constraint.
- [Multipoint CST optimization](examples/cst_range_optimization/README.rst): minimize equally weighted drag at `Cl = 0.55` and `0.70`, with guards at `0.30` and `0.80`. Includes analytic gradients, independent optimality checks, and geometry and convergence plots.

Run the multipoint example from the repository root:

```bash
python examples/cst_range_optimization/run.py \
    --binary scirocco/build-philote/src/philote/scirocco_server \
    --output pa_output/cst_multipoint
```

## API Reference

### Scirocco Class

| Method | Description |
|--------|-------------|
| `Scirocco(naca, npanel)` | Create airfoil from NACA designation |
| `Scirocco(coords, npanel)` | Create airfoil from coordinates |
| `SetOper(alpha, cl, Re, Ma)` | Set operating conditions |
| `Solve()` | Run analysis, returns convergence status |
| `GetCl(), GetCd(), GetCm()` | Get force coefficients |
| `GetCdf(), GetCdp()` | Get friction and pressure drag |
| `GetCp(), GetCf()` | Get surface distributions |
| `GeomFlap(hinge, angle)` | Deploy trailing edge flap |
| `GeomAddCamber(camber)` | Add camber line modification |
| `GeomDerotate()` | Remove angle of attack from geometry |

### SetOper Parameters

| Parameter | Description |
|-----------|-------------|
| `alpha` | Angle of attack in degrees |
| `cl` | Target lift coefficient (use -1.0 for alpha mode) |
| `Re` | Reynolds number (0 for inviscid) |
| `Ma` | Mach number |

## Project Structure

```
scirocco/
├── include/scirocco/
│   ├── core/           # Solver modules
│   │   ├── geometry.h      # Airfoil generation
│   │   ├── spline.h        # Cubic spline interpolation
│   │   ├── paneling.h      # Panel discretization
│   │   ├── inviscid.h      # Potential flow solver
│   │   ├── boundary_layer.h # BL equations
│   │   ├── transition.h    # Transition prediction
│   │   ├── coupled_solver.h # Viscous solver
│   │   └── forces.h        # Force calculation
│   ├── data/           # Data structures
│   └── util/           # Utilities
├── src/                # Implementation
├── tests/              # Test suite
└── examples/           # Example programs
```

## Testing

The default suite (`SCIROCCO_BUILD_TESTS=ON`, all other test flags off) builds four gtest executables that together cover every `src/*.cpp` file:

| Executable | Purpose |
|---|---|
| `scirocco_tests` | Spline, geometry, and the mfoil primal reference suite |
| `scirocco_forward_ref_tests` | Forward `R_x` / `R_U` partials vs mfoil |
| `scirocco_ping_tests` | Analytic Jacobians vs mfoil-CS references (BL, residuals, BuildGlobSys, WakeSys, forces) |
| `scirocco_sensitivity_ref_tests` | End-to-end adjoint `dJ/d{x, α, Ma, Re}` vs complex-step references |

Run it:
```bash
cmake -S scirocco -B scirocco/build-default -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build-default -j
ctest --test-dir scirocco/build-default --output-on-failure
```

Release-mode runtime is ~10 s for 90 tests. The default suite is measured at **93.4% line / 94.0% function** coverage (see `SCIROCCO_COVERAGE` below).

### Static analysis

Run clang-tidy after configuring a build directory so `compile_commands.json`
is available:

```bash
run-clang-tidy -p scirocco/build-default -header-filter='scirocco/include/' 'scirocco/src/.*\.cpp$'
```

### Sensitivity-debug suite (opt-in)

When the end-to-end `sensitivity_ref` tests regress, turn on `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` to build four chain-piece test executables that isolate which link is broken: `scirocco_sensitivity_tests` (per-panel influence derivs), `scirocco_viscous_chain_tests` (geometric chain), `scirocco_viscous_residual_tests` (residual partials), `scirocco_viscous_output_ping_tests` (force output partials). These exercise the same source files as the end-to-end test, so they don't add coverage — they add diagnostic resolution.

### Probe binaries (opt-in)

`SCIROCCO_BUILD_PROBES=ON` builds 14 standalone `probe_*` diagnostic binaries (no gtest, not part of `ctest`) used for ad-hoc adjoint debugging.

### Coverage measurement

```bash
cmake -S scirocco -B scirocco/build-cov -DCMAKE_BUILD_TYPE=Release -DSCIROCCO_COVERAGE=ON
cmake --build scirocco/build-cov -j
cd scirocco/build-cov && LLVM_PROFILE_FILE="cov-%p.profraw" ctest
xcrun llvm-profdata merge -sparse cov-*.profraw -o cov.profdata
xcrun llvm-cov report \
    scirocco_tests \
    -object scirocco_forward_ref_tests \
    -object scirocco_sensitivity_ref_tests \
    -object scirocco_ping_tests \
    -instr-profile=cov.profdata \
    $(ls ../src/*.cpp)
```

## Solver Details

### Inviscid Solver
- Linear vortex panel method with Kutta condition
- Aerodynamic influence coefficient (AIC) matrix formulation
- Stagnation point tracking

### Viscous Solver
- Integral boundary layer equations
- Thwaites initialization for laminar flow
- e^N transition prediction (default n_crit = 9)
- Coupled Newton iteration with the inviscid solution
- Wake modeling with displacement thickness effects

### Compressibility
- Karman-Tsien correction for subsonic flow
- Sutherland viscosity model

## License

Apache License 2.0 — see [LICENSE](LICENSE) for details.

## Acknowledgments

Scirocco is a C++ port of [mfoil](https://websites.umich.edu/~kfid/codes.html) by Krzysztof J. Fidkowski. The original mfoil is a Python implementation of subsonic airfoil analysis and is licensed under the MIT License. See [NOTICE](NOTICE) for details.
