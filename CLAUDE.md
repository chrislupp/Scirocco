# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Scirocco** is a modern C++ library for inviscid and viscous aerodynamic analysis of airfoils using panel methods and integral boundary layer theory. The codebase includes a primal solver, complete adjoint chain (derivatives with respect to geometry, angle of attack, Reynolds number, Mach number), extensive validation against complex-step reference data, and optional Philote integration for multidisciplinary optimization.

## Build Commands

The project uses CMake 3.18+ with a single source tree (`scirocco/`) but flexible build directories.

### Basic builds
```bash
# Default Release build (primal solver + default tests + examples)
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j

# Debug build (includes sensitivity debug tests)
cmake -S scirocco -B scirocco/build-debug -DCMAKE_BUILD_TYPE=Debug \
      -DSCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON \
      -DSCIROCCO_BUILD_PROBES=ON
cmake --build scirocco/build-debug -j
```

### Build options
- `SCIROCCO_BUILD_TESTS=ON` (default) — 5 gtest executables: `scirocco_tests`, `scirocco_forward_ref_tests`, `scirocco_ping_tests`, `scirocco_sensitivity_ref_tests`, `scirocco_viscous_chain_tests`
- `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` (opt-in) — 3 chain-piece isolators: `scirocco_sensitivity_tests`, `scirocco_viscous_residual_tests`, `scirocco_viscous_output_ping_tests`
- `SCIROCCO_BUILD_PHILOTE=ON` (default) — Philote-Cpp gRPC discipline interface for multidisciplinary optimization frameworks
- `SCIROCCO_BUILD_PROBES=ON` (opt-in) — 14 diagnostic binaries for ad-hoc adjoint debugging
- `SCIROCCO_BUILD_EXAMPLES=ON` (default) — `naca_example`, `viscous_example`
- `SCIROCCO_COVERAGE=OFF` — `-fprofile-instr-generate -fcoverage-mapping` for `llvm-cov` reporting
- `SCIROCCO_NATIVE_ARCH=ON` (default) — `-march=native` in Release builds

Build directories are local and disposable; delete and regenerate if CMake cache or test discovery state gets stale.

## Testing

### Default test suite
```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j
ctest --test-dir scirocco/build --output-on-failure
```
Runtime ~10 s for 90 tests. Achieves **93.4% line / 94.0% function** coverage (see coverage command below).

### Test categories

| Test Suite | Purpose | Required Dependencies |
|---|---|---|
| `scirocco_tests` | Spline, geometry, primal reference vs mfoil | HDF5 |
| `scirocco_forward_ref_tests` | Forward `R_x` / `R_U` partials vs mfoil complex-step | HDF5 |
| `scirocco_ping_tests` | Analytic Jacobians (BL, residuals, BuildGlobSys, WakeSys, forces) vs mfoil-CS | HDF5 |
| `scirocco_sensitivity_ref_tests` | End-to-end adjoint `dJ/d{x, α, Ma, Re}` vs CS references | HDF5 |
| `scirocco_viscous_chain_tests` | Closed-form forward geometric chains validated against CS | HDF5 |
| `scirocco_sensitivity_tests` (debug) | Per-panel influence derivatives | HDF5 |
| `scirocco_viscous_residual_tests` (debug) | Residual partial derivatives (phase 2) | HDF5 |
| `scirocco_viscous_output_ping_tests` (debug) | Force output partials (phase 3) | HDF5 |

### Run a single test
```bash
# Run one executable
scirocco/build/scirocco_ping_tests

# Run one test case within an executable
scirocco/build/scirocco_sensitivity_ref_tests --gtest_filter="SensitivityRefTest.AdjointAlpha"

# List all test cases in an executable
scirocco/build/scirocco_sensitivity_ref_tests --gtest_list_tests
```

### Static analysis
```bash
run-clang-tidy -p scirocco/build -header-filter='scirocco/include/' 'scirocco/src/.*\.cpp$'
```

### Coverage measurement

Configure with `SCIROCCO_COVERAGE=ON`, run ctest with `LLVM_PROFILE_FILE="cov-%p.profraw"`, then merge and report with `llvm-profdata` / `llvm-cov`. Platform-specific invocation (e.g. `xcrun` on macOS) goes in `CLAUDE.local.md`.

## Code Architecture

### Module organization

The library is split into **4 conceptual layers** across `src/` and `include/scirocco/`:

1. **Geometry & Paneling** — airfoil generation, spline interpolation, panel discretization
   - `geometry.{h,cpp}` — NACA 4/5-digit, flap deployment, camber addition, derotation
   - `spline.{h,cpp}` — cubic spline with clamped/natural BC
   - `paneling.{h,cpp}` — panel node layout and topology

2. **Inviscid Solver** — linear vortex panel method with Kutta condition
   - `inviscid.{h,cpp}` — AIC matrix, influence coefficient assembly, stagnation point tracking
   - `panel_influence.{h,cpp}` — aerodynamic influence coefficients (primal)
   - `panel_influence_deriv.{h,cpp}` — AIC derivatives w.r.t. geometry
   - `wake.{h,cpp}` — wake closure and displacement thickness feedback

3. **Boundary Layer Solver** — integral BL equations with e^N transition
   - `boundary_layer.{h,cpp}` — Thwaites laminar initialization, e^N predictor, closure equations
   - `transition.{h,cpp}` — transition detection and criteria
   - `bl_properties.{h,cpp}` — BL property closures (momentum thickness, shape factor)
   - `coupled_solver.{h,cpp}` — Newton coupling between inviscid and viscous

4. **Adjoint & Sensitivity** — reverse-mode automatic differentiation
   - `adjoint.{h,cpp}` — reverse-mode accumulation for all states
   - `inviscid_deriv.{h,cpp}` — inviscid system partials
   - `viscous_geom_deriv.{h,cpp}` — geometric chain (d(foil.s)/d(foil.x), wake closure, etc.)
   - `viscous_residual_deriv.{h,cpp}` — residual partials (phase 2 of adjoint)
   - `forces.{h,cpp}` — force/moment calculation (used in primal and adjoint)

### Reference test data

Sensitivity validation uses HDF5 reference files generated by mfoil with complex-step (CS) finite differences:
- `mfoil_reference.h5` — primal outputs and reference solutions
- `forward_reference.h5` — forward partials `R_x` / `R_U`
- `ping_reference.h5` — individual Jacobian terms (BL, residuals, BuildGlobSys, WakeSys, forces)
- `sensitivity_reference.h5` — end-to-end sensitivities `dJ/d{x, α, Ma, Re}`
- `geom_chain_reference.h5` — forward geometric chain terms
- `comprehensive_reference.h5` — comprehensive viscous residual partials
- `viscous_output_ping_reference.h5` — force output partials

Located in `scirocco/tests/`, used by CMake `target_compile_definitions()` to bake file paths at build time.

### Solver state & parameters

All solver state is encapsulated in `scirocco::State` (defined in `data/state.h`):
- Operating conditions: `alpha`, `cl_target`, `Re`, `Ma`
- Geometry: foil coordinates, paneling, wake
- Primal: `cp`, `cf`, `theta`, `dstar`, inviscid/viscous outputs
- Adjoint: lambda vectors (reverse-mode multipliers), sensitivity seeds
- Parameters: tolerances, iteration limits, verbosity (`param.verb`)

The main entry point `Scirocco` class (`scirocco.h`, `scirocco.cpp`) wraps `State` and orchestrates solve phases.

### Adjoint structure

The adjoint is organized in **three phases** (reflected in test executables and probe binaries):

1. **Phase 1: Influence derivatives** — d(AIC, cp, circ)/d(foil.x) via `panel_influence_deriv.cpp`
2. **Phase 2: Residual partials** — d(residual)/d(state) via `viscous_residual_deriv.cpp` (boundary layer equations)
3. **Phase 3: Output partials** — d(Cl, Cd, Cm)/d(state) via `forces.cpp`

Reverse-mode accumulation walks backward from output adjoint seeds through these phases via lambdas in `adjoint.cpp`. Key insight: **forward geometric chain** (`viscous_geom_deriv.cpp`) feeds d(foil.s)/d(foil.x) to phase 1, which cascades through the rest.

### Reference data and CS validation

The canonical reference for all derivative correctness is **mfoil complex-step (CS)**. Never use finite differences (FD) as the reference — CS is exact to machine precision. If an adjoint test fails, verify the primal against mfoil first, then compare each phase against its CS reference before modifying Scirocco code.

## Changelog

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/) format and **must be updated with every feature addition or bug fix** before committing.

### Structure

```
## [Unreleased]        ← all in-progress work goes here
### Added              ← new capabilities
### Changed            ← behaviour changes to existing features
### Fixed              ← bug fixes

## [0.2.0] - 2026-05-29
...
```

Only the `[Unreleased]` section is edited during development. On release, it is retitled with the version number and date, a new empty `[Unreleased]` block is opened, and the comparison links at the bottom of the file are updated.

### Entry style

- One bullet per logical change, wrapped at 80 columns.
- Lead with the affected symbol or subsystem in backticks when helpful.
- Be specific enough that a user can understand the impact without reading the diff.
- Omit entries for pure documentation, formatting, or test-only changes unless they affect observable behaviour.

### Sections to use

| Section | Use for |
|---------|---------|
| `Added` | New public API, new CMake options, new binaries, new examples |
| `Changed` | Modified behaviour, renamed symbols, tightened tolerances, default value changes |
| `Fixed` | Incorrect results, crashes, stale documentation that misrepresented behaviour |

## Development Notes

### Compilers & dependencies

- **C++17** required; tested on Clang 14+ and GCC 9+
- **Eigen 3.4+** for linear algebra (dense matrices/vectors)
- **Google Test** (optional) for test discovery and assertions
- **HDF5** (optional) for reference data I/O in validation tests
- **Protobuf + gRPC + Philote-Cpp** (optional) for multidisciplinary optimization

### Code style

- **#pragma once** for header guards (not `#ifndef`)
- **Doxygen-style comments** for public APIs (/// or /** */)
- **Clang-format** configured in `.clang-format` for consistent indentation
- **Clang-tidy** configured in `.clang-tidy` for static analysis

Run both before every commit and push:

```bash
# Format all changed source files
clang-format -i scirocco/src/*.cpp scirocco/include/scirocco/**/*.h

# Static analysis (requires a compile_commands.json from CMake)
run-clang-tidy -p scirocco/build -header-filter='scirocco/include/' 'scirocco/src/.*\.cpp$'
```

### Workflow for adjoint regressions

When `scirocco_sensitivity_ref_tests` fails:

1. **Narrow down**: Build `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON`
2. **Isolate**: Run `scirocco_sensitivity_tests`, `scirocco_viscous_chain_tests`, `scirocco_viscous_residual_tests`, `scirocco_viscous_output_ping_tests` to find which phase broke
3. **Debug**: Use probe binaries (`SCIROCCO_BUILD_PROBES=ON`) for ad-hoc inspection (e.g., `probe_rv_compare`, `probe_drv_dma_cs`)
4. **Validate**: Compare primal solves against mfoil before assuming adjoint error

## Directory Structure

```
scirocco/
  include/scirocco/
    core/              # Solver modules (geometry, paneling, inviscid, BL, adjoint)
    data/              # State structures and enum types
    properties/        # BL property closures (separate module)
    io/                # Optional I/O utilities
    util/              # Helper utilities
    *.h                # Public API headers (scirocco.h, config.h)
  src/                 # Implementation files
  tests/               # Google Test suite + HDF5 reference data
  examples/            # Example programs (naca_example, viscous_example)
  CMakeLists.txt       # Primary CMake configuration
  docs/                # Docusaurus documentation site (separate build)

examples/              # Standalone example applications

docs/                  # Docusaurus site (Node.js build, separate from scirocco/)

plans/                 # Developer documentation and design notes
  adjoint.md           # Adjoint formulation and chain structure

mfoil/                 # Reference mfoil Python implementation (for comparison)

Testing/               # Validation data and reference generation scripts

validation/            # Post-run validation and diagnostics
```

@CLAUDE.local.md
