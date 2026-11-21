# AGENTS.md

Guidance for coding agents working in this repository. Claude Code reads it
through `CLAUDE.md`, which imports this file.

## Local overrides

If `AGENTS.local.md` exists at the repository root, read it as well. It holds
one developer's machine-specific setup (build directories, paths, tools) and
takes precedence over this file. It is gitignored, so it stays on that
machine.

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

### fauXFOIL
`fauxfoil/` holds fauXFOIL, an XFOIL-style command line and Qt plot window
built on Scirocco's public API. It is a separate CMake project that compiles
the sibling `scirocco/` tree (or, with `-DFAUXFOIL_SCIROCCO_SOURCE=`, finds
an installed Scirocco):
```bash
cmake -S fauxfoil -B fauxfoil/build -DCMAKE_BUILD_TYPE=Release
cmake --build fauxfoil/build -j
ctest --test-dir fauxfoil/build --output-on-failure
```
Options: `FAUXFOIL_GRAPHICS=OFF` (opt in for the Qt 6 Widgets plot window;
when ON, configuration fails if Qt is not found, never falling back to
text-only), `FAUXFOIL_BUILD_TESTS=ON`. The Qt tests run without a display
(`QT_QPA_PLATFORM=offscreen`). Nothing builds fauXFOIL unless `fauxfoil/`
is configured; the Scirocco build, Docker images and release artifacts do
not include it.

Homebrew ships it as its own formula. `chrislupp/homebrew-tools` has
`scirocco` (the library) and `fauxfoil`, which builds `fauxfoil/` with Qt
from the same release tarball against the installed `scirocco`. Library
users never get fauXFOIL or Qt. Both formulae carry the same `url` and
`sha256`, and the release job updates them together (see
[Releases](#releases)). Keep fauXFOIL buildable against an installed
Scirocco (`-DFAUXFOIL_SCIROCCO_SOURCE=`); the formula depends on it.

fauXFOIL's version is Scirocco's: CMake reads it from
`project(scirocco VERSION ...)` (or the installed package) and passes it in
as `SCIROCCO_VERSION`, which the banner, polar file headers and plot tags
print. Never hard-code it in `fauxfoil/`; the `fauXFOIL.Version` test
fails if the binary reports anything else. The screenshots in
`docs/static/img/fauxfoil-*.png` show the version too; see
[Releases](#releases) for regenerating them.

### Scirocco Viewer
`viewer/` holds Scirocco Viewer, a Qt 6 Quick (QML) program for looking at
the HDF5 files `scirocco::Hdf5Writer` writes (`SCIROCCO_HDF5_OUTPUT=ON`,
see `scirocco/examples/hdf5_output_example.cpp`). It is a separate CMake
project like fauXFOIL, but reads only the HDF5 file and does not link
Scirocco, so its reader (`viewer/src/hdf5_reader.cpp`) must follow any
change to the writer's group and dataset names. It needs Qt 6.5+ (Core,
Gui, Qml, Quick, QuickControls2, Widgets, Charts) and HDF5 (C, C++):
```bash
cmake -S viewer -B viewer/build -DCMAKE_BUILD_TYPE=Release
cmake --build viewer/build -j
viewer/build/SciroccoViewer.app/Contents/MacOS/SciroccoViewer alpha_sweep.h5
```
On macOS it builds as an app bundle. Nothing builds it unless `viewer/` is
configured; the Scirocco build, Docker images and release artifacts do not
include it.

### Build options
- `SCIROCCO_BUILD_TESTS=ON` (default) — solver, geometry, Jacobian, and adjoint suites; see CONTRIBUTING.md for dependency requirements
- `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` (opt-in) — 3 chain-piece isolators: `scirocco_sensitivity_tests`, `scirocco_viscous_residual_tests`, `scirocco_viscous_output_ping_tests`
- `SCIROCCO_BUILD_PHILOTE=OFF` (default) — Philote-Cpp gRPC discipline interface for multidisciplinary optimization frameworks
- `SCIROCCO_BUILD_PROBES=ON` (opt-in) — 14 diagnostic binaries for ad-hoc adjoint debugging
- `SCIROCCO_BUILD_EXAMPLES=ON` (default) — `quickstart`, `naca_example`, `viscous_example`, `bench`, `docs_figure_data`
- Coverage uses compiler flags; there is no `SCIROCCO_COVERAGE` option.
- `SCIROCCO_NATIVE_ARCH=ON` (default) — `-march=native` in Release builds
- `SCIROCCO_USE_BLAS=ON` (default) — BLAS/LAPACK for the dense coupled solve (Eigen fallback)

Build directories are local and disposable; delete and regenerate if CMake cache or test discovery state gets stale.

## Testing

### Default test suite
```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j
ctest --test-dir scirocco/build --output-on-failure
```
Test counts and runtime depend on features and dependencies. Inspect the suite
with `ctest --test-dir scirocco/build -N`. The CI `Coverage` job uploads to
Codecov.

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

Use the compiler flags and measurement commands in CONTRIBUTING.md. There is no
CMake coverage option.

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

## [0.1.0] - 2026-03-18
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

## Contributing

### Branching and PR workflow

| Branch | Purpose |
|--------|---------|
| `main` | Stable, tagged releases only |
| `develop` | Integration branch — all PRs target this |
| `feature/<name>` | New features |
| `fix/<name>` | Bug fixes |
| `docs/<name>` | Documentation-only changes |

Keep branches focused — one logical change per PR. If a fix requires a
refactor, split them into separate PRs.

Stacked PRs (one PR based on another's branch) land one at a time, bottom
first, and every layer lands in `develop`: merge the bottom PR into
`develop`, retarget the next PR to `develop`, update its branch, and repeat.
Merging a PR into the branch below it leaves its commits outside `develop`;
that is how the gradient speedups from #26–#28 were lost.

Before opening a PR, document any new feature on the documentation site
(see [Documenting features](#documenting-features)).

Label every PR when you open it (`gh pr create --label ...`), matching its
commit category:

| PR | Label |
|----|-------|
| `Docs:` | `documentation` |
| `Add:`, `Refactor:`, `CI:` | `enhancement` |
| `Fix:` | `bug` (and `documentation` when it fixes the docs) |
| Release, `develop` → `main` | `release` and `major`, `minor` or `patch` (see [Releases](#releases)) |

Before marking a PR ready for review:
- All CI checks pass (build + full test suite).
- `clang-format` and `clang-tidy` produce no new warnings.
- `CHANGELOG.md` is updated for any feature addition or bug fix.

Squash-merge is preferred for feature branches; merge commits are used for
release integrations (`develop` → `main`).

### Releases

A release is a PR from `develop` to `main` labelled `release` and one of
`major`, `minor` or `patch`, merged with a merge commit. On merge,
`.github/workflows/release.yml` makes every versioned change: it bumps
`scirocco/CMakeLists.txt`, retitles `CHANGELOG.md`, snapshots the docs
(`scripts/release_docs.py`: `docs/versioned_docs/`,
`docs/versioned_sidebars/`, `docs/versions.json`, `lastVersion`), tags,
and publishes. Its `homebrew` job then points the `scirocco` and
`fauxfoil` formulae in `chrislupp/homebrew-tools` at a stable release
(`scripts/homebrew_formula.py`; `RELEASE_TOKEN` must be able to push there).
Those changes come only from the job: agents never commit a docs version
snapshot or a version bump made by hand.

Prepare a release on a `docs/` branch into `develop`, merged before the
release PR:

1. Run `scripts/fauxfoil_screenshots.sh <bump>` with the release's bump
   label and commit the two PNGs it writes. It builds fauXFOIL with the
   coming version in a temporary worktree; the pictures come from this
   script, never from hand-picked or edited plots.
2. Build the docs as they will look after the release: in a throwaway
   worktree, run `scripts/release_docs.py <version> false`, then
   `npm run build` in `docs/`. Delete the worktree; commit nothing from it.
   Links in `docs/src/pages/` and redirects resolve against the default
   docs version, so a page moved since the last release needs
   `createRedirects` in `docs/docusaurus.config.ts` (a static redirect to
   it fails until the release) and a homepage link that looks the page up.

### Commit messages

```
Category: short imperative summary (≤72 chars)

Optional body explaining why, not what. Wrap at 72 columns.
```

Common categories: `Fix`, `Add`, `Refactor`, `Docs`, `CI`, `Test`.
Put issue references in the body, not the subject line.

### Adding tests

New solver features must include at least one test that would fail without the
change. Derivative correctness must be validated against mfoil complex-step
(CS) — never finite differences. See [Reference data and CS validation](#reference-data-and-cs-validation).

Do not suppress clang-tidy warnings with `NOLINT` without a comment on the
same line explaining why.

### Documenting features

A feature is finished only when the documentation site (`docs/docs/`, a
Docusaurus site) covers it, on the same branch, before the PR is opened:

- **API:** new public symbols, options and defaults go on the matching
  `api/` page (for example, a new `Param` field goes in
  `api/solver-params.md`).
- **Usage:** what the feature does and when to use it goes in `concepts/` or
  `guides/`, with a short example.
- **Validation:** comparisons with data, accuracy and known limits go on a
  docs page, with figures drawn from `docs/src/data/` (see `docs/README.md`).
  A README next to validation scripts carries only the steps to reproduce
  them.

Build the site with `npm run build` in `docs/`; it must finish with no broken
links or anchors.

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
fauxfoil/              # fauXFOIL: XFOIL-style command line on Scirocco (own CMake project)
  src/                 # Session, plots, PostScript, font reader; qt/ holds the plot window
  fonts/               # FauxfoilPlot.ttf (compiled into the program) and make_font.py, its source
  tests/               # Google Test suites (text and offscreen Qt)

viewer/                # Scirocco Viewer: Qt Quick viewer for HDF5 output (own CMake project)
  src/                 # HDF5 reader and video exporter exposed to QML
  qml/                 # Main window, charts and airfoil canvas; components/ holds panels

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
