# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Release readiness
- Break Philote RPC service ownership cycles so server shutdown and startup
  failures release the discipline and its solver state.
- Fix stack overwrites in reference-test scalar readers by converting HDF5
  storage values to the destination C++ memory type. Initialize wake coordinates
  and wake-inclusive velocities consistently in debug test fixtures.
- Add ASan/UBSan CI and an opt-in instrumented CMake build, including installed
  consumers and Philote tests.
- Verify published source archives and full-version container images after
  release, with manual rechecks of existing tags.
- Correct test-suite and coverage documentation and record main branch
  protection requirements and the release-token administrator exception.
- Reject empty, malformed, nonfinite, or unordered camber data before changing
  geometry, and handle server startup failures without dereferencing null.
- Validate the prepared release commit, installed consumer, Python references,
  documentation, and both container images before publication. Publish the tested
  images without rebuilding; refuse publication if `main` advances.
- Reject nonfinite operating conditions, invalid panel counts, malformed coordinate
  options, and nonfinite solver results instead of reporting successful analyses.
- Clear previous boundary-layer and wake results when switching to inviscid mode.
- Default the optional Philote interface to OFF; server CI and Docker builds
  explicitly enable it.
- Keep prerelease target versions fixed through candidate progression and stable
  promotion; validate release labels and tags before publishing atomically.

### Added
- `quickstart`, `bench`, and `docs_figure_data` example programs. `bench`
  times paneling, inviscid and viscous solves, and the full gradient at
  several panel counts; `docs_figure_data` writes the data plotted on the
  documentation site.
- Redesigned documentation site: homepage, C++ quickstart and guides,
  performance page, local search, and redirects from moved pages.
- Automatically snapshot Docusaurus documentation for stable releases, select
  the release as the default docs version, and validate it before publication.
- Relocatable CMake package configuration and an installed-consumer smoke test.
- Solver input, mode-switch, RPC recovery, and release-version regression tests.
- Docker multi-stage images for `scirocco-server` (port 50051) and
  `kulfan-server` (port 50052), both running on `debian:stable-slim` as
  a non-root user. A CI workflow validates both builds on every PR/push
  and publishes to Docker Hub as `chrislupp/scirocco` and
  `chrislupp/scirocco-kulfan` through the validated release workflow.
- Docusaurus 3 documentation site scaffolded at `docs/` with KaTeX math
  support, versioned docs, and a GitHub Actions workflow that deploys to
  GitHub Pages on push to `develop`.
- Kulfan/CST (Class-Shape Transformation) airfoil geometry discipline
  (`KulfanDiscipline`) with analytic partials of node coordinates w.r.t.
  Bernstein coefficients and a built-in second-difference smoothness
  output for waviness constraints.
- Pseudo-transient continuation (PTC) in the coupled BL Newton solver:
  `R_V` Jacobian augmented with `(1/τ)·I` on BL rows,
  switched-evolution-relaxation `τ` ramp in `SolveCoupled`. Stops the
  undamped Newton iteration from diverging on rough intermediate
  geometries while leaving the converged R = 0 state unchanged. New
  `Param` fields: `use_ptc`, `tau_init`, `tau_min`, `tau_max`,
  `tau_grow_clamp`, `tau_shrink_clamp`.
- `RefreshPanelGeometry` helper (in `paneling.cpp`) that rebuilds
  `panel.s` and `panel.t` from the current `panel.x` via
  `Spline2DCreate`. Called from `SiroccoDiscipline::Compute` after
  overwriting `foil.x` from an upstream geometry discipline, so the BL
  solver doesn't run on new coordinates with stale tangents.
- `RefreshPanelGeometry` declaration in `scirocco/core/paneling.h`.
- Three airfoil-optimization examples:
  - `examples/kulfan_inviscid_optimization/` — minimize `Cm²` at target
    `Cl` (inviscid).
  - `examples/kulfan_viscous_optimization/` — minimize `Cd` at target
    `Cl` with area floor and CST smoothness constraint.
  - `examples/kulfan_multipoint_optimization/` — weighted `Σ wᵢ·Cdᵢ`
    across N viscous design points sharing one geometry, with one
    analysis server per point.
- `scirocco_kulfan_server` binary hosting `KulfanDiscipline` for use
  with `RemoteExplicitComponent` clients.
- Runnable Philote 0.8 Python tutorial suite under `examples/tutorials/`,
  covering inviscid and viscous analysis, geometry changes, custom
  coordinates, analytic sensitivities, and HDF5 histories.

### Changed
- Viscous `ComputeSensitivities()` is about 13x faster (11.2 s to 0.88 s at
  199 panels on an M1 Pro). The geometric Jacobians in the `R_ue` chain,
  `d(gamma)/d(foil.x)` and the wake sweeps, were rebuilt 16 times per gradient;
  they are now built once in a `RueVJPContext` and shared by all four outputs.
  Results are bitwise identical.
- `DGammaDFoilX` assembles the AIC residual Jacobian once
  (`AICResidualJacobianX`) and does one LU solve, instead of one O(N^2) VJP per
  row. The viscous gradient drops from 0.88 s to 0.21 s at 199 panels.
  `AICResidualSensVJP` and `AICMatVecVJP` now share this kernel
  (`AICMatVecJacobianX`) instead of two copies of the panel chain. Results
  change at round-off level only (about 1e-14 of the largest entry).
- `ComputeRueVJP` accumulates adjoints on wake points and applies the wake
  Jacobians once per call, instead of a dense row update for every foil-point,
  wake-panel pair. Terms that do not depend on the adjoint vector move into
  `RueVJPContext`. The viscous gradient drops from 0.21 s to about 0.08 s at
  199 panels, about 140x faster than before these changes.
- `RueVJPContext` assembles the source-panel rows of `ComputeRueVJP`
  (`dsrc_dx`, plus the wake-source terms in `wake_rows`) once, so the panel
  derivatives are evaluated once per gradient rather than once per output. Each
  output's VJP is now a few matrix-vector products (0.1 ms). The viscous
  gradient drops from 0.08 s to 0.057 s at 199 panels.
- The wake-geometry and gamma chains of the viscous adjoint run in reverse
  mode. `LinearizeWakeSweep` and `WakeSweepVJP` propagate the four outputs back
  through `BuildWake`'s predictor-corrector sweep, replacing the forward wake
  Jacobians (which carried all 2N coordinate directions) and the 2N-column
  `d(gamma)/d(foil.x)` solve. `ComputeRueVJPBatch` takes several adjoint
  vectors at once. The viscous gradient drops from 57 ms to 36 ms at 199 panels
  and grows more slowly with panel count. The forward Jacobians remain for
  validation.
- Remove the unreleased 0.2.0 documentation scaffold so its first stable release
  can create a complete snapshot instead of colliding with placeholder pages.
- Temporarily disable the CI coverage job and Codecov upload.
- Multipoint Kulfan optimization example switched from OpenMDAO's
  `ScipyOptimizeDriver`/SLSQP to IPOPT via `modopt.OpenMDAOProblem`.
  gRPC server ports shifted to 50060–50063 to avoid collisions. With
  Philote 0.8, the example initializes each angle of attack at its lift
  target and rejects unsuccessful solves. Verified: `Solve_Succeeded`
  at iter 222, Σw·Cd = 4.601e-3.
- Philote tutorial, Docker, and CI builds now pin the protocol-0.8 C++
  implementation and `philote-mdo==0.8.0`; CI executes the complete
  Python tutorial sequence against a Release build.
- All C++ headers converted to Doxygen `/** */` style with
  `@brief`/`@param`/`@return` tags; `#include` guards replaced with
  `#pragma once` across all 32 headers. Internal `_dMa` closure
  derivatives and `ThermoMaDeriv` merged from the separate
  `bl_chain_internals.h` into `bl_properties.h`.
- `Param::niglob` default raised from 50 to 200 to give the SER `τ`
  ramp headroom on rough intermediate iterates.
- `SiroccoDiscipline::Compute` now refreshes panel-derived state and
  clears solver state after each `foil.x` overwrite, so it works
  correctly when chained downstream of any upstream geometry discipline
  (not just one whose coords match the construction-time NACA paneling).
- `boundary_layer.cpp::ResidualTransition` and
  `transition.cpp::MarchAmplification` cleaned up to match
  `mfoil.py:2662` and `mfoil.py:2602` exactly (off-by-one in the
  implicit-xt Newton loop guard; non-convergence warning print).
- `scirocco_viscous_chain_tests` (geometric forward-chain CS validation:
  `DFoilSDFoilXCumulative`, `DSstagDFoilX`, `DWakeSDFoilX`,
  `DXiDFoilX`, `DWgapDFoilX`, `DWakeXDFoilX`) promoted from opt-in
  (`SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON`) into the default test
  suite. Locks the wake-chain VJP correctness in CI at `1e-10` rel /
  `2e-11` abs over 60k+ Jacobian entries.
- `TOL_VISC_SHAPE` tightened from `1e-4` to `1e-8` in
  `sensitivity_reference_test.cpp::ViscousShapeSensitivity`. End-to-end
  adjoint dcl/dcd/dcm vs mfoil complex-step now agrees to ~3e-11
  relative.

### Fixed
- `hdf5_output_example` now records non-converged sweep points without
  requesting unavailable sensitivities, so the sweep continues and the
  output file retains every attempted operating condition.
- Preserve HDF5 reference-test discovery when the output writer has already
  created its library target.
- Use the multipoint example's NACA 0012 CST fit and finite trailing-edge gap
  throughout the Philote tutorials so their viscous sweeps converge. Add a
  Reynolds-sweep regression and refresh the documented results.
- Handle failed continuation and trust-region trial exceptions explicitly,
  preserving state rollback while satisfying Clang-Tidy.
- Normalize the Python reference pressure forces and lift partials once after
  panel integration, fixing nonunit-chord outputs and lift trim. Remove manual
  force normalization from the gradient reference generators.
- Remove the remaining distant-source streamfunction cancellation in Python
  and C++, including endpoint, evaluation-point, and panel-frame derivatives.
  Add independent source-integral checks and an 80-digit full Python gradient
  reference with a tighter regression for the disputed nonunit-chord case.
- Avoid cancellation in distant linear-vortex streamfunction coefficients and
  coordinate derivatives by integrating the same log-distance kernel directly
  in Python and C++. Validate against independent 60-digit integrals and remove
  the nonunit-chord gradient tolerance exception.
- Keep wake-spacing derivative iterations consistent with the primal, refresh
  stagnation/transition data through convergence, and retain nonzero sparse
  Jacobian entries.
- Constructor-time arc labels now have an explicit fixed mode, so standalone
  shape gradients match Python's unchanged arc labels. Spline and cumulative
  geometry modes retain their matching derivative chains.
- Lift-target gradients now include the implicit angle adjustment and target
  derivatives, with the original Python frozen-wake behavior differentiated
  through its initial angle. Inviscid lift targets are now actually solved.
- Near-field pressure-drag gradients are computed for viscous flow; separate
  angle, Mach, and Reynolds derivatives are available for `cdpi`.
- The Python CS update preserves upstream post-update projections and wake
  rebuilding. Its real update is checked directly against bundled Python.
- Cached gradient getters reject stale results after a solve, input change,
  or mutable state access.
- Viscous shape gradients now differentiate Python mfoil's spline arc
  lengths throughout the boundary-layer and wake chains, including
  leading-edge and trailing-edge nodes. Explicit cumulative-length callers
  retain their own matching derivative convention.
- Python complex-step references preserve the trailing-edge absolute-value
  derivative and support complex knots in the original `spline2d` routine.
  Independent spline identities and real-primal comparisons validate the
  reference extension before it is used to check C++ gradients.
- Inviscid `cd` angle and Mach derivatives are zero, matching the Python
  output, instead of incorrectly returning derivatives of `cdpi`.
- Unconverged analyses, stale operating conditions, and singular/inaccurate
  adjoint solves now report errors.
  Removed `Param::adjoint_ridge` and the Philote stale/zero-gradient fallback.
- Lift-target iteration refreshes force linearizations as Python mfoil does
  and applies the correct lower-surface velocity sign. Failed linear solves
  no longer access an unallocated update matrix.
- Gradient regression tests now use fresh Python references for the actual
  spline geometry input path, including full-shape directions and leading-
  and trailing-edge nodes, without relaxing tolerances or skipping nodes.
- Stale "PHASE 1 STATUS: STUB (returns zeros)" docstring on
  `DWakeXDFoilX` in `viscous_geom_deriv.h`. The function is a real
  120-line linearized predictor-corrector wired into the wake-side
  `ueinv` chain (`ComputeRueVJP:263`) and the B2c wake-source-panel
  chain (`ComputeRueVJP:401`); the docstring now describes the
  implementation accurately.

## [0.2.0] - 2026-05-29

### Added
- NACA 4-digit airfoil optimization example demonstrating Scirocco
  integration with optimization frameworks.
- Philote-Cpp discipline interface for coupling Scirocco with
  multidisciplinary optimization tools.
- HDF5 iteration history writer for saving optimization/analysis results
  to HDF5 format.
- `SCIROCCO_COVERAGE` CMake option that instruments the library and
  tests for `llvm-cov` source-based coverage. Default suite measured at
  93.4% line / 94.0% function coverage.
- README documentation for all test-suite build flags and the local
  coverage workflow.

### Changed
- Split the test suite into a lean default (`scirocco_tests`,
  `scirocco_forward_ref_tests`, `scirocco_ping_tests`,
  `scirocco_sensitivity_ref_tests`) and an opt-in sensitivity-debug
  suite. Chain-piece tests (`sensitivity_ping`, `viscous_chain`,
  `viscous_residual_deriv`, `viscous_output_ping`) and the 14 `probe_*`
  diagnostic binaries now build only when
  `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` or
  `SCIROCCO_BUILD_PROBES=ON`.
- Extended CI: GCC + Clang Release matrix runs the default suite, a
  separate job exercises the opt-in flags, a Clang coverage job uploads
  lcov to Codecov, and a Valgrind memcheck job fails on any
  definite/indirect leak.

## [0.1.0] - 2026-03-18

Initial public release.

### Added
- Inviscid analysis using linear vortex panel method with Kutta
  condition.
- Viscous analysis with coupled boundary layer solver (Newton
  iteration).
- Laminar-to-turbulent transition prediction (e^N method).
- Karman-Tsien compressibility correction for subsonic flow.
- NACA 4-digit and 5-digit airfoil generation.
- Custom airfoil coordinate input.
- Geometry modification: trailing edge flaps, camber addition,
  derotation.
- Force coefficient output: Cl, Cd, Cm, Cdf, Cdp.
- Surface distribution output: Cp, cf, theta, delta*, Hk, Ret, ue.
- Boundary layer property functions with analytical Jacobians.
- 73 reference validation tests against mfoil.
- 8 finite-difference ping tests for Jacobian verification.
- CMake build system with Eigen3 dependency.
- Apache 2.0 license.

[Unreleased]: https://github.com/chrislupp/Scirocco/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/chrislupp/Scirocco/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.1.0
