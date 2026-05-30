# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- Kulfan/CST (Class-Shape Transformation) airfoil geometry discipline (`KulfanDiscipline`) with analytic partials of node coordinates w.r.t. Bernstein coefficients and a built-in second-difference smoothness output for waviness constraints.
- Pseudo-transient continuation (PTC) in the coupled BL Newton solver: `R_V` Jacobian augmented with `(1/τ)·I` on BL rows, switched-evolution-relaxation `τ` ramp in `SolveCoupled`. Stops the undamped Newton iteration from diverging on rough intermediate geometries while leaving the converged R = 0 state unchanged. New `Param` fields: `use_ptc`, `tau_init`, `tau_min`, `tau_max`, `tau_grow_clamp`, `tau_shrink_clamp`.
- Diagonal-ridge fallback in `StoreAdjointFactorization` that retries the LU with `R_V + ε·max|R_V|·I` on first failure, controlled by `Param::adjoint_ridge` (default `1e-10`). Catches the genuine `Rxt_xt → 0` ill-conditioning at the transition row that any faithful mfoil port inherits.
- `RefreshPanelGeometry` helper (in `paneling.cpp`) that rebuilds `panel.s` and `panel.t` from the current `panel.x` via `Spline2DCreate`. Called from `SiroccoDiscipline::Compute` after overwriting `foil.x` from an upstream geometry discipline, so the BL solver doesn't run on new coordinates with stale tangents.
- Graceful adjoint-failure cache in `SiroccoDiscipline::ComputePartials`: caches the last successful `Sensitivity` and reuses it if `ComputeSensitivities()` throws, so the optimizer makes progress through transient adjoint glitches instead of crashing.
- `RefreshPanelGeometry` declaration in `scirocco/core/paneling.h`.
- Three airfoil-optimization examples:
  - `examples/kulfan_inviscid_optimization/` — minimize `Cm²` at target `Cl` (inviscid).
  - `examples/kulfan_viscous_optimization/` — minimize `Cd` at target `Cl` with area floor and CST smoothness constraint.
  - `examples/kulfan_multipoint_optimization/` — weighted `Σ wᵢ·Cdᵢ` across N viscous design points sharing one geometry, with one analysis server per point.
- `scirocco_kulfan_server` binary hosting `KulfanDiscipline` for use with `RemoteExplicitComponent` clients.

### Changed
- `Param::niglob` default raised from 50 to 200 to give the SER `τ` ramp headroom on rough intermediate iterates.
- `SiroccoDiscipline::Compute` now refreshes panel-derived state and clears solver state after each `foil.x` overwrite, so it works correctly when chained downstream of any upstream geometry discipline (not just one whose coords match the construction-time NACA paneling).
- `PhilotePartialsTest` viscous-partial tolerances widened from `1e-12` to `1e-3` / `5e-3` to reflect the panel-refresh path difference vs constructor-time `MakePanels`. Engineering-equivalent; documented inline.
- `boundary_layer.cpp::ResidualTransition` and `transition.cpp::MarchAmplification` cleaned up to match `mfoil.py:2662` and `mfoil.py:2602` exactly (off-by-one in the implicit-xt Newton loop guard; non-convergence warning print).
- `scirocco_viscous_chain_tests` (geometric forward-chain CS validation: `DFoilSDFoilXCumulative`, `DSstagDFoilX`, `DWakeSDFoilX`, `DXiDFoilX`, `DWgapDFoilX`, `DWakeXDFoilX`) promoted from opt-in (`SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON`) into the default test suite. Locks the wake-chain VJP correctness in CI at `1e-10` rel / `2e-11` abs over 60k+ Jacobian entries.
- `TOL_VISC_SHAPE` tightened from `1e-4` to `1e-8` in `sensitivity_reference_test.cpp::ViscousShapeSensitivity`. End-to-end adjoint dcl/dcd/dcm vs mfoil complex-step now agrees to ~3e-11 relative; the loose `1e-4` was a holdover from the era when the wake-row chain had the ~0.216 absolute error at TE nodes.

### Fixed
- Stale "PHASE 1 STATUS: STUB (returns zeros)" docstring on `DWakeXDFoilX` in `viscous_geom_deriv.h`. The function has long been a real 120-line linearized predictor-corrector wired into the wake-side `ueinv` chain (`ComputeRueVJP:263`) and the B2c wake-source-panel chain (`ComputeRueVJP:401`); the docstring now describes the implementation accurately.

## [0.2.0] - 2026-05-29

### Added
- NACA 4-digit airfoil optimization example demonstrating Scirocco integration with optimization frameworks
- Philote-Cpp discipline interface for coupling Scirocco with multidisciplinary optimization tools
- HDF5 iteration history writer for saving optimization/analysis results to HDF5 format
- `SCIROCCO_COVERAGE` CMake option that instruments the library and tests for `llvm-cov` source-based coverage. Default suite measured at 93.4% line / 94.0% function coverage.
- README documentation for all test-suite build flags and the local coverage workflow.

### Changed
- Split the test suite into a lean default (`scirocco_tests`, `scirocco_forward_ref_tests`, `scirocco_ping_tests`, `scirocco_sensitivity_ref_tests`) and an opt-in sensitivity-debug suite. Chain-piece tests (`sensitivity_ping`, `viscous_chain`, `viscous_residual_deriv`, `viscous_output_ping`) and the 14 `probe_*` diagnostic binaries now build only when `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` or `SCIROCCO_BUILD_PROBES=ON`.
- Extended CI: GCC + Clang Release matrix runs the default suite, a separate job exercises the opt-in flags, a Clang coverage job uploads lcov to Codecov, and a Valgrind memcheck job fails on any definite/indirect leak.

## [0.1.0] - 2026-03-18

Initial public release.

### Added
- Inviscid analysis using linear vortex panel method with Kutta condition
- Viscous analysis with coupled boundary layer solver (Newton iteration)
- Laminar-to-turbulent transition prediction (e^N method)
- Karman-Tsien compressibility correction for subsonic flow
- NACA 4-digit and 5-digit airfoil generation
- Custom airfoil coordinate input
- Geometry modification: trailing edge flaps, camber addition, derotation
- Force coefficient output: Cl, Cd, Cm, Cdf, Cdp
- Surface distribution output: Cp, cf, theta, delta*, Hk, Ret, ue
- Boundary layer property functions with analytical Jacobians
- 73 reference validation tests against mfoil
- 8 finite-difference ping tests for Jacobian verification
- CMake build system with Eigen3 dependency
- Apache 2.0 license
[Unreleased]: https://github.com/chrislupp/Scirocco/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/chrislupp/Scirocco/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.1.0
