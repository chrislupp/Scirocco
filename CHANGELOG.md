# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- Homebrew formulae: `brew install chrislupp/tools/scirocco` installs the
  core library, headers and CMake package from the release source, and
  `brew install chrislupp/tools/fauxfoil` installs `fauXFOIL` with its Qt
  plot window. The release workflow's new `homebrew` job points both
  formulae in `chrislupp/homebrew-tools` at each stable release.

## [0.3.1] - 2026-11-21

### Fixed
- `fauXFOIL` Cp plot: the wake displacement edges on the airfoil sketch now
  start at the trailing-edge displacement points and share the wake delta*
  between the sides as XFOIL does, instead of splitting it evenly about the
  wake centreline, which left a visible break at the trailing edge. Where
  the boundary layer is thick (separated flow at high alpha), wake edge
  points that would fold back behind the trailing edge are left out.

## [0.3.0] - 2026-11-21

### Added
- `fauXFOIL`: an interactive command-line program that looks and behaves like
  XFOIL 6.99 (menus, prompts, abbreviations, per-iteration printout, polar
  and coordinate file formats) and runs Scirocco. XFOIL batch scripts run
  unchanged. It is a separate CMake project in `fauxfoil/` that builds
  against the Scirocco source tree beside it or an installed Scirocco.
- `fauXFOIL` plots: with Qt 6, a native plot window (no X11 needed on
  macOS) shows XFOIL's plots: Cp(x) after each point and for sequences,
  pressure vectors (`CPV`), boundary-layer variables (`VPLO`), polars
  (`PPLO`) and the `GDES` geometry. `HARD` writes PostScript hardcopies to
  `plot.ps`. Lettering uses FauxfoilPlot, a plotter-style font compiled
  into the program; hardcopies draw its letters as outlines, so they need
  no fonts. The window is opt-in with `FAUXFOIL_GRAPHICS=ON`, which
  requires Qt 6; otherwise `fauXFOIL` is text-only.
- `fauXFOIL` cursor input and the rest of XFOIL's plot commands: zoom
  (`Z`, `U`), plot-object size (`SIZE`), annotation (`ANNO`), boundary-layer
  velocity profiles (`BL`, `BLC`, `BLWT`), reference Cp, force and polar
  overlays (`CREF`, `FREF`, `RGET`, `RDEL`), polar axis limits (`PPAX`),
  polar airfoils (`APLO`), `VPLO` blowup, symbols, labels, clipping, dump
  and overlay, and the `GDES` blowup, window, tick, inch-scale, overlay,
  curvature and distance commands.
- `Glob::iteration_observer`: optional read-only callback that the coupled
  Newton solve calls at every iteration, and `Glob::relax`, the
  under-relaxation factor of the last update. Neither changes the solve.

## [0.2.0] - 2026-11-15

### Added
- Docker images: each push to `develop` publishes `chrislupp/scirocco:next`
  and `chrislupp/scirocco-kulfan:next` to Docker Hub, once the image passes
  the smoke test. Release images keep their version and `latest` tags.
- `SCIROCCO_USE_BLAS` CMake option (default on): the dense product and LU of
  the coupled solve's Schur complement go through BLAS/LAPACK (Accelerate,
  OpenBLAS, MKL) when CMake finds one, Eigen otherwise. Viscous solves are
  10-15% faster at 200 panels and 15% at 284 on an M1 Pro.
- `Param::warm_start` (default on): `Scirocco::Solve` starts a viscous solve
  from the last converged boundary layer of the same object, as XFOIL does,
  and falls back to a cold start (after at most 30 iterations) if that does
  not converge. The previous boundary layer is shifted along each surface
  to the new stagnation point, as XFOIL's STMOVE does, and its edge
  velocity is predicted from the new inviscid solution at the old mass
  defect, so a 0.5-degree polar step mostly takes 5-6 Newton iterations.
  This also applies across geometry updates through the Philote discipline.
- `BoundaryLayerMarch` (new header `core/bl_march.h`) and
  `Param::bl_march` (default on): before each coupled Newton step the
  boundary layer is re-solved by a mixed-mode march at the current ue, as
  XFOIL's MRCHDU does. Each station holds a ue-Hk combination fixed (direct
  mode in attached flow, inverse near separation), transition is placed with
  `UpdateTransition`'s own test, the stagnation pair and the first wake
  station are re-solved too, and failed stations are extrapolated. Over 144
  cases (NACA 0012/2412/4412, 99-284 panels, alpha -2 to 12, Re 5e5 and
  3e6) the median Newton iteration count drops from 19 to 10 and the mean
  from 47 to 15, and 142 cases converge without continuation (128 before);
  converged solutions are unchanged. On the XFOIL comparison cases Scirocco
  now takes 8-11 steps, like XFOIL.
- `Param::step_limiter` and `StepLimiter::Relative`: XFOIL's
  under-relaxation rule as an alternative to the default mfoil rule.
- `Param::rtol_floor` (1e-8): the coupled Newton solve also stops once the
  residual is below this and a step no longer halves it, instead of running
  to `niglob` at the round-off floor.
- `Param::transition` and `TransitionModel::Captured` (experimental): the
  amplification factor and shear stress live at every node and the laminar
  and turbulent closures are blended smoothly over a transition length
  (`capt_L`), so no station switches equations when transition moves. On
  NACA 2412 at Re 1e6 the adjoint dCd/dalpha ripple drops from 25% to 1.6%,
  and it converges under panel refinement. Against measured Eppler 387 drag
  (UIUC Vol 5, Re 1e5-4.6e5) the RMS error is 5.3%, against 16% for the cut
  cell; `capt_L` = 60 momentum thicknesses is calibrated on that data. The
  primal matches the mfoil prototype
  (`mfoil/mfoil_captured.py`) to 1e-10. The adjoint is fully analytic: shape
  sensitivities match its complex step to 1e-8, and alpha, Re and Mach to
  2e-11. Prescribed angle of attack only (no lift-target mode).
- `Param::onset` and `OnsetModel::Envelope`: opt-in amplification onset
  ramp for the e^N model, fitted to spatial Orr-Sommerfeld envelope N-factor
  curves of Falkner-Skan profiles (H = 2.36-3.48). Growth starts at the
  critical Re_theta and keeps the default asymptote, so far-downstream
  amplification is unchanged. Shape gradients near the onset become two to
  four times smaller and lose the node-to-node oscillation. The default stays
  `OnsetModel::Cubic` (XFOIL/mfoil). Adjoint derivatives are validated
  against mfoil complex step (`param.onset = 'envelope'`); the stability
  solver, data, and fit are in `scirocco/validation/onset_envelope`.
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
  solver doesn't run on new coordinates with stale tangents. Declared in
  `scirocco/core/paneling.h`.
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
- Runnable Philote-MDO 0.8.0 Python tutorial suite under `examples/tutorials/`,
  covering inviscid and viscous analysis, geometry changes, custom
  coordinates, analytic sensitivities, and HDF5 histories.
- Optional Philote-Cpp discipline interface for coupling Scirocco with
  multidisciplinary optimization tools (`SCIROCCO_BUILD_PHILOTE`, OFF by
  default; server CI and Docker builds enable it).
- HDF5 iteration history writer (`SCIROCCO_HDF5_OUTPUT`) for saving
  analysis and optimization results.
- `examples/naca_optimization/`: NACA 4-digit airfoil optimization through
  an optimization framework.
- `examples/cst_range_optimization/`: CST multipoint optimization for low
  drag at Cl = 0.55 and 0.70, with guard conditions at Cl = 0.30 and 0.80.
- `SCIROCCO_SANITIZERS` CMake option (AddressSanitizer and
  UndefinedBehaviorSanitizer) and a sanitizer CI job, including installed
  consumers and Philote tests.
- Release validation: the release workflow checks the prepared release
  commit, installed consumer, Python references, documentation, and both
  container images before publishing; publishes the tested images without
  rebuilding; and refuses to publish if `main` advances. Published source
  archives and full-version images are verified after release, and existing
  tags can be rechecked manually.

### Changed
- `ResidualStation`, `ResidualTransition` and `BuildGlobSys` take an optional
  `scalar_derivatives` flag (default true). Newton iterations and the
  boundary-layer march no longer compute the Re and Mach residual partials,
  which only the adjoint uses; the solver assembles them once at the accepted
  solution. Viscous solves are about 15% faster (NACA 0012, 200 panels:
  28-32 ms cold, 12-15 ms for a 0.5-degree polar step).
- The influence-matrix LU from the inviscid solve is stored
  (`Isol::AIC_lu`, `SolveAIC`) and reused by `CalcUeM` and the adjoint
  instead of refactoring, and `BuildGlobSys` adds the stagnation-point
  chain as triplets instead of inserting into the compressed matrix.
- Geometry-only work is cached while the panel nodes are unchanged:
  `BuildGamma` reuses the influence matrix, its LU and the reference
  solutions (`Isol::AIC_x`), and `CalcUeM` reuses the airfoil-panel part of
  the source influence (`Isol::Bp_foil`), as XFOIL caches the airfoil part
  of its DIJ matrix. `RebuildUeM`, called when the stagnation point moves to
  another panel, now only reapplies the surface signs. A 0.5-degree polar
  step at 200 panels takes 15-18 ms instead of 22-25 ms.
- The coupled viscous Newton system and the viscous adjoint are factored
  through the mass-defect Schur complement (`CoupledFactorization`, new
  header `core/coupled_linear.h`) instead of a sparse LU of the assembled
  `R_V`, whose dense edge-velocity block filled in. Substituting the mass
  defect for displacement thickness leaves an N x N dense system after a
  sparse boundary-layer solve. Steps and adjoints are unchanged to round-off;
  the linear solve per Newton iteration drops from 15 ms to 1.9 ms at 200
  panels. The boundary-layer block is solved by block forward substitution
  over its station groups (`BlockTriangularSolver`). A sparse LU of the assembled system remains as a fallback, and the
  trust-region step still assembles `R_V`.
- `Glob::lu_factorization` is replaced by `Glob::coupled_factorization`, and
  `StoreAdjointFactorization` no longer fills `Glob::R_V`; call
  `coupled_factorization->Assemble()` for the matrix. A viscous `Solve()` no
  longer factors the adjoint system; `ComputeSensitivities` does, on demand.
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
- The CI coverage job runs the full default test suite in parallel, instead
  of the spline, geometry, mfoil reference, and Philote tests only. A
  15-minute job limit and a 60-second per-test limit keep a slow test from
  stalling the pipeline. Eigen, the standard library, GoogleTest, gRPC, and
  protobuf are not instrumented (`.github/coverage-profile-list.txt`); their
  inlined counters made instrumented tests about 17x slower than Release and
  pushed `CalcForce_Viscous` past the limit. The reported coverage of
  Scirocco's sources is unchanged.
- Nonfinite operating conditions, invalid panel counts, malformed coordinate
  options, and nonfinite solver results now raise errors instead of being
  reported as successful analyses.
- Prerelease target versions stay fixed through candidate progression and
  stable promotion; release labels and tags are validated before an atomic
  publish.
- Multipoint Kulfan optimization example switched from OpenMDAO's
  `ScipyOptimizeDriver`/SLSQP to IPOPT via `modopt.OpenMDAOProblem`.
  gRPC server ports shifted to 50060–50063 to avoid collisions. With
  Philote-Cpp 0.6.0, the example initializes each angle of attack at its lift
  target and rejects unsuccessful solves. Verified: `Solve_Succeeded`
  at iter 222, Σw·Cd = 4.601e-3.
- Philote tutorial, Docker, and CI builds now pin Philote-Cpp 0.6.0
  and `philote-mdo==0.8.0`; CI executes the complete
  Python tutorial sequence against a Release build.
- Philote-Cpp pin moved from a pre-release commit (`cd19666`) to the
  v0.6.0 release (`49e1b42`) in CI, release, sanitizer, and Docker builds
  and the install guide. The C++ sources are unchanged, so
  `philote-mdo==0.8.0` still applies.
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
- The test suite is split into a default suite and opt-in diagnostics. The
  default suite includes `scirocco_viscous_chain_tests`, which checks the
  geometric forward chains (`DFoilSDFoilXCumulative`, `DSstagDFoilX`,
  `DWakeSDFoilX`, `DXiDFoilX`, `DWgapDFoilX`, `DWakeXDFoilX`) against
  complex step at `1e-10` rel / `2e-11` abs. The chain-piece isolators build
  only with `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON` and the 14 `probe_*`
  binaries only with `SCIROCCO_BUILD_PROBES=ON`. CI runs the default suite
  with GCC 14 and Clang 18 in Release, plus a job for the opt-in flags.
- `TOL_VISC_SHAPE` tightened from `1e-4` to `1e-8` in
  `sensitivity_reference_test.cpp::ViscousShapeSensitivity`. End-to-end
  adjoint dcl/dcd/dcm vs mfoil complex-step now agrees to ~3e-11
  relative.

### Fixed
- NACA 0012 at alpha 0 (stagnation point on a node) no longer runs all 200
  Newton iterations and a continuation restart: the residual stalls near
  1e-10, above `rtol`, and the round-off floor test now ends the solve.
- Pseudo-transient continuation added `1/tau` at matrix position `(i, i)`,
  but boundary-layer rows are ordered `3 * node + equation` and columns
  `4 * node + variable`, so the term coupled equations to other nodes'
  variables. Each equation now gets the term on its own variable.
- Viscous analysis rejects a closed trailing edge instead of returning wrong
  results (#18). The wake model is singular when the trailing-edge gap
  closes: on NACA 0012 a gap of 1e-9 chord or less never converged, and an
  exactly closed edge could report convergence at a wrong lift (Cl 0.03
  instead of 0.21). `BuildWake` now calls `CheckViscousTrailingEdge`, which
  throws `std::invalid_argument` when the trailing-edge thickness is below
  `kMinViscousTEGap` (1e-6 chord). Inviscid solves on closed trailing edges
  are unchanged. The `KulfanDiscipline` defaults (`te_upper = te_lower = 0`)
  give a closed edge; viscous workflows must set a gap, as all examples do.
- Philote RPC services no longer hold ownership cycles, so server shutdown
  and startup failures release the discipline and its solver state.
- Server startup failures no longer dereference null.
- Empty, malformed, nonfinite, or unordered camber data is rejected before
  the geometry is changed.
- Switching to inviscid mode clears the previous boundary-layer and wake
  results.
- Reference-test scalar readers no longer overwrite the stack: HDF5 values
  are converted to the destination C++ type. Debug test fixtures initialize
  wake coordinates and wake-inclusive velocities consistently.
- Test-suite and coverage documentation now matches the build, and the main
  branch protection requirements and release-token administrator exception
  are recorded.
- Reynolds-number derivatives were wrong by up to 2e-3 (relative) when
  either surface stays laminar to the trailing edge, e.g. the lower surface
  of NACA 2412 at α = 4°, Re = 10⁶. The wake-start shear stress then comes
  from the transition closure `cttr`, which depends on Re, but `WakeSys` had
  no analytic Re partial. It now provides one (`lin_Re_times_Re`, via the new
  `detail::GetCttr_dRet`). `dcl/dcd/dcm_dRe` now match mfoil complex step
  to about 1e-12 with the default and envelope onset models.
- `Spline2DCreate` runs all ten arc-length passes again. The early exit
  added with the tridiagonal spline solve stopped after about three, while
  `SplineArcLengthDeriv` differentiates ten, so shape gradients on re-paneled
  geometry missed the high-precision reference by ~6e-13
  (`RemainingGradientTest/nonunit_chord`). Paneling stays O(n): 1.2 ms at
  99–199 panels.
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

[0.1.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.1.0
[0.2.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.2.0
[0.3.0]: https://github.com/chrislupp/Scirocco/releases/tag/v0.3.0
[Unreleased]: https://github.com/chrislupp/Scirocco/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/chrislupp/Scirocco/releases/tag/v0.3.1
