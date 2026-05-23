# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Changed
- Split the test suite into a lean default (`cppfoil_tests`, `cppfoil_forward_ref_tests`, `cppfoil_ping_tests`, `cppfoil_sensitivity_ref_tests`) and an opt-in sensitivity-debug suite. Chain-piece tests (`sensitivity_ping`, `viscous_chain`, `viscous_residual_deriv`, `viscous_output_ping`) and the 14 `probe_*` diagnostic binaries now build only when `CPPFOIL_BUILD_SENSITIVITY_DEBUG_TESTS=ON` or `CPPFOIL_BUILD_PROBES=ON`.
- Extended CI: GCC + Clang Release matrix runs the default suite, a separate job exercises the opt-in flags, a Clang coverage job uploads lcov to Codecov, and a Valgrind memcheck job fails on any definite/indirect leak.

### Added
- `CPPFOIL_COVERAGE` CMake option that instruments the library and tests for `llvm-cov` source-based coverage. Default suite measured at 93.4% line / 94.0% function coverage.
- README documentation for all test-suite build flags and the local coverage workflow.

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
[Unreleased]: https://github.com/chrislupp/ProjectTrysil/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/chrislupp/ProjectTrysil/releases/tag/v0.1.0
