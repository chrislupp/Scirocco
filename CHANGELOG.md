# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/).

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
