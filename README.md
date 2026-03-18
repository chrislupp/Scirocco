# cppfoil

A modern C++ library for inviscid and viscous aerodynamic analysis of airfoils.

## Overview

cppfoil provides computational tools for analyzing airfoil performance across different operating conditions. It predicts aerodynamic coefficients (lift, drag, moment) and boundary layer distributions using panel methods coupled with integral boundary layer equations.

Key capabilities:
- Inviscid (potential flow) analysis using linear vortex panel methods
- Viscous analysis with coupled boundary layer solver
- Laminar-to-turbulent transition prediction (e^N method)
- Compressibility corrections for subsonic flow
- NACA 4-digit and 5-digit airfoil generation
- Geometry modification (flaps, camber, derotation)

## Requirements

- CMake 3.18+
- C++17 compatible compiler
- Eigen 3.4+ (linear algebra library)

Optional:
- Google Test (for unit tests)
- HDF5 (for reference validation tests)

## Building

```bash
cd cppfoil
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

Build options:
```bash
cmake -DCPPFOIL_BUILD_TESTS=ON -DCPPFOIL_BUILD_EXAMPLES=ON ..
```

## Usage

### Basic Inviscid Analysis

```cpp
#include <cppfoil/cppfoil.h>

cppfoil::Cppfoil foil("0012", 199);  // NACA 0012, 199 panels
foil.SetOper(5.0, -1.0, 0.0, 0.0);  // alpha=5°, no Re (inviscid)
foil.Solve();

std::cout << "Cl = " << foil.GetCl() << std::endl;
std::cout << "Cm = " << foil.GetCm() << std::endl;
```

### Viscous Analysis

```cpp
#include <cppfoil/cppfoil.h>

cppfoil::Cppfoil foil("2412", 129);  // NACA 2412, 129 panels
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
Eigen::MatrixXd coords(n, 2);  // n points, columns are x and z
// ... fill coordinates ...

cppfoil::Cppfoil foil(coords, 199);
```

### Geometry Modification

```cpp
cppfoil::Cppfoil foil("0012", 199);

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

## API Reference

### Cppfoil Class

| Method | Description |
|--------|-------------|
| `Cppfoil(naca, npanel)` | Create airfoil from NACA designation |
| `Cppfoil(coords, npanel)` | Create airfoil from coordinates |
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
cppfoil/
├── include/cppfoil/
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

```bash
cd build
ctest
# or
./cppfoil_tests
```

The test suite validates against reference data from mfoil to ensure numerical accuracy.

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

MIT License

## Acknowledgments

This library is validated against [mfoil](https://github.com/), a Python aerodynamic analysis tool, ensuring numerical accuracy and consistency.
