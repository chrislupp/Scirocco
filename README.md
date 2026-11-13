# Scirocco

[![Build](https://github.com/chrislupp/Scirocco/actions/workflows/build.yml/badge.svg)](https://github.com/chrislupp/Scirocco/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/chrislupp/Scirocco/branch/main/graph/badge.svg)](https://codecov.io/gh/chrislupp/Scirocco/branch/main)

A C++17 library for inviscid and viscous analysis of 2D airfoils, with
adjoint gradients for design optimization.

Give Scirocco a shape and an operating point, and it returns lift, drag,
moment, and the boundary layer. One more call returns the derivatives of
those results with respect to the shape, angle of attack, Reynolds number,
and Mach number.

- Linear-vortex panel method coupled by Newton's method to an integral
  boundary layer and wake
- e^N transition prediction and Kármán–Tsien compressibility correction
- NACA 4- and 5-digit sections, custom coordinates, flaps, and CST (Kulfan)
  geometry
- Discrete adjoint for Cl, Cd, and Cm, tested against complex-step
  references
- Philote gRPC servers for OpenMDAO and other Python MDO frameworks

**Documentation: <https://chrislupp.github.io/Scirocco/>**

## Quick start

You need CMake 3.18+, a C++17 compiler, and Eigen 3.4+.

```bash
git clone https://github.com/chrislupp/Scirocco.git
cd Scirocco
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j --target quickstart
scirocco/build/quickstart
```

That builds and runs this program, a viscous solve of a NACA 2412:

```cpp
#include <scirocco/scirocco.h>

scirocco::Scirocco foil("2412", 160);  // NACA 2412, 160 panels
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);   // alpha = 4 deg, Re = 1e6, Mach 0
foil.Solve();

double cd = foil.GetCd();
foil.ComputeSensitivities();           // adjoint gradients
const auto& dcd_dx = foil.GetDcdDx();  // dCd / d(node coordinates)
```

From Python, Scirocco runs as a Philote server that OpenMDAO connects to.
There is no in-process `import scirocco` module.

## Where to go next

| If you want to | Read |
| --- | --- |
| Install, including the optional HDF5 and Philote dependencies | [Installation](https://chrislupp.github.io/Scirocco/docs/getting-started/installation) |
| Run analyses, change geometry, or compute gradients in C++ | [Guides](https://chrislupp.github.io/Scirocco/docs/guides/analysis) |
| Use Scirocco from Python and OpenMDAO | [Python overview](https://chrislupp.github.io/Scirocco/docs/python) |
| Optimize an airfoil | [Optimization examples](https://chrislupp.github.io/Scirocco/docs/optimization/cst-design) and [`examples/`](examples/) |
| Look up a class, option, or CMake flag | [API reference](https://chrislupp.github.io/Scirocco/docs/api) |
| Understand the method, transition models, and their validation | [Concepts](https://chrislupp.github.io/Scirocco/docs/concepts/panel-method) |
| Fix a convergence or gradient failure | [Troubleshooting](https://chrislupp.github.io/Scirocco/docs/getting-started/troubleshooting) |
| Build the tests or contribute a change | [CONTRIBUTING.md](CONTRIBUTING.md) |

Changes between releases are in [CHANGELOG.md](CHANGELOG.md).

## License

Apache License 2.0. See [LICENSE](LICENSE).

## Acknowledgments

Scirocco is a C++ port of [mfoil](https://websites.umich.edu/~kfid/codes.html)
by Krzysztof J. Fidkowski, a Python code for subsonic airfoil analysis
released under the MIT License. See [NOTICE](NOTICE).
