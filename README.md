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
- Philote gRPC servers, so any MDO framework with a Philote client can drive
  the analysis

**Documentation: <https://chrislupp.github.io/Scirocco/>**

## Quick start

### C++

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

### Python and MDO frameworks, with Docker

From Python or an MDO framework, Scirocco runs as two Philote gRPC servers
that any Philote client connects to; the examples use OpenMDAO. There is no
in-process `import scirocco` module. Both servers come as Docker images, so
there's nothing to compile:

| Image | Port | Server |
| --- | --- | --- |
| [`chrislupp/scirocco`](https://hub.docker.com/r/chrislupp/scirocco) | 50051 | Analysis: shape and operating point to Cl, Cd, Cm, with gradients |
| [`chrislupp/scirocco-kulfan`](https://hub.docker.com/r/chrislupp/scirocco-kulfan) | 50052 | Geometry: CST coefficients to panel coordinates |

```bash
docker run -d --rm --name scirocco -p 50051:50051 chrislupp/scirocco:latest
docker run -d --rm --name scirocco-kulfan -p 50052:50052 chrislupp/scirocco-kulfan:latest

git clone https://github.com/chrislupp/Scirocco.git
cd Scirocco
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r examples/tutorials/requirements.txt
python examples/tutorials/quickstart.py
```

`latest` is the newest release. Pin a version tag, such as `0.2.0`, for
results that don't change under you; `next` tracks the `develop` branch. See
[Docker](https://chrislupp.github.io/Scirocco/docs/python/docker) for the
tags, and stop the servers with `docker stop scirocco scirocco-kulfan`.

## Documentation

Everything else is on the
[documentation site](https://chrislupp.github.io/Scirocco/). To build the
tests or contribute a change, see [CONTRIBUTING.md](CONTRIBUTING.md).
Changes between releases are in [CHANGELOG.md](CHANGELOG.md).

## Citation

Scirocco implements the methods of the papers below. If you use it in
published work, please cite them:

- K. J. Fidkowski, "A Coupled Inviscid–Viscous Airfoil Analysis Solver,
  Revisited," *AIAA Journal*, Vol. 60, No. 5, 2022, pp. 2961–2971.
  [doi:10.2514/1.J061341](https://doi.org/10.2514/1.J061341). This is the
  mfoil solver that Scirocco ports: the coupled Newton method, the
  stagnation-point treatment, and the closure relations.
- M. Drela, "XFOIL: An Analysis and Design System for Low Reynolds Number
  Airfoils," in T. J. Mueller (ed.), *Low Reynolds Number Aerodynamics*,
  Lecture Notes in Engineering, Vol. 54, Springer, 1989, pp. 1–12.
  [doi:10.1007/978-3-642-84010-4_1](https://doi.org/10.1007/978-3-642-84010-4_1).
  This is the source of the integral boundary-layer formulation, the e^N
  transition model, and the boundary-layer march between Newton steps.

```bibtex
@article{Fidkowski2022mfoil,
  author  = {Fidkowski, Krzysztof J.},
  title   = {A Coupled Inviscid--Viscous Airfoil Analysis Solver, Revisited},
  journal = {AIAA Journal},
  volume  = {60},
  number  = {5},
  pages   = {2961--2971},
  year    = {2022},
  doi     = {10.2514/1.J061341}
}

@incollection{Drela1989xfoil,
  author    = {Drela, Mark},
  title     = {{XFOIL}: An Analysis and Design System for Low {R}eynolds
               Number Airfoils},
  booktitle = {Low {R}eynolds Number Aerodynamics},
  editor    = {Mueller, Thomas J.},
  series    = {Lecture Notes in Engineering},
  volume    = {54},
  pages     = {1--12},
  publisher = {Springer},
  address   = {Berlin, Heidelberg},
  year      = {1989},
  doi       = {10.1007/978-3-642-84010-4_1}
}
```

## License

Apache License 2.0. See [LICENSE](LICENSE).

## Acknowledgments

Scirocco is a C++ port of [mfoil](https://websites.umich.edu/~kfid/codes.html)
by Krzysztof J. Fidkowski, a Python code for subsonic airfoil analysis
released under the MIT License. See [NOTICE](NOTICE).
