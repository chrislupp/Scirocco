---
sidebar_position: 1
---

# Scirocco

A modern C++ library for inviscid and viscous aerodynamic analysis of 2D
airfoils, with Philote servers and Python bindings for scripted analysis and
optimization.

Scirocco predicts force coefficients ($C_l$, $C_d$, $C_m$) and boundary-layer
distributions using a linear-vortex panel method coupled with integral
boundary-layer equations and $e^N$ transition prediction. It is a C++17 port and
reformulation of [mfoil](https://websites.umich.edu/~kfid/codes.html), with
analytic adjoints, HDF5 output, and a Philote-Cpp discipline interface for use
inside MDO frameworks.

## Capabilities

| Area | What's in the box |
| --- | --- |
| Inviscid | Linear-vortex panel method, Kutta condition, AIC formulation |
| Viscous | Coupled boundary-layer / inviscid Newton solve, wake displacement |
| Transition | $e^N$ method with configurable $n_\text{crit}$ |
| Compressibility | Kármán–Tsien correction, Sutherland viscosity |
| Geometry | NACA 4- and 5-digit generation, CST/Kulfan parameterization, flap, camber, derotation |
| Sensitivities | Full reverse-mode adjoint w.r.t. coordinates, $\alpha$, $Re$, $M_\infty$ |
| I/O | HDF5 iteration output (optional) |
| MDO | Philote-Cpp gRPC discipline; OpenMDAO examples for NACA and CST design |

## Who this is for

- Aerodynamicists who want a fast, scriptable 2D analysis tool with derivatives
- MDO practitioners who need an airfoil discipline with analytic gradients
- Researchers working on airfoil optimization

## Where to next

- [Install Scirocco](./getting-started/installation) and confirm a clean build
- Run the 5-minute [Python quickstart](./getting-started/quickstart) through
  the Philote-MDO bindings
- Browse the [API Reference](./api) for the public surface


## License & attribution

Apache License 2.0. Scirocco is a C++ port of mfoil (MIT) by Krzysztof J.
Fidkowski; see the project [`NOTICE`](https://github.com/chrislupp/Scirocco/blob/develop/NOTICE)
for details.
