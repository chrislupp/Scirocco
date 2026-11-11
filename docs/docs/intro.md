---
title: Overview
sidebar_label: Overview
description: What Scirocco computes, how it works, and where to start.
---

# Scirocco

Scirocco is a C++17 library for analyzing 2D airfoils. Give it a shape and an
operating point, and it returns lift, drag, moment, and the boundary layer.
One more call returns the derivatives of those results with respect to every
input, which is what a gradient-based optimizer needs.

<CpPlot />

## What it does

- **Inviscid and viscous analysis.** A linear-vortex panel method, coupled by
  Newton's method to an integral boundary layer and wake when `Re > 0`.
- **Transition.** $e^N$ amplification, with `ncrit = 9` by default.
- **Compressibility.** Kármán–Tsien correction for subsonic Mach numbers.
- **Geometry.** NACA 4- and 5-digit sections, your own coordinates, flaps,
  camber changes, and a CST (Kulfan) parameterization.
- **Gradients.** A discrete adjoint for $C_l$, $C_d$, and $C_m$ with respect
  to node coordinates, $\alpha$, $Re$, and $M_\infty$, tested against
  complex-step reference data.
- **Output and integration.** Optional HDF5 history files, and Philote gRPC
  servers that plug into OpenMDAO.

Scirocco is a C++ port of [mfoil](https://websites.umich.edu/~kfid/codes.html)
by Krzysztof J. Fidkowski, extended with analytic adjoints. The Python mfoil
implementation is the reference for every result and derivative in the test
suite.

## Pick your path

<Next items={[
  { to: "/docs/getting-started/quickstart", title: "C++ quickstart", body: "Build the library and solve your first airfoil in a few minutes." },
  { to: "/docs/python", title: "Python and OpenMDAO", body: "Drive Scirocco from Python through its Philote servers." },
  { to: "/docs/concepts/panel-method", title: "How it works", body: "The panel method, boundary layer, transition, and adjoint." },
]} />

Scirocco is licensed under Apache 2.0. mfoil is MIT-licensed; see
[`NOTICE`](https://github.com/chrislupp/Scirocco/blob/develop/NOTICE).
