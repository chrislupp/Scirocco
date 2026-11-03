---
sidebar_position: 2
---

# Viscous coupling

A positive Reynolds number enables the integral boundary-layer model on both
airfoil surfaces and in the wake. Boundary-layer displacement changes the outer
flow; the outer-flow edge velocity, in turn, changes the boundary layer.
Scirocco solves this interaction as a coupled nonlinear system.

## State and coupling

The global viscous state has four values at each airfoil or wake station:

| State | Meaning |
| --- | --- |
| `theta` | Momentum thickness |
| `delta_star` | Displacement thickness |
| Amplification / shear state | Laminar amplification or turbulent shear-stress variable |
| `ue` | Edge velocity |

The displacement effect enters through the mass-source quantity
$m = u_e\delta^*$. An influence matrix maps changes in this quantity to changes
in edge velocity. Nodal vortex strengths belong to the inviscid solution; they
are not additional rows in the four-row viscous state array.

Initialization builds an inviscid solution, wake geometry, stagnation/surface
indexing, and a boundary-layer starting state. Coupled iterations assemble the
residual and Jacobian, update the state, and track changes in stagnation and
transition. They do not restart the entire boundary-layer march on every step.

## Convergence

The solver uses a sparse coupled linear system with scaling and step control.
[Pseudo-transient continuation](./ptc) adds damping when needed; continuation
retries can help recover failed cold starts. Success requires the residual to
meet the configured tolerance. Inspect the return value from `Solve()` before
using forces or sensitivities.

The converged, unregularized Jacobian is used by the
[adjoint](../theory). A failed primal or adjoint does not provide usable
gradients. See [Troubleshooting](../getting-started/troubleshooting).

## Outputs

| Quantity | Source |
| --- | --- |
| `cl`, `cm` | Surface pressure integration |
| `cd` | Far-wake momentum deficit with Squire–Young extrapolation |
| `cdf` | Integrated skin-friction contribution |
| `cdp` | `cd - cdf` |
| `cdpi` | Separate near-field pressure integration |

Boundary-layer distributions include airfoil and wake nodes; pressure contains
only airfoil nodes. See [Conventions](../api/conventions) for sizes and reference
quantities. Integral boundary-layer modeling has limits in strongly separated
flow; numerical convergence alone does not establish physical accuracy.

Implementation references:

- [`boundary_layer.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/boundary_layer.cpp)
- [`coupled_solver.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/coupled_solver.cpp)
- [`forces.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/forces.cpp)
