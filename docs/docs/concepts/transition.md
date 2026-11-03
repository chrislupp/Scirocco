---
sidebar_position: 3
---

# Transition

Scirocco uses an amplification model to predict laminar-to-turbulent transition
in its integral boundary layer. The default threshold is `ncrit = 9`.

## Amplification threshold

The laminar state tracks an amplification quantity along each airfoil surface.
`MarchAmplification` integrates its growth using the local boundary-layer
state. `UpdateTransition` updates the laminar/turbulent station classification
as the coupled solution changes.

Lower `ncrit` values generally trigger earlier transition; higher values allow
more amplification before transition. Select the threshold to represent the
flow environment and compare against appropriate validation data. It is not a
universal material property or a direct turbulence-intensity input.

In C++, set it before solving:

```cpp
foil.GetState().param.ncrit = 9.0;
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
bool converged = foil.Solve();
```

The current Philote discipline does not expose `ncrit` as an option, and the
public sensitivity interface does not provide derivatives with respect to it.

## Transition and convergence

A change in the transition station changes the equations used at that station.
This can complicate convergence and finite-difference gradient checks near a
transition change. Start gradient verification from a converged, repeatable
operating point; see [Troubleshooting](../getting-started/troubleshooting).

The initial turbulent shear state is computed by `GetCttr` from the local
boundary-layer state and the `CtauC` and `CtauE` closure constants. It is not a
log-linear function of the accumulated amplification.

## Forced transition

There is no supported forced-transition input on the public `Solve()` or
Philote interface. `StoreTransition` records internal transition information;
calling it does not establish a persistent trip constraint, because the solver
updates transition during iteration.

Implementation references:

- [`transition.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/transition.cpp)
- [`bl_properties.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/bl_properties.cpp)
