---
title: Transition
description: The e^N amplification model and the ncrit parameter.
---

# Transition

Scirocco predicts where the laminar boundary layer turns turbulent with the
$e^N$ method. Along each surface, it tracks the amplification $\tilde n$ of
the most unstable disturbance. Where $\tilde n$ reaches the critical value
$n_\text{crit}$, the flow becomes turbulent.

## Choosing `ncrit`

The default is `ncrit = 9`, the standard value for an average wind tunnel.
Quieter environments, such as clean free flight, call for higher values,
which move transition back. Noisier ones call for lower values, which move
it forward.

```cpp
foil.GetState().param.ncrit = 9.0;   // set before solving
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.Solve();
```

`ncrit` is a modelling choice, not a measured property. Pick it to match
your flow environment, and validate against data where you can. It isn't an
option on the Philote server, and there are no derivatives with respect to it.

## How transition moves during a solve

Transition is part of the Newton iteration: as the solution changes, the
transition point can move from one panel to the next. When it does, the
equations at that station switch from laminar to turbulent. This is the
main source of convergence trouble and of noisy finite-difference
gradients. Check gradients at points where transition sits still.

At transition, the turbulent shear stress starts from a value set by the
local boundary-layer state through the closure constants `CtauC` and `CtauE`.

## Limits

There is no way to force transition (a trip strip) through the public API
or the Philote server. The amplification model is the standard one from
XFOIL and mfoil, with the same strengths and weaknesses: it is good for
smooth airfoils at moderate Reynolds numbers, and does not model bypass
transition from high free-stream turbulence.

Source: [`transition.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/transition.cpp),
[`bl_properties.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/bl_properties.cpp).
