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

## Onset of amplification

Disturbances only grow above a critical $Re_\theta$, which depends on the
shape factor. The default model, from XFOIL and mfoil, switches growth on
with a short cubic ramp around that point. You can opt in to a ramp fitted to
linear stability theory instead:

```cpp
foil.GetState().param.onset = scirocco::OnsetModel::Envelope;
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.Solve();
```

<OnsetPlot />

The envelope ramp comes from spatial Orr–Sommerfeld analysis of the
Falkner–Skan profiles, for $H$ from 2.36 to 3.48. Growth starts at the
critical Reynolds number and then approaches the asymptotic rate slowly. Far
downstream the asymptote is the same as the default ramp's, so the
amplification factor there is unchanged. On NACA 2412 at α = 4°,
Re = 10⁶, $C_d$ changes by 0.06% and transition moves 0.0007c.

The ramp also matters for gradients. The default ramp is sharp, so the drag
gradient with respect to the shape near the onset is large and oscillates.
For a smooth bump 0.01c wide, the envelope ramp makes that gradient two to
four times smaller and removes the oscillation. Use it for shape
optimization. Leave the default when you need to match XFOIL or mfoil
results. The data, solver, and fit are in
[`scirocco/validation/onset_envelope`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/onset_envelope).

## How transition moves during a solve

Transition is part of the Newton iteration: as the solution changes, the
transition point can move from one panel to the next. When it does, the
equations at that station switch from laminar to turbulent. This is the
main source of convergence trouble and of noisy finite-difference
gradients. Check gradients at points where transition sits still.

At transition, the turbulent shear stress starts from a value set by the
local boundary-layer state through the closure constants `CtauC` and `CtauE`.

## Limits

The model turns amplification on over a narrow band of $Re_\theta$, a
fixed numerical ramp inherited from XFOIL and mfoil. Shape gradients are
exact for the discrete model, but they can spike from node to node where the
amplification starts to grow, and refining the paneling does not remove the
spikes. See
[where the shape gradient is rough](../guides/gradients#where-the-shape-gradient-is-rough).

There is no way to force transition (a trip strip) through the public API
or the Philote server. The amplification model is the standard one from
XFOIL and mfoil, with the same strengths and weaknesses: it is good for
smooth airfoils at moderate Reynolds numbers, and does not model bypass
transition from high free-stream turbulence.

Source: [`transition.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/transition.cpp),
[`bl_properties.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/bl_properties.cpp).
