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
results. [Transition validation](../validation/transition) compares both
ramps with the stability results. The data, solver, and fit are in
[`scirocco/validation/onset_envelope`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/onset_envelope).

## Captured transition (experimental)

By default, transition sits inside one panel: stations upstream use laminar
equations and stations downstream turbulent ones. When transition moves
across a node, the station switches equations. Outputs stay continuous, but
their derivatives jump. On NACA 2412 at Re = 10⁶, dCd/dα saw-tooths by
±25% as α changes.

The captured model avoids the switch. It carries the amplification factor
and the shear stress at every node and blends the laminar and turbulent
closures over a transition length:

```cpp
foil.GetState().param.transition = scirocco::TransitionModel::Captured;
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.Solve();
```

- **Smoothness:** dCd/dα becomes smooth. The remaining ripple is 0.2–8% at
  160 panels and shrinks as the panels are refined. The cut-cell ripple
  stays at 20–30%.
- **Drag:** checked against measured Eppler 387 drag (UIUC, Re = 1–4.6×10⁵,
  laminar separation bubbles throughout). The RMS error is 5.3%, against 16%
  for the default model, which under-predicts bubble drag by up to 28%. At
  Re ≥ 10⁶ the two models agree within 4% on NACA airfoils.
- **Transition length:** set by `capt_L`, in momentum thicknesses. The
  default of 60 is calibrated on the Eppler 387 data. The zone should span
  several panels; at high Reynolds numbers, use more panels.
- **Gradients:** `ComputeSensitivities()` works as usual: shape, angle of
  attack, Reynolds and Mach number derivatives all match complex step
  through the mfoil prototype: shape to 1e-8, the scalars to 1e-10.
- **Limits:** prescribed angle of attack only; a lift target throws.

[Transition validation](../validation/transition) has the full comparison
with the Eppler 387 measurements, the gradient sweeps and the derivative
checks. The formulation and the dead ends are in
[`scirocco/validation/captured_transition`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/captured_transition).

## Forced transition

A trip strip fixes where the boundary layer turns turbulent. Set its
position on each surface as a fraction of the chord:

```cpp
scirocco::Scirocco foil("2412", 160);
foil.SetOper(4.0, -1.0, 1e6, 0.0);  // alpha = 4 deg, Re = 1e6
foil.GetState().param.xtrip_upper = 0.2;
foil.GetState().param.xtrip_lower = 0.6;
foil.Solve();

const auto& vsol = foil.GetState().vsol;
std::printf("cl = %.4f  cd = %.5f  upper: x/c = %.3f %s  lower: x/c = %.3f %s\n",
            foil.GetCl(), foil.GetCd(),
            vsol.Xt(1, 1), vsol.tripped[1] ? "(trip)" : "(free)",
            vsol.Xt(0, 1), vsol.tripped[0] ? "(trip)" : "(free)");
```

```output
cl = 0.6671  cd = 0.00927  upper: x/c = 0.200 (trip)  lower: x/c = 0.600 (trip)
```

Without the trips the same case gives cl = 0.7202 and cd = 0.00693, with
transition at x/c = 0.397 on the upper surface and a lower surface that
stays laminar to the trailing edge.

The rule is XFOIL's: transition happens at the trip or where the
amplification reaches `ncrit`, whichever comes first. A trip behind free
transition changes nothing, and the default of 1 leaves a surface free.
After a solve, `vsol.Xt(si, 1)` holds the transition x on the lower
(`si = 0`) and upper (`si = 1`) surface, and `vsol.tripped[si]` tells
whether the trip set it.

Use it to model a roughness strip or a turbulator, or to compare with
tripped wind-tunnel data: see
[the tripped NACA 0012](../validation/experiments#tripped-drag).
Derivatives with respect to shape, angle of attack, Reynolds number and
Mach number account for the trip and match complex-step references. The
trip stays at its x/c, so it does not move with the shape.

Trips work in both transition models. With
[captured transition](#captured-transition-experimental) the trip acts as
onset does: the amplification factor jumps to `ncrit` at the trip, and the
closures then blend from laminar to turbulent over `capt_L` momentum
thicknesses behind it. The case above gives cd = 0.00925 with the captured
model. A trip behind the end of the transition zone changes nothing; one
inside the zone shortens it slightly. `vsol.Xt` and `vsol.tripped` are
those of the cut-cell solution the captured solve starts from.

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

The amplification model is the standard one from XFOIL and mfoil, with the
same strengths and weaknesses: it is good for smooth airfoils at moderate
Reynolds numbers, and does not model bypass transition from high
free-stream turbulence.

[Forced transition](#forced-transition) has these limits:

- **Trips close behind the stagnation point may not converge.** There the
  boundary layer is too thin for the turbulent closure ($Re_\theta$ below
  about 100), and `Solve()` returns `false`. This rules out a trip at
  x/c = 0 as a way to get a fully turbulent solution, and can affect the
  lower surface at high angles of attack, where the stagnation point sits
  just ahead of the trip. Move the trip a few percent of chord downstream.
  XFOIL converges in these cases.
- **A trip the stagnation point has passed is ignored.** At a high angle of
  attack the stagnation point can move behind a lower-surface trip near the
  leading edge. That surface then transitions freely, as in XFOIL.
- **The trip is placed by x, not along a chord line.** It is at
  x = `xtrip` × chord, interpolated linearly between panel nodes, which
  assumes the leading edge at x = 0.
- **Not on the Philote server**, which does not expose solver parameters.

Source: [`transition.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/transition.cpp),
[`bl_properties.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/src/bl_properties.cpp).
