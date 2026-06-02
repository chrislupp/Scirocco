---
sidebar_position: 3
---

# Transition

Scirocco predicts laminar-to-turbulent transition with the $e^N$ method.
This page describes the model, the one knob you typically care about
($n_\text{crit}$), and what the solver does internally with the result.
Implementation in
[`core/transition.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/transition.h).

## The $e^N$ method in one paragraph

Small perturbations of a laminar boundary layer amplify according to local
linear stability theory. The integrated growth rate from the neutral point
along a streamwise station $s$ is the amplification factor

$$
\tilde n(s) = \int_{s_\text{neutral}}^{s} \frac{d\tilde n}{ds}\, ds'
$$

Tollmien–Schlichting waves are observed to break down to turbulence when
$\tilde n$ reaches a problem-dependent threshold $n_\text{crit}$. The
classical mfoil/XFOIL closure gives $d\tilde n / ds$ as a function of
shape factor $H$ and momentum-thickness Reynolds number $Re_\theta$.

That's it: march $\tilde n$ from the stagnation point downstream; wherever
$\tilde n$ first crosses $n_\text{crit}$, declare the transition station
and hand off to the turbulent closure for the rest of the surface.

## The one knob: $n_\text{crit}$

```cpp
foil.GetState().param.ncrit = 9.0;   // default
```

Rule of thumb:

| Environment | $n_\text{crit}$ |
| --- | --- |
| Very clean wind tunnel (sailplane, low Tu) | 13 – 14 |
| Average wind tunnel | 9 (default) |
| Dirty free flight | 7 – 8 |
| Heavy free-stream turbulence (helicopter rotor) | 4 – 5 |

Lower values transition earlier (more turbulent surface, higher friction
drag, less laminar-bucket performance). Higher values transition later.
For sensitivity studies, $n_\text{crit}$ is a meaningful design knob — but
it is not a coordinate, so adjoint w.r.t. $n_\text{crit}$ is not exposed
in the current API.

## What `MarchAmplification` does

For each surface (upper, lower, and each wake side), `MarchAmplification`:

1. Initializes $\tilde n = 0$ at the stagnation point.
2. Steps station-by-station along the laminar segment, computing
   $d\tilde n / ds$ from the local $H$ and $Re_\theta$.
3. Returns the index of the **last laminar station** before $\tilde n$
   exceeds $n_\text{crit}$.

That index then drives `StoreTransition` and `UpdateTransition`, which
mark the transition location for the coupled Newton system.

## Forced transition

To pin transition to a specific location (e.g. simulating a boundary-layer
trip), set the transition index directly via `StoreTransition` before
calling `SolveCoupled`. This bypasses the $e^N$ march on that surface.
This is intentionally a lower-level intervention: the public `Solve()`
path leaves $e^N$ in charge.

## Initial $C_\tau$ after transition

The first turbulent station inherits a finite shear-stress coefficient
seeded from the $C_\tau$ closure:

$$
\log_{10} C_\tau\big|_{0,\text{turb}} = C_{\tau,C} + C_{\tau,E} \,\tilde n_\text{trans}
$$

The constants `CtauC = 1.8` and `CtauE = 3.3` (see
[`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h))
give the standard mfoil seed.
