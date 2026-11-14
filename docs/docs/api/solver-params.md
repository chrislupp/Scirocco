---
sidebar_position: 7
---

# Solver parameters

Every field of the `Param` struct
([`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h)).
Set them through `foil.GetState().param` before solving.

```cpp
foil.GetState().param.verb = 1;
foil.GetState().param.ncrit = 11.0;
foil.GetState().param.use_ptc = false;
```

The defaults work for most cases. These fields aren't available on the
Philote server.

## Newton solver

| Field | Default | Notes |
| --- | --- | --- |
| `verb` | `0` | `0` silent, `1` summary per iter, `2` very verbose. Set to `1` when convergence trouble starts. |
| `rtol` | `1e-11` | Newton residual tolerance. Below `1e-12`, round-off gets in the way. |
| `rtol_floor` | `1e-8` | Also converged once the residual is below this and a step no longer halves it: the round-off floor, which on some cases sits above `rtol` |
| `niglob` | `200` | Maximum iterations for a coupled solve |
| `bl_march` | `true` | Before each Newton step, re-solve the boundary layer by a [march](../concepts/viscous-coupling#the-solve) at the current `ue`, as XFOIL does. `false` gives mfoil's iteration |
| `step_limiter` | `StepLimiter::Absolute` | Under-relaxation rule: `Absolute` (mfoil) caps the amplification change at 2, the shear-stress change at 0.05 and the `ue` change at 0.2; `Relative` (XFOIL) scales the step so no change exceeds +150% or −50% of its value |
| `warm_start` | `true` | Start a viscous `Solve()` from the last converged boundary layer of the same object; cold start if that fails |
| `scale_coupled_system` | `true` | Equilibrate the assembled system of the trust-region step |
| `use_continuation` | `true` | Retry failed cold starts along an angle/lift path |
| `continuation_max_attempts` | `24` | Maximum continuation attempts |
| `continuation_max_iterations` | `800` | Total extra coupled iterations |
| `use_trust_region` | `false` | Experimental dogleg step control |

## Pseudo-transient continuation

See [PTC concepts](../concepts/ptc) for the algorithm.

| Field | Default | Notes |
| --- | --- | --- |
| `use_ptc` | `true` | Master switch |
| `tau_init` | `1e6` | Start nearly Newton |
| `tau_min` | `1e-2` | Maximum damping floor |
| `tau_max` | `1e8` | Effectively "PTC off" ceiling |
| `tau_grow_clamp` | `4.0` | Per-iter growth cap |
| `tau_shrink_clamp` | `0.1` | Per-iter shrink cap |

## Adjoint

There's nothing to tune. The adjoint uses the plain converged Jacobian, and
reports a failed or inaccurate solve as an error. (An older `adjoint_ridge`
option was removed: regularizing the Jacobian gives the wrong derivatives.)

## Boundary-layer closures

Constants of the Drela/mfoil closure relations. Leave them alone unless you
mean to change the model.

| Field | Default | Role |
| --- | --- | --- |
| `ncrit` | `9.0` | $e^N$ transition threshold (see [Transition](../concepts/transition)) |
| `transition` | `TransitionModel::CutCell` | `TransitionModel::Captured` blends laminar and turbulent closures smoothly, for smooth gradients; experimental, prescribed angle of attack only (see [Captured transition](../concepts/transition#captured-transition-experimental)) |
| `capt_L`, `capt_ws`, `capt_wc`, `capt_b` | 60, 1, 1, 0 | Captured model: transition length in momentum thicknesses, onset and blend widths in amplification factor, blend offset |
| `onset` | `OnsetModel::Cubic` | Amplification onset ramp; `OnsetModel::Envelope` is fitted to stability theory (see [Onset of amplification](../concepts/transition#onset-of-amplification)) |
| `Cuq` | `1.0` | Scales the $u_q$ term in the shear lag equation |
| `Dlr` | `0.9` | Wall/wake dissipation length ratio |
| `SlagK` | `5.6` | Shear lag constant |
| `CtauC` | `1.8` | Transition shear-state multiplier |
| `CtauE` | `3.3` | Transition shear-state exponent |
| `GA`, `GB`, `GC` | 6.7, 0.75, 18.0 | G–$\beta$ locus constants |

## Thermodynamics

Most of these are computed from the operating point during the solve, so
treat them as outputs.

| Field | Default | Notes |
| --- | --- | --- |
| `Minf` | `0.0` | Freestream Mach |
| `Vinf` | `0.0` | Freestream speed (the wrapper uses normalized speed 1) |
| `muinf` | `0.0` | Freestream dynamic viscosity (derived) |
| `mu0` | `0.0` | Stagnation dynamic viscosity (derived) |
| `rho0` | `1.0` | Stagnation density |
| `H0` | `0.0` | Stagnation enthalpy (derived) |
| `Tsrat` | `0.35` | Sutherland $T_s / T_\text{ref}$ |
| `gam` | `1.4` | Ratio of specific heats — air |
| `KTb`, `KTl`, `cps` | derived | Kármán–Tsien constants and sonic $C_p$, initialized during solving |

Speed and density are normalized, so these fields aren't a way to describe a
different fluid. Changing `gam` or `Tsrat` departs from the validated
model.

## Station flags

Internal flags for the station being evaluated. Don't set them.

| Field | Meaning |
| --- | --- |
| `simi` | True at a similarity (stagnation-region) station |
| `turb` | True at a turbulent station |
| `wake` | True at a wake station |
