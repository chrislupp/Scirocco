---
sidebar_position: 7
---

# Solver parameters

A field-by-field reference for the `Param` struct
([`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h)).
Access it through `foil.GetState().param`.

```cpp
foil.GetState().param.verb = 1;
foil.GetState().param.ncrit = 11.0;
foil.GetState().param.use_ptc = false;
```

Start with the default solver parameters. Adjust them when diagnosing
convergence problems. These native state fields are not options on the current
Philote discipline.

## Newton solver

| Field | Default | Notes |
| --- | --- | --- |
| `verb` | `0` | `0` silent, `1` summary per iter, `2` very verbose. Set to `1` when convergence trouble starts. |
| `rtol` | `1e-11` | Newton residual tolerance. Below `1e-12` you start fighting roundoff. |
| `niglob` | `200` | Maximum iterations for a coupled solve |
| `scale_coupled_system` | `true` | Equilibrate equations and variables |
| `use_continuation` | `true` | Retry failed cold starts along an angle/lift path |
| `continuation_max_attempts` | `24` | Maximum continuation attempts |
| `continuation_max_iterations` | `800` | Total extra coupled iterations |
| `use_trust_region` | `false` | Experimental dogleg step control |

## Pseudo-transient continuation

See [PTC concepts](./ptc) for the algorithm.

| Field | Default | Notes |
| --- | --- | --- |
| `use_ptc` | `true` | Master switch |
| `tau_init` | `1e6` | Start nearly Newton |
| `tau_min` | `1e-2` | Maximum damping floor |
| `tau_max` | `1e8` | Effectively "PTC off" ceiling |
| `tau_grow_clamp` | `4.0` | Per-iter growth cap |
| `tau_shrink_clamp` | `0.1` | Per-iter shrink cap |

## Adjoint

Adjoints use the unregularized converged Jacobian. Factorization failures and
inaccurate transpose solves are reported as errors. The former
`adjoint_ridge` option has been removed because a regularized Jacobian does
not differentiate the Python mfoil equations.

## Boundary-layer closures

These shape the Drela/mfoil closure relations. Don't touch them unless you
are intentionally departing from the standard closure.

| Field | Default | Role |
| --- | --- | --- |
| `ncrit` | `9.0` | $e^N$ transition threshold (see [Transition](./transition)) |
| `Cuq` | `1.0` | Scales the $u_q$ term in the shear lag equation |
| `Dlr` | `0.9` | Wall/wake dissipation length ratio |
| `SlagK` | `5.6` | Shear lag constant |
| `CtauC` | `1.8` | Transition shear-state multiplier |
| `CtauE` | `3.3` | Transition shear-state exponent |
| `GA`, `GB`, `GC` | 6.7, 0.75, 18.0 | G–$\beta$ locus constants |

## Thermodynamics

These fields are initialized during solving from the operating conditions.
`SetOper(..., Ma)` stores the Mach input; it does not refresh thermodynamics
immediately. Treat derived fields as solver outputs.

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

`rho0` is derived and overwritten during initialization. The wrapper uses
normalized freestream density and speed; these fields are not a general fluid
input interface. Advanced changes to `gam` or `Tsrat` require separate validation.

## Station flags

These are bookkeeping fields that mark the current boundary-layer station
type while the solver iterates. Read-only from a user's perspective.

| Field | Meaning |
| --- | --- |
| `simi` | True at a similarity (stagnation-region) station |
| `turb` | True at a turbulent station |
| `wake` | True at a wake station |
