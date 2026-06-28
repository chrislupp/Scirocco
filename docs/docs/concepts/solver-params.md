---
draft: true
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

The defaults are tuned to be sane out-of-the-box; the table below explains
when each field is worth touching.

## Newton solver

| Field | Default | Notes |
| --- | --- | --- |
| `verb` | `0` | `0` silent, `1` summary per iter, `2` very verbose. Set to `1` when convergence trouble starts. |
| `rtol` | `1e-11` | Newton residual tolerance. Below `1e-12` you start fighting roundoff. |
| `niglob` | `200` | Max Newton iterations. Raised from 50 to give PTC room. |

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

| Field | Default | Notes |
| --- | --- | --- |
| `adjoint_ridge` | `1e-10` | Diagonal ridge added to $R_V$ if SparseLU fails. Biases the adjoint by $O(\varepsilon \cdot \|R_V\|)$ but rescues near-singular configurations. |

## Boundary-layer closures

These shape the Drela/mfoil closure relations. Don't touch them unless you
are intentionally departing from the standard closure.

| Field | Default | Role |
| --- | --- | --- |
| `ncrit` | `9.0` | $e^N$ transition threshold (see [Transition](./transition)) |
| `Cuq` | `1.0` | Scales the $u_q$ term in the shear lag equation |
| `Dlr` | `0.9` | Wall/wake dissipation length ratio |
| `SlagK` | `5.6` | Shear lag constant |
| `CtauC` | `1.8` | $\log_{10}\,C_\tau$ constant after transition |
| `CtauE` | `3.3` | $\tilde n$ coefficient in the $C_\tau$ seed |
| `GA`, `GB`, `GC` | 6.7, 0.75, 18.0 | G–$\beta$ locus constants |

## Thermodynamics

Most of these are filled by `InitThermo` (called inside `SetOper`) — you set
`Minf` indirectly through `SetOper(..., Ma)` and the rest follows.

| Field | Default | Notes |
| --- | --- | --- |
| `Minf` | `0.0` | Freestream Mach |
| `Vinf` | `0.0` | Freestream speed (derived from `Re` and `mu0`/`rho0`) |
| `muinf` | `0.0` | Freestream dynamic viscosity (derived) |
| `mu0` | `0.0` | Stagnation dynamic viscosity (derived) |
| `rho0` | `1.0` | Stagnation density |
| `H0` | `0.0` | Stagnation enthalpy (derived) |
| `Tsrat` | `0.35` | Sutherland $T_s / T_\text{ref}$ |
| `gam` | `1.4` | Ratio of specific heats — air |
| `KTb`, `KTl`, `cps` | derived | Kármán–Tsien constants and sonic $C_p$, recomputed each `SetOper` |

If you want non-air gas properties, change `rho0`, `gam`, and possibly
`Tsrat` **before** calling `SetOper`; `SetOper` runs `InitThermo` which
respects what you've set.

## Station flags

These are bookkeeping fields that mark the current boundary-layer station
type while the solver iterates. Read-only from a user's perspective.

| Field | Meaning |
| --- | --- |
| `simi` | True at a similarity (stagnation-region) station |
| `turb` | True at a turbulent station |
| `wake` | True at a wake station |
