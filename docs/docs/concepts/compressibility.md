---
sidebar_position: 4
---

# Compressibility

Scirocco extends its incompressible panel method to subsonic flows with the
**Kármán–Tsien** pressure correction, and accounts for temperature-dependent
viscosity in the boundary layer with the **Sutherland** model. Both are
parameterized in
[`data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h)
and applied by
[`properties/bl_properties.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/properties/bl_properties.h).

## Kármán–Tsien correction

The incompressible result $C_{p,0}$ from the panel solve is mapped to a
compressible $C_p$ via

$$
C_p = \frac{C_{p,0}}{\sqrt{1 - M_\infty^2} + \tfrac{1}{2} M_\infty^2 \,
              \big( 1 + \sqrt{1 - M_\infty^2} \big)^{-1}\, C_{p,0} }
$$

This is equivalent to a Prandtl–Glauert correction with a quadratic
nonlinear term — better than pure PG near the suction-peak, fast to apply
once per integration. The implementation precomputes two scalars

$$
\beta = \sqrt{1 - M_\infty^2}, \qquad \lambda = \frac{M_\infty^2}{(1 + \beta)^2}
$$

stored as `Param::KTb` and `Param::KTl` and updated whenever `Ma` changes.

### Range of validity

Kármán–Tsien is solid through approximately $M_\infty \le 0.7$ on
attached flow. Above that, shock formation breaks the underlying linearized
assumption and you should use a code that resolves compressibility properly
(Euler, RANS). Scirocco will run at higher Mach but the results stop being
physically defensible.

A **sonic $C_p$** is precomputed and stored as `Param::cps`. If the
corrected $C_p$ at any node goes below this, you're crossing the sonic line
locally; the solver does not flag this but it's a useful diagnostic.

## Sutherland viscosity

The viscous boundary-layer system needs the kinematic viscosity at every
station. Scirocco evaluates Sutherland's law:

$$
\mu(T) = \mu_\text{ref} \left(\frac{T}{T_\text{ref}}\right)^{3/2}
                       \frac{T_\text{ref} + T_s}{T + T_s}
$$

with $T_s / T_\text{ref} = 0.35$ as the default (`Param::Tsrat`). The
freestream value `Param::muinf` is set at thermodynamic initialization
(`InitThermo` in `core/inviscid.h`) along with the freestream speed
`Vinf`, stagnation density `rho0`, and stagnation viscosity `mu0`.

## What the user touches

You never set the Kármán–Tsien constants directly — they're derived from
`Ma`:

```cpp
foil.SetOper(alpha, -1.0, Re, 0.3);   // M_infty = 0.3, KT correction active
```

If you want to disable compressibility for a comparison run, just pass
`Ma = 0`. If you want to tune the Sutherland temperature ratio for a
specific working fluid, override it manually:

```cpp
foil.GetState().param.Tsrat = 0.4;
foil.SetOper(alpha, -1.0, Re, Ma);   // SetOper recomputes thermo with the new ratio
```
