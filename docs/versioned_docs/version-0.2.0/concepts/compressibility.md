---
title: Compressibility
description: The Kármán–Tsien correction and Sutherland viscosity.
---

# Compressibility

Set a Mach number in `SetOper` and Scirocco corrects the incompressible
pressure for compressibility. The boundary layer also uses
temperature-dependent viscosity. There's nothing else to configure.

```cpp
foil.SetOper(2.0, -1.0, 3.0e6, 0.5);   // Mach 0.5
```

## Kármán–Tsien correction

The panel method gives the incompressible pressure $C_{p,0}$. The
Kármán–Tsien rule maps it to a compressible value:

$$
C_p = \frac{C_{p,0}}{\beta + \dfrac{M_\infty^2}{1 + \beta}\,\dfrac{C_{p,0}}{2}},
\qquad \beta = \sqrt{1 - M_\infty^2}
$$

For small $C_{p,0}$ this reduces to the Prandtl–Glauert factor $1/\beta$.
The extra term makes suction peaks grow faster than Prandtl–Glauert
predicts, which matches experiment better at higher Mach numbers.

The sonic pressure coefficient is stored in `GetState().param.cps`. If the
surface $C_p$ drops below it, the flow is locally supersonic and the result
is outside what this model can represent.

## Sutherland viscosity

Viscosity in the boundary layer follows Sutherland's law,

$$
\mu(T) = \mu_\text{ref}\left(\frac{T}{T_\text{ref}}\right)^{3/2}\frac{T_\text{ref} + T_s}{T + T_s},
$$

with $T_s / T_\text{ref} = 0.35$ for air (`param.Tsrat`). Change it only if
you have a reason to model another gas.

## Limits

The correction doesn't capture shocks. It's appropriate for subsonic flow
where the surface stays below sonic speed. Mach numbers up to just under 1
are accepted as input, but treat results with supersonic pockets as
unreliable. For transonic flow, use a solver that resolves shocks.

Source: [`bl_properties.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/properties/bl_properties.h),
[`param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h).
