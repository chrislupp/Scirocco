# Plan: replace DRVDRe / DRVDMa FD with analytic chains

**Status:** not yet implemented. Current `DRVDRe` (`adjoint.cpp:1189`) and
`DRVDMa` (`adjoint.cpp:DRVDMa`) use 4th-order central FD on `mu0` and
`Ma` respectively. This violates the "no FD ever" mandate and caps
Re/Ma sensitivities at ~1e-3 / 1e-4 rel.

**Caveat:** end-to-end Re/Ma test floors will still hit the 1e-5 lambda
bug ([BUG_rv_ue_istag_column.md](BUG_rv_ue_istag_column.md)) until that
is resolved. Implementing analytic Re/Ma alone won't fix the end-to-end
tests; it will fix the *partials* (which currently aren't accurate
enough to land below the lambda noise).

## Re chain (smaller scope)

Re enters R_V *only* through `mu0` (set in `InitThermo`), and `mu0`
enters *only* through `Ret = rho0 * ue * th / mu0` in BL closures.

Chain rule:
```
∂R_V/∂Re = ∂R_V/∂mu0 · ∂mu0/∂Re
∂mu0/∂Re = -mu0/Re             (analytic, from InitThermo formula)
∂R_V/∂mu0 = Σ_station ∂R_station/∂Ret_station · ∂Ret_station/∂mu0
∂Ret/∂mu0 = -Ret/mu0
```

So the only piece to implement is `∂R_station/∂Ret_station` for every
station.

### Closures that use Ret (need `_dRet` companions)

From `bl_properties.cpp`:
- `GetCf`         (lam & turb branches)
- `GetCDi_Lam`
- `GetCDi_TurbWall`
- `GetCDi_LamStress`     (CDi_outer doesn't use Ret directly)
- `GetHs`         (turb branch)
- `GetCteq`
- `GetUq`
- `GetUs`         (turb branch)
- `GetDamp`       (amplification rate, laminar transition)

Each closure's Ret-derivative is a simple inspection of its existing
formula. Examples:

**GetCf laminar** (`cf = num/Ret`):
```cpp
double GetCf_dRet(const Vector4d& U, const Param& param) {
    if (param.wake) return 0.0;
    if (param.turb) {
        // turbulent branch (see below)
    } else {
        auto [Hk, _] = GetHk(U, param);
        auto [Ret, __] = GetRet(U, param);
        double num; /* ... computed same as in GetCf ... */
        return -num / (Ret * Ret);  // = -cf / Ret
    }
}
```

**GetCf turbulent**: chain through bb = log(Ret/Fc). `∂bb/∂Ret = 1/Ret`
(when bb >= 3). All other intermediates depend on Hk only, no Ret. So:
```cpp
// cf0 = 0.3 * exp(aa) * bb^cc
// ∂cf0/∂Ret = 0.3 * exp(aa) * cc * bb^(cc-1) * (1/Ret)
//           = cf0 * cc / (Ret * bb)
// cf = (cf0 + 1.1e-4*(dd-1))/Fc → ∂cf/∂Ret = (∂cf0/∂Ret) / Fc
```

**GetCDi_Lam** (`cDi = num/Ret`): `∂cDi/∂Ret = -cDi/Ret` — same as laminar Cf.

**GetCDi_TurbWall**: log-based, follow chain.

**GetHs turbulent**: explicit log(...) dependence on Ret. Inspect & derive.

**Others**: similar pattern. Each is 1-5 lines of derivative formula.

### Aggregation in ResidualStation

Per station, accumulate `∂R/∂Ret_1` (Vector3d) and `∂R/∂Ret_2` (Vector3d)
where Ret_1, Ret_2 are Ret at the two adjacent BL nodes.

Each appearance of `cf`, `cDi`, `Hs`, `cteq`, `uq`, `us`, `damp`, `Ret`
in the residual formulas contributes:
```
∂R/∂Ret_1 += (∂R/∂cf) * (∂cf/∂Ret * ∂Ret_1/∂Ret_1) * (upwind weight for cf_1)
```

This is essentially a Vector3-valued "second linearization" pass through
the residual formulas. ~150-250 new lines of code in ResidualStation,
mostly mechanical chain rule.

### Then: assembly

In `BuildGlobSys`, after the existing R_U / R_x triplet construction,
add a `dR_dRet` dense vector of length `3*Nsys` (or sparse — only
stations that actually use Ret contribute). For each station residual
contribution to global row `Ig`, add the station's lin_Ret values into
`dR_dRet[Ig:Ig+3]` indexed by the involved station's Ret value
(multiplied by `Ret_station / Re` since Ret depends linearly on Re).

### Then: replace DRVDRe

```cpp
VectorXd DRVDRe(FoilState& state) {
    BuildGlobSys(state);
    return state.glob.dR_dRet;  // already includes ∂Ret/∂Re scaling
}
```

### Estimated effort

- `_dRet` companion functions for 9 closures: ~200 lines
- `ResidualStation` lin_Ret pass: ~150 lines
- `ResidualTransition` lin_Ret pass: ~80 lines
- `BuildGlobSys` aggregation + new `glob.dR_dRet` field: ~50 lines
- New `DRVDRe` analytic + ping test against 4th-order FD: ~50 lines

**Total ~530 lines**, ~4-6 hours focused work plus testing.

## Ma chain (larger scope)

Ma enters R_V through MULTIPLE thermodynamic parameters set by
`InitThermo`:
- `KTb = sqrt(1 - Ma²)`  — Karman-Tsien beta
- `KTl = Ma²/(1+KTb)²`   — Karman-Tsien lambda
- `H0 = (1 + 0.5·gmi·Ma²)·Vinf²/(gmi·Ma²)` — stagnation enthalpy
- `cps = 2/(gam·Ma²)·(ratio^(gam/gmi) - 1)` — sonic Cp
- `mu0 = muinf / finf(Tr, Tsrat)`   — stagnation viscosity (Tr depends on H0)
- `rho0 = rhoinf · (1 + 0.5·gmi·Ma²)^(1/gmi)` — stagnation density

Each closure that uses one or more of these parameters needs its
own `_dMa` companion. Many more closures touch Ma than touch Ret:
- `GetUk` (Karman-Tsien correction)
- `GetMach2`
- `GetCp` (forces output)
- `GetH`, `GetHk`, `GetHs`, `GetHss` (compressibility corrections via Hk)
- `GetRet` (via mu0, rho0)
- All closures that take `param` and use these fields

### Strategy options

1. **Per-closure `_dMa`** (like Re chain, larger): ~1000-1500 lines.
   Each closure needs careful derivative through all 6 Ma-dependent
   params.

2. **Templated closures on scalar type T** (cleanest, most invasive):
   Templated `T = double` or `T = std::complex<double>`. ∂/∂Ma via CS:
   ```cpp
   state.param.Ma = state.oper.Ma + std::complex<double>(0, 1e-30);
   // ... compute all derived params with complex arithmetic ...
   // ... call templated closures + ResidualStation + assembly ...
   ∂R/∂Ma = Im(R) / 1e-30;
   ```
   Requires templating: Param struct, all closures, ResidualStation,
   ResidualTransition, BuildGlobSys, and StagnationState. Plus careful
   handling of `std::abs`, `std::sqrt`, `std::log`, comparisons. ~50%
   of bl_properties.cpp + coupled_solver.cpp touched (signature
   changes), but the body of each function is mostly unchanged. Mech.
   ~2000-3000 lines.

3. **Hybrid**: implement Re chain analytically (small scope) and
   keep DRVDMa as FD with documentation. Re/Ma both still capped at
   ~1e-5 by bug #1 anyway. Re chain demonstrates the pattern; Ma can
   be done once bug #1 is fixed and the demand is clearer.

## Recommendation

Defer until [BUG_rv_ue_istag_column.md](BUG_rv_ue_istag_column.md) is
resolved. That bug gates the entire viscous sensitivity test suite at
1e-5; analytic Re/Ma chains alone won't improve any end-to-end test
result. Once the 1e-5 floor is gone, the precision gains from analytic
Re/Ma become real and the work is more clearly justified.

If "no FD" mandate must be enforced *before* bug #1 is fixed, do the Re
chain first (smaller, contained — proves the pattern). Then decide on
Ma based on whether per-closure or templated approach fits the project
better.
