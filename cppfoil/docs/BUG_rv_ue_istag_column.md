# BUG: R_V differs from mfoil at the ue-Istag column (1e-5 gradient floor)

**Status:** open, gating machine-precision viscous gradients.
**Discovered:** 2026-05-22 (see commits 98d3c00, 70bde84).
**Affects:** all viscous sensitivities (shape, alpha, Mach, Re) — they hit
a 1e-5 rel floor against mfoil-CS, even with primal at 1e-13 rel.

## Symptom

`probe_lambda_cs` reports cppfoil's `lambda_cl` differs from mfoil-CS-derived
`lambda_cl` at **3.5e-5 rel** (worst entry: i=414 = sa slot at wake node 2),
while `lambda_cd` differs at 3.6e-6 rel. The 1.19e-5 floor in
`ViscousAlphaSensitivity` (`dcl/dalpha` rel = 1.19e-5) is exactly this
propagating through the lambda chain.

End-to-end test status with this bug:
| Test | Floor |
|---|---|
| `ViscousShapeSensitivity` (dcl/dcd) | ~1e-5 rel |
| `ViscousShapeSensitivity` (dcm)     | ~1e-7 rel (untouched) |
| `ViscousAlphaSensitivity` (dcl/dcd) | ~1e-5 rel |
| `ViscousMachSensitivity` (dcl/dcd)  | ~1e-4 rel (FD chain noise dominates) |
| `ViscousReSensitivity`              | ~1e-3 rel (FD chain noise dominates) |

## What's ruled out (verified to 1e-12 or better)

- **Primal**: cl/cd/cm match mfoil at 1e-13 after the cumulative-foil.s
  fix (commit e7431af).
- **Phase 1 chain partials**: DXi, DWakeS, DWakeX, DWgap, DSstag all pass
  at 1e-10 (`viscous_chain_test`).
- **Phase 2 residual partials**: ∂R_BL/∂x_foil and ∂R_ue/∂x_foil pass at
  1e-12 rel (`viscous_residual_deriv_test`).
- **Output partials**: cl_ue verified against FD on cl(U[3,i]) at
  2.6e-11 abs (`probe_internal_fd`). cm_ue, cd_U also clean
  (`viscous_output_ping_test`).
- **DRVDAlpha**: matches 4th-order FD at 1e-7 rel
  (`viscous_output_ping_test::DRVDAlphaMatchesFD`).
- **Forward/reverse LU consistency**: `b^T (R_V^{-1} g) == lambda^T g` to
  1e-12 (`viscous_output_ping_test::ForwardReverseConsistency`).

## Localization

`probe_rv_compare` shows cppfoil's R_V vs mfoil's R_V at the
`comprehensive_reference` linearization point:

```
Per-row max abs error (top 20 rows): all concentrated at R_BL rows
44–54 (stations adjacent to stagnation Istag=[47,48])

Worst entry: R_V[144, 195] cpp=102.213 mfoil=102.101 abs=0.111835
   row 144 = R_BL momentum at node 48 (Istag[1])
   col 195 = 4*48+3 = U[ue, node 48] (Istag[1])

Max abs error by column variable type:
   col th: max=8.6e-7 at node 48
   col ds: max=7.6e-7 at node 48
   col sa: max=1.1e-12 at node 49
   col ue: max=0.111835 at node 48   <-- THE BUG
```

**The discrepancy is concentrated entirely in the ue columns at Istag.**
This is exactly where the xi-stag-chain augmentation
`R_U[:,4*Istag[k]+3] += R_st * sstag_ue[k]` lands.

## What didn't fix it

1. **Refresh `sstag_ue` from override-U** before BuildGlobSys: cppfoil's
   `sstag_ue` was stale (from cppfoil's own Newton-converged ue) instead of
   reflecting the override `U_ref`. Refreshed via the explicit formula
   `sstag_ue[0] = ue1*(S1-S0)/(den^2)`. Result: R_V[144,195] error reduced
   from 0.111835 → 0.110646. Only 1% of the discrepancy.

2. **Move LU factorization after xi swap** (`StoreAdjointFactorization`
   after `StagPointFind/ExtendXiToWake/SetWakeGap`): no change in gradient
   accuracy.

3. **Remove `if (rnorm > 1e-10)` guard** around `StagPointMove`/
   `UpdateTransition`: no change in primal cl/cd/cm. Newton still reaches
   1.87e-14.

## Likely cause

The xi-stag-chain augmentation in `BuildGlobSys` (`coupled_solver.cpp:417-441`)
or the similarity-station residual lin_U in `ResidualStation` is producing
a slightly different value than mfoil's. The discrepancy is structural
(0.1 abs, not roundoff), and lives in the rows immediately adjacent to
stagnation.

Possible specific causes to investigate:
1. `sstag_ue` formula or which-state-it's-evaluated-at differs in some
   edge case.
2. `StagnationState` linearization (Ust_U) differs at the simi=true
   evaluation path.
3. mfoil's `build_glob_sys` includes some term in lin_U that cppfoil's
   `ResidualStation` doesn't — particularly the simi-branch's xi=xst
   special handling.
4. The `R_st = R_x * x_st` product differs because cppfoil's R_x in the
   override-U scenario doesn't quite match mfoil's R_x at the same state.

## How to reproduce

```bash
cmake --build build --target probe_rv_compare probe_lambda_cs probe_sstag_ue
./build/probe_rv_compare    # diffs R_V entry by entry
./build/probe_lambda_cs     # diffs lambda
./build/probe_sstag_ue      # checks stagnation linearization values
```

`comprehensive_reference.h5` must have `/comprehensive/R_V` — regenerate
with `python3 tests/add_rv_to_reference.py` if needed.

## Next steps to try

1. **Dump ResidualStation lin_U at simi=true for node 48** in both cppfoil
   and mfoil. Diff entry by entry to localize which formula differs.
2. **Test without the xi-stag augmentation**: zero it out and see if the
   rest of R_V matches mfoil. If yes, the bug is in the augmentation. If
   not, the bug is in ResidualStation lin_U.
3. **Replicate mfoil's `R_st += 1e-15` hack**: shouldn't matter but worth
   ruling out as a sparse-storage edge case.
4. **Check that cppfoil's `R_x` at the override-U state matches mfoil's
   `R_x` at the same state row 144**: print row 144 of cppfoil R_x and
   compare to what `dR_BL_dx[144, :, :]` in HDF5 collapses to.

## Acceptance for fix

`probe_rv_compare` should report max abs error < 1e-10 across all R_V
entries. `probe_lambda_cs` should then show max rel < 1e-10 for both
lambda_cl and lambda_cd.
