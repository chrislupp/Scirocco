# Viscous sensitivity status (2026-05-23)

## Final test results

| Sensitivity | dcl rel | dcd rel | dcm rel | Status |
|---|---|---|---|---|
| `ViscousShape` (foil.x) | 1.95e-11 | 6.18e-12 | 1.27e-11 | **machine precision** ✓ |
| `ViscousAlpha` (α) | 7.01e-12 | 1.01e-10 | 3.97e-11 | **machine precision** ✓ |
| `ViscousMach` (Ma) | 8.65e-05 | 2.07e-08 | 4.25e-05 | FD-limited |
| `ViscousRe` (Re) | 8.29e-11 | 6.23e-11 | 3.20e-11 | **machine precision** ✓ |

| Validation test | Status |
|---|---|
| `viscous_chain_test` (Phase 1 geom chains) | 0 failures, ~1e-13 max abs |
| `viscous_residual_deriv_test` (Phase 2 ∂R/∂x) | 0 failures, max abs 2.5e-12 |

## Completed work (this session)

1. **R_V ue-Istag column bug fix** (commit 4051c2d): removed post-Newton
   `StagPointFind` swap in `SolveViscous`. Production now uses ue-based xi
   throughout. Resolved the 0.11 abs error in R_V at the xi-stag-chain
   augmentation column.
2. **Analytic chain → ue-based** (commit b7359c4): `DSstagDAlpha`,
   `DSstagDFoilX` switched from gamma-based to ue-based formulas, matching
   the new production state. Shape and Alpha jumped from 1e-5 to machine
   precision.
3. **Analytic Re chain** (commits 2c28b00, 09262f1):
   - `Vector3d lin_Ret` field in `ResidualLin`
   - `ComputeStationLinRet` mirrors `ResidualStation` for the Ret chain
   - 13 `_dRet` helpers (laminar Cf, CDi sub-functions, Hs, Us, Hk, Cfxt,
     CDixt, Cteq, Damp, Uq via direct Ret arg + cf chain)
   - `ResidualTransition` extended with xt(Re) implicit FN chain
   - `BuildGlobSys` accumulator `glob.dR_dRe_times_Re`
   - `DRVDRe` analytic via `dR_dRe_times_Re / Re`
4. **Comprehensive reference regen at ue-based xi** (commit 35ebcf0):
   `extract_comprehensive_reference.py` now calls `stagpoint_move`
   (not `stagpoint_find`) after `recompute_geom`, matching cppfoil's
   production state. `viscous_residual_deriv_test` fixture also refreshes
   sstag via `StagPointMove(s)` after the U override.
5. **Ma chain infrastructure** (commit 555e658): full scaffold in place
   (Vector3d lin_Ma, dR_dMa_times_Ma, ComputeStationLinMa, per-closure
   _dMa helpers) but not yet correct end-to-end — `DRVDMa` falls back to
   4th-order FD until the analytic chain is debugged.

## Remaining work: Ma chain debugging

The Ma chain infrastructure is fully laid down but `probe_drv_dma` shows:
  - Worst entry: station 71 R_lag (upper turbulent) — analytic 0.715 vs FD -0.014
  - Laminar Rshape mismatches at ~3x off

Likely bug locations (in `ComputeStationLinMa`, `cppfoil/src/boundary_layer.cpp`):
  - The turbulent Rlag chain — Clag/Us, cteq, uq Ma-chains have many paths;
    likely an arithmetic error in one of them
  - Possible double-counting in `cfxtu_dMa` (upwind expansion has the cross-product
    `-upw_dMa*cfxt1 + upw_dMa*cfxt2` which should be `upw_dMa*(cfxt2 - cfxt1)`)
  - Direct ∂J/∂Ma analytic (Karman-Tsien + Squire-Young) also needs proper
    debugging before being put back in `ViscousSensitivities`

To finish Ma chain:
  1. Use `probe_drv_dma` to localize bugs station-by-station, equation-by-equation
  2. Compare against FD on each closure individually (write `probe_get_X_dMa`)
  3. Fix bugs in `ComputeStationLinMa`
  4. Re-enable analytic `DRVDMa` (already a one-line change in `adjoint.cpp`)
  5. Re-enable analytic `dJ_dMa` in `ViscousSensitivities` (the FD block
     near line 632 of adjoint.cpp)

Estimated effort: ~6-10 hours focused debugging.

## Why ViscousMach hasn't improved despite all the lambda-bug fixes

ViscousMach was at ~1e-5 rel before any of this session's work, and stays at
~1e-5 now. With `DRVDMa` on FD and `dJ_dMa` on FD, ViscousMach precision is
capped by FD truncation/cancellation on the Ma direction. Once the analytic
Ma chain is debugged and enabled, ViscousMach should drop to machine
precision (1e-11 or so), matching Shape/Alpha/Re.
