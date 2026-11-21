---
title: Verification and validation
description: How Scirocco is checked against mfoil, XFOIL and wind-tunnel data.
---

# Verification and validation

**Verification** asks whether Scirocco solves its equations correctly.
**Validation** asks whether those equations describe real airfoils. They
need different references:

| Question | Reference | Result | Page |
| --- | --- | --- | --- |
| Does Scirocco implement its model correctly? | mfoil, the Python code whose model Scirocco implements | Identical forces and transition locations wherever both converge | [XFOIL and mfoil](./xfoil-mfoil) |
| How close is it to the standard tool? | XFOIL 6.99 | Cl within 0.3–0.5% on average, Cd within 0.5–1% | [XFOIL and mfoil](./xfoil-mfoil) |
| Does it predict measured flow? | Published wind-tunnel data | Lift slope 3.5% high, Cd 9–17% low, Cl,max 8% high | [Wind-tunnel data](./experiments) |
| Are the transition options accurate? | Stability theory, Eppler 387 drag | Captured transition is about three times closer to measured drag | [Transition validation](./transition) |
| How fast is it? | XFOIL, CMPLXFOIL, finite differences | Faster viscous solves at 200+ panels; gradients 6x faster than complex step | [Performance](./performance) |

## Test suite

The C++ tests run on every pull request. `scirocco_tests` checks geometry,
splines, panel influences, boundary-layer closures, residuals and the Newton
solve function by function against reference data generated with mfoil,
stored as HDF5 files in `scirocco/tests/`. Most functions agree to 1e-10.
[Contributing](../contributing) covers how to build and run the tests.

## What isn't covered

- **Forced transition.** Scirocco has no trip strip (see
  [Transition limits](../concepts/transition#limits)), so tripped
  measurements are compared only where the trip barely matters, such as
  lift.
- **Separated flow.** Near and past stall, the codes agree with each other
  but not with measurements. See [Wind-tunnel data](./experiments).
- **Transonic flow.** The compressibility correction doesn't capture
  shocks (see [Compressibility](../concepts/compressibility#limits)). The
  compressible case here stays subsonic.
