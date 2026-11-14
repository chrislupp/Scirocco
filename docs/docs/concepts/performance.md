---
title: Performance
description: Measured timings for paneling, inviscid and viscous solves, and gradients.
---

# Performance

These are measured times for one cold start: a new object, paneling, a
solve, and the full gradient. Each phase is timed separately; the table
shows the median of 11 runs.

<BenchTable />

NACA 0012, α = 5°, Re = 10⁶, on an Apple M1 Pro with Apple Clang, a Release
build, and `-march=native`. Run it yourself:

```bash
cmake --build scirocco/build -j --target bench
scirocco/build/bench
```

## Reading the numbers

- **Paneling** costs about a millisecond whatever the panel count. It runs
  once per geometry, in the constructor and after each geometry change.
  Changing the operating point doesn't re-panel.
- **Inviscid solves** take a millisecond or two. Sweeping angle of attack
  in potential flow is essentially free.
- **Viscous solves** take 20-65 ms. Time depends mostly on how many Newton
  steps the case needs; each, with its boundary-layer march, costs about
  4 ms at 200 panels. The 199-panel case here takes 14 Newton steps, where
  the NACA 2412 case in [Viscous coupling](./viscous-coupling) takes 6 (see
  [PTC](./ptc)).
  A second solve on the same object, such as the next point of a polar,
  starts from the previous solution and is usually faster.
- **The full gradient** covers $C_l$, $C_d$, and $C_m$ with respect to every
  node coordinate, $\alpha$, $Re$, and $M_\infty$. All outputs share one pass:
  the panel-influence derivatives are evaluated once, and the chains through
  the vortex strengths and the wake run in reverse for every output together.
  Cost grows a little slower than the square of the panel count.

## Compared with finite differences

A finite-difference gradient needs one extra viscous solve per input, of
which there are $2N + 3$. At 99 panels that's 203 solves, about 24 s; the
adjoint takes 11 ms. At 149 panels it's about 29 s against 21 ms. And the
adjoint result is exact, with no step size to choose.

## Compared with XFOIL

XFOIL, built from the CMPLXFOIL sources with `gfortran -O2`, on the same
machine and coordinates (NACA 0012, Re = 10⁶, 161 input points repaneled by
each code). Cold starts are fresh solves; warm starts follow a solve
0.5° away. Times in milliseconds, Newton steps in parentheses.

| Panels | α | XFOIL cold | Scirocco cold | XFOIL warm | Scirocco warm |
| --- | --- | --- | --- | --- | --- |
| 160 | 2° | 22 (10) | 23 (8) | 12 | 18 |
| 160 | 8° | 18 (8) | 27 (9) | 10 | 25 |
| 200 | 2° | 37 (10) | 39 (10) | 22 | 26 |
| 200 | 5° | 36 (10) | 42 (10) | 20 | 30 |
| 200 | 8° | 31 (8) | 46 (11) | 16 | 27 |
| 284 | 5° | 86 (10) | 89 (10) | 49 | 57 |

Both codes now take about the same number of Newton steps, and Scirocco
matches XFOIL on cold starts at the lower angles. A Scirocco step costs about
20% more at 200 panels, mostly its boundary-layer march, and warm starts
save less than XFOIL's.
XFOIL's tolerance
(`EPS1 = 1e-10` on the RMS relative change) is about as tight as Scirocco's.
