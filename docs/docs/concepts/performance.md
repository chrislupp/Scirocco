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

- **Paneling** costs less than a millisecond whatever the panel count. It runs
  once per geometry, in the constructor and after each geometry change.
  Changing the operating point doesn't re-panel.
- **Inviscid solves** take a millisecond or two. Sweeping angle of attack
  in potential flow is essentially free.
- **Viscous solves** take about 0.1 s at 99 and 149 panels. Time depends
  mostly on how many Newton steps the case needs. The 199-panel case here is
  a hard one: it takes 53 Newton steps, where the NACA 2412 case in
  [Viscous coupling](./viscous-coupling) takes 8 (see [PTC](./ptc)).
- **The full gradient** covers $C_l$, $C_d$, and $C_m$ with respect to every
  node coordinate, $\alpha$, $Re$, and $M_\infty$. All outputs share one pass:
  the panel-influence derivatives are evaluated once, and the chains through
  the vortex strengths and the wake run in reverse for every output together.
  Cost grows a little faster than the square of the panel count.

## Compared with finite differences

A finite-difference gradient needs one extra viscous solve per input, of
which there are $2N + 3$. At 99 panels that's 203 solves, about 24 s; the
adjoint takes 134 ms. At 149 panels it's about 29 s against 405 ms. And the
adjoint result is exact, with no step size to choose.
