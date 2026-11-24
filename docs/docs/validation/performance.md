---
title: Performance
description: Measured timings for paneling, inviscid and viscous solves, and gradients, against XFOIL, flexfoil and CMPLXFOIL.
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
- **Viscous solves** take 12-40 ms. Time depends mostly on how many Newton
  steps the case needs; each, with its boundary-layer march, costs about
  2.5 ms at 200 panels. The 199-panel case here takes 14 Newton steps, where
  the NACA 2412 case in [Viscous coupling](../concepts/viscous-coupling) takes 6 (see
  [PTC](../concepts/ptc)). A second solve on the same object, such as the next point
  of a polar, starts from the previous solution and reuses everything that
  depends only on the geometry: a 0.5° step takes 8-17 ms at 200 panels.
- **The full gradient** covers $C_l$, $C_d$, and $C_m$ with respect to every
  node coordinate, $\alpha$, $Re$, and $M_\infty$. All outputs share one pass:
  the panel-influence derivatives are evaluated once, and the chains through
  the vortex strengths and the wake run in reverse for every output together.
  Cost grows a little slower than the square of the panel count.

## Compared to finite differences

A finite-difference gradient needs one extra viscous solve per input, of
which there are $2N + 3$. At 99 panels that's 203 solves, about 2.5 s even
with warm starts; the adjoint takes 9 ms. At 149 panels it's about 7 s
against 18 ms. And the adjoint result is exact, with no step size to choose.

## Compared to XFOIL

These are viscous analyses alone, without gradients. XFOIL here is the
Fortran from the CMPLXFOIL sources, built with `gfortran -O2
-fdefault-real-8` and called through a small Fortran driver, so there is no
Python overhead. Both codes ran on the same machine with the same
coordinates: NACA 0012 at Re = 10⁶, 321 input points, repaneled by each code
to the panel count shown. Cold starts are fresh solves. Warm starts follow
a solve 0.5° away. Times are in milliseconds, with Newton steps in
parentheses.

| Panels | α | XFOIL cold | Scirocco cold | XFOIL warm | Scirocco warm |
| --- | --- | --- | --- | --- | --- |
| 160 | 2° | 22 (10) | 17 (8) | 12 | 8 |
| 160 | 8° | 18 (8) | 20 (9) | 10 | 11 |
| 200 | 2° | 37 (10) | 28 (10) | 22 | 11 |
| 200 | 5° | 36 (10) | 29 (10) | 20 | 17 |
| 200 | 8° | 31 (8) | 31 (11) | 16 | 13 |
| 284 | 5° | 86 (10) | 53 (10) | 49 | 26 |

Both codes take about the same number of Newton steps. Scirocco is faster
on most of these cases, and the margin grows with the panel count: its
Newton step costs O(N²) plus one dense N × N factorization done with BLAS.
Warm starts reuse the geometry-only work and predict the new edge velocity,
so most take 5-6 steps. XFOIL's convergence tolerance (`EPS1 = 1e-10` on the
RMS relative change) is about as tight as Scirocco's.

## Compared to flexfoil

[flexfoil](./flexfoil) is a Rust rewrite of XFOIL. These are viscous
analyses on the same nodes: mfoil panels each airfoil, Scirocco takes the
nodes unchanged, and flexfoil 1.1.6 builds its panels from them, so both
codes solve the same discrete problem. flexfoil runs through its Python
package, single-threaded. Times are medians of 11 runs in milliseconds,
with Newton steps in parentheses. Warm starts follow a Scirocco solve 0.5°
away; flexfoil's interface has no warm start.

| Case | Panels | flexfoil | Scirocco cold | Scirocco warm |
| --- | --- | --- | --- | --- |
| NACA 0012, Re = 10⁶, α = 2° | 159 | 33 (7) | 17 (8) | 7.5 |
| NACA 0012, Re = 10⁶, α = 8° | 159 | 28 (6) | 15 (6) | 9.0 |
| NACA 0012, Re = 10⁶, α = 2° | 199 | 53 (7) | 24 (8) | 12 |
| NACA 0012, Re = 10⁶, α = 5° | 199 | 70 (11) | 37 (15) | 16 |
| NACA 0012, Re = 10⁶, α = 8° | 199 | 39 (5) | 48 (15) | 12 |
| NACA 0012, Re = 10⁶, α = 5° | 283 | 156 (12) | 53 (10) | 26 |
| NACA 2412, Re = 3×10⁶, α = 4° | 199 | 68 (10) | 26 (9) | 12 |
| NACA 4412, Re = 2×10⁵, α = 4° | 199 | 44 (6) | 25 (9) | 13 |

From a cold start Scirocco is about twice as fast at 160–200 panels and
three times as fast at 283. flexfoil takes as many Newton steps or fewer;
each costs about 5.5 ms at 199 panels, against about 2 ms for Scirocco, and
its cost grows faster with the panel count. The one case flexfoil wins,
NACA 0012 at 8°, is one where Scirocco needs 15 steps to flexfoil's 5.

A whole polar widens the gap, because Scirocco warm-starts each point
from the last. NACA 0012 at Re = 10⁶ with 199 panels, α = 0…15°:

| | Time | Converged |
| --- | --- | --- |
| flexfoil, one point after another | 1.33 s | 15 / 16 |
| flexfoil, `analyze_faithful_batch` on 10 cores | 0.66 s | 15 / 16 |
| Scirocco, one sweep, including process start | 0.28 s | 16 / 16 |

flexfoil's batch solves points in parallel, but the point at 1° runs all
100 iterations without converging (the limit cycle XFOIL also has, see
[flexfoil](./flexfoil#viscous-polars)) and sets the wall-clock time.

The PyPI wheel is built for size (`opt-level = "s"`). The times above use
flexfoil built from the same source with `opt-level = 3` and
`target-cpu=native`, to match Scirocco's `-march=native`; that is only
about 5% faster than the wheel. The scripts are in
[`scirocco/validation/vv`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/vv)
(`flexfoil_timing.py`).

## Compared to CMPLXFOIL

In a gradient-based optimization every evaluation is a solve plus a
gradient. [CMPLXFOIL](https://github.com/mdolab/CMPLXFOIL) computes its
gradient by forward complex step, which costs one complex-valued XFOIL solve
per design variable. Scirocco's adjoint costs one solve however many design
variables there are: it gives $C_l$, $C_d$, and $C_m$ with respect to every
node coordinate, $\alpha$, $Re$, and $M_\infty$, and any shape
parameterization is a chain-rule product on top.

**Setup.** CMPLXFOIL 2.1.1 through its Python API, with a fresh geometry per
evaluation (`setCoordinates`, then a solve, then one
`setCoordinatesComplex` and complex solve per design variable, as
`evalFunctionsSens` does). The case is NACA 0012 at α = 4° and Re = 10⁶,
from 160 input points that CMPLXFOIL repanels to 284 nodes. Scirocco follows
its [Philote discipline](../python/servers): it writes the node
coordinates, clears the solution, solves, and calls
`ComputeSensitivities()`. Times are medians over 8 nearby geometries.

### CST shape variables

These use 4 CST coefficients per surface plus α, 9 design variables in all,
as in CMPLXFOIL's own optimization example:

| | Function | Gradient | Function + gradient |
| --- | --- | --- | --- |
| CMPLXFOIL, complex step | 94 ms | 1.43 s | 1.53 s |
| CMPLXFOIL, finite differences | 96 ms | 0.85 s | 0.95 s |
| Scirocco, 285 nodes | 38 ms | 0.21 s | 0.24 s |
| Scirocco, 201 nodes | 22 ms | 0.09 s | 0.11 s |

At the same resolution Scirocco is about 6x faster than CMPLXFOIL's complex
step. Each additional design variable costs CMPLXFOIL one more complex
solve, about 160 ms, and costs Scirocco nothing.

### Full shape gradient

These are derivatives with respect to every coordinate, for one geometry:

| | Gradient |
| --- | --- |
| CMPLXFOIL, complex step (321 solves: x and y of 160 points, plus α) | 51.4 s |
| Scirocco, 285 nodes (all node coordinates, α, Re, Ma) | 0.21 s |
| Scirocco, 201 nodes | 0.09 s |

That is about 250x faster at the same resolution.

### Notes

- **Robustness.** One of the 8 geometries did not converge in CMPLXFOIL,
  and neither did any of its gradient solves. Scirocco converged on all 8.
- **Making CMPLXFOIL's complex step work.** Built with gfortran 15, NumPy 2,
  and f2py's meson backend, CMPLXFOIL needed three fixes:
  - Add `-fdefault-double-8` to the f2py step of the complex build. Without
    it, `-fdefault-real-8` widens the wrapper's `double complex` to 32 bytes,
    so Python reads and writes the Fortran common blocks at the wrong
    offsets.
  - Pass `h=1e-30j` to `evalFunctionsSens(..., mode="CS")`. The default
    step, `1e-200j`, underflows inside the boundary-layer code and returns
    a derivative close to the inviscid one.
  - Replace `np.NaN` with `np.nan` for NumPy 2.

  `AeroProblem` also derives Re from Mach and altitude, so the CMPLXFOIL
  cases ran at M = 10⁻⁴ rather than 0.

All timings are on an Apple M1 Pro, with Release builds, `-march=native`, and
BLAS/LAPACK from the Apple Accelerate framework.
