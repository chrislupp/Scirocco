---
title: flexfoil
description: Scirocco against flexfoil, Flexcompute's Rust rewrite of XFOIL, on identical panel nodes.
---

# Comparison with flexfoil

[flexfoil](https://github.com/flexcompute/flexfoil) is Flexcompute's
rewrite of XFOIL in Rust, available as a Python package and in the
browser. It is a second, independent port of XFOIL's panel method and
integral boundary layer, written in a compiled language like Scirocco. On
the same nodes, flexfoil 1.1.6 reproduces XFOIL 6.99 to about 0.01% in Cl
and 0.001% in Cd, and fails to converge at the same angles. Scirocco
therefore differs from flexfoil by the same amount as from XFOIL: on
average 0.3–0.5% in Cl and 0.5–1% in Cd.

| Case | Converged (Scirocco / flexfoil / XFOIL) | Cl vs flexfoil, mean / max | Cd vs flexfoil, mean / max | Transition vs flexfoil, max | flexfoil vs XFOIL, Cl / Cd max |
| --- | --- | --- | --- | --- | --- |
| NACA 2412, inviscid, α = −4…12° | 9 / 9 / 9 | 0.005% / 0.012% | — | — | 0.011% / — |
| NACA 0012, Re = 10⁶, α = 0…15° | 16 / 15 / 15 | 0.31% / 0.88% | 0.94% / 3.8% | 0.0066c | 0.009% / 0.001% |
| NACA 2412, Re = 3×10⁶, α = −4…14° | 19 / 18 / 18 | 0.49% / 1.6% | 0.54% / 1.0% | 0.0026c | 0.015% / 0.003% |
| NACA 4412, Re = 2×10⁵, α = −2…12° | 15 / 15 / 15 | 0.35% / 0.84% | 0.84% / 2.8% | 0.0082c | 0.006% / 0.004% |
| NACA 2412, Re = 6×10⁶, M = 0.4, α = −2…6° | 8 / 0 / 9 | — | — | — | — |

Differences are relative, over the angles where both codes converged; Cl
differences skip |Cl| < 0.05. The cases and nodes are those of
[XFOIL and mfoil](./xfoil-mfoil), so the two pages can be read side by side.

## Setup

- **Same nodes.** mfoil panels each airfoil with 199 panels. Scirocco takes
  those nodes unchanged, XFOIL loads them without repaneling (`PCOP`), and
  flexfoil builds its panels from the coordinates it is given, without
  repaneling.
- **Same physics.** ncrit = 9, free transition on both surfaces, at most
  100 Newton iterations per point. Scirocco runs mfoil's closure, as on the
  [XFOIL and mfoil](./xfoil-mfoil#the-skin-friction-floor) page.
- **How each code runs.** Scirocco and XFOIL sweep each polar from the
  angle nearest 0°, each point warm-starting from the last (Scirocco with
  `warm_start = true`). flexfoil's
  Python interface has no warm start, so it solves every point from a cold
  start, like mfoil.
- **Version.** flexfoil 1.1.6 from PyPI (`pip install flexfoil==1.1.6`),
  through its `analyze_faithful` and `get_bl_distribution` functions.

## Viscous polars

<CodePolarPlot name="2412_re3e6" third="flexfoil" />

flexfoil's markers sit on the XFOIL curve throughout, including the miss at
6°, where both stop short of convergence.

<CodePolarPlot name="0012_re1e6" third="flexfoil" />

XFOIL's limit cycle at 1° on NACA 0012 (see
[Convergence](./xfoil-mfoil#convergence)) recurs in flexfoil: both settle
on Cl = 0.1075 and stop there after 100 iterations, which suggests that
flexfoil follows XFOIL's Newton iteration step by step. Scirocco converges at that point.

<CodePolarPlot name="4412_re2e5" third="flexfoil" />

## Transition

<CodeTransitionPlot name="4412_re2e5" third="flexfoil" />

flexfoil and XFOIL put transition within 0.0001c of each other in this
case. The only differences above 0.0005c anywhere are two points, NACA 0012
at 6° and NACA 2412 at 5°, where XFOIL places lower-surface transition
within 0.006c of the trailing edge and flexfoil reports none (x/c = 1).
That edge case sets the 0.0066c maximum against Scirocco for NACA 0012;
Scirocco and XFOIL agree there to 0.0014c.

## Surface distributions

NACA 2412 at Re = 3×10⁶ and α = 4°. The markers are every sixth flexfoil
node.

<CodeDistributionPlot name="2412_re3e6" quantity="cp" third="flexfoil" />

<CodeDistributionPlot name="2412_re3e6" quantity="cf" surface="upper" yMax={0.01} third="flexfoil" />

flexfoil's Cp matches XFOIL's to 4e-5 and its Cf to 6e-6 at every node,
in all three viscous cases, so the Scirocco–flexfoil differences are those
described for XFOIL in [Surface distributions](./xfoil-mfoil#surface-distributions).

## Compressible flow

flexfoil 1.1.6 converged at none of the nine points at M = 0.4, nor at
M = 0.1, within 100 iterations; at M = 0 the same points converge in at
most nine. Its unconverged Cl stays at the incompressible value instead of rising
by the Kármán–Tsien factor, about 9% at M = 0.4. No compressible comparison
is possible with this version.

## Inviscid

flexfoil's inviscid routine (`analyze_inviscid`) is separate from its
viscous solver. Its Cl agrees with Scirocco's and XFOIL's to 0.012%, but
its Cm drifts away from theirs as lift grows: by 0.003 at 4° and 0.03 at
12°. Its viscous solver doesn't share the discrepancy; its Cm matches
XFOIL's to 5e-6.

<CodeDistributionPlot name="inviscid_2412" quantity="cp" third="flexfoil" />

flexfoil reports inviscid Cp at panel midpoints rather than at the nodes.

## Speed

Scirocco solves a single point about twice as fast as flexfoil at 200
panels, and a warm-started polar 2.4 to 4.7 times as fast, depending on
whether flexfoil solves the points in parallel; see
[Performance](./performance#compared-to-flexfoil).

## Reproducing

[`scirocco/validation/vv`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/vv)
holds the scripts. `vv.py codes` runs flexfoil alongside Scirocco, XFOIL
and mfoil and writes all four into `docs/src/data/vv_codes.json`; it needs
the `flexfoil` Python package as well as XFOIL. flexfoil adds a few
seconds; the cold-start mfoil solves dominate the run time.

flexfoil's boundary-layer output normalizes Cf by the local edge velocity,
so `vv.py` multiplies it by ue² to get XFOIL's freestream-based Cf, and
computes Cp from ue with XFOIL's Kármán–Tsien formula.
