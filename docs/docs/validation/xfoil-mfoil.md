---
title: XFOIL and mfoil
description: Scirocco against XFOIL 6.99 and mfoil on identical panel nodes, inviscid and viscous, from low Reynolds number to M = 0.4.
---

# Comparison with XFOIL and mfoil

Scirocco implements the model of [mfoil](https://websites.umich.edu/~kfid/codes/mfoil/),
Fidkowski's Python reimplementation of XFOIL's panel method and integral
boundary layer. Wherever both converge, Scirocco with mfoil's closure
selected reproduces mfoil's forces and transition locations to every
stored digit. Both differ from XFOIL by the same small amount: on average
0.3–0.5% in Cl and 0.5–1% in Cd.

| Case | Converged (Scirocco / XFOIL / mfoil) | Cl vs XFOIL, mean / max | Cd vs XFOIL, mean / max | Transition vs XFOIL, max | vs mfoil |
| --- | --- | --- | --- | --- | --- |
| NACA 2412, inviscid, α = −4…12° | 9 / 9 / 9 | 0.00% / 0.00% | — | — | identical |
| NACA 0012, Re = 10⁶, α = 0…15° | 16 / 15 / 12 | 0.31% / 0.89% | 0.94% / 3.8% | 0.0014c | identical |
| NACA 2412, Re = 3×10⁶, α = −4…14° | 19 / 18 / 14 | 0.50% / 1.6% | 0.54% / 1.0% | 0.0023c | identical |
| NACA 4412, Re = 2×10⁵, α = −2…12° | 15 / 15 / 13 | 0.35% / 0.85% | 0.84% / 2.8% | 0.0082c | identical |
| NACA 2412, Re = 6×10⁶, M = 0.4, α = −2…6° | 8 / 9 / 7 | 0.47% / 0.68% | 0.48% / 0.77% | 0.0024c | identical |

Differences are relative, over the angles where both codes converged; Cl
differences skip |Cl| < 0.05. "Identical" means equal to the six
significant digits stored, which leaves a largest transition difference of
1e-4c from rounding.

## Setup

- **Same nodes.** mfoil panels each airfoil with 199 panels. Scirocco takes
  those nodes, arc lengths and tangents unchanged, and XFOIL loads the
  nodes without repaneling (`PCOP`). Paneling differences are therefore
  ruled out; what remains is the solver.
- **Same physics.** ncrit = 9, free transition on both surfaces, Kármán–Tsien
  correction at M = 0.4. Scirocco runs mfoil's closure
  (`cf_laminar_floor = false`); see
  [The skin-friction floor](#the-skin-friction-floor) for its default.
- **How each code runs.** Scirocco and XFOIL sweep each polar the way a user
  would, starting at the angle nearest 0° and continuing up, then down, each
  point warm-starting from the last (Scirocco with `warm_start = true`; see
  [warm starts](../concepts/viscous-coupling#warm-starts)). mfoil solves every point from a cold
  start. XFOIL is version 6.99 with 100 iterations allowed per point.
- **Precision.** XFOIL's forces and distributions are read from its binary
  polar dump, in double precision, not from the rounded text polar.

## Inviscid

On the same nodes the three linear-vorticity panel methods give the same
lift and moment to five or more significant digits. The inviscid Cp is
the starting point for everything else, so the viscous differences below
don't come from the panel method.

<CodeDistributionPlot name="inviscid_2412" quantity="cp" />

## Viscous polars

<CodePolarPlot name="2412_re3e6" />

The three curves lie on top of each other until stall approaches. The
largest XFOIL differences, 1–4%, come at the highest angles.

<CodePolarPlot name="0012_re1e6" />

At low Reynolds number a laminar separation bubble forms ahead of
transition. The codes agree on it as closely as at the higher Reynolds
numbers:

<CodePolarPlot name="4412_re2e5" />

At M = 0.4 the Kármán–Tsien correction matches XFOIL's to the same level
as the incompressible cases. Scirocco and mfoil both fail to converge at
4°, from the same model, where XFOIL converges.

<CodePolarPlot name="2412_re6e6_m04" />

## Transition

The transition point is where the codes could disagree most, because it
depends on the amplification-factor march and on where in the interval
the cut cell places it. Scirocco and XFOIL put it within 0.0025c of each
other at Re ≥ 10⁶, and within 0.01c at Re = 2×10⁵:

<CodeTransitionPlot name="4412_re2e5" />

## Surface distributions

NACA 2412 at Re = 3×10⁶ and α = 4°. The markers are every sixth mfoil node.

<CodeDistributionPlot name="2412_re3e6" quantity="cp" />

<CodeDistributionPlot name="2412_re3e6" quantity="cf" surface="upper" yMax={0.01} />

Away from transition the skin friction from Scirocco and XFOIL agrees to
about 1e-6 (the median difference). The largest differences, up to 6e-4,
are in the first few turbulent stations, where XFOIL's Cf rises a little
faster. The pressure differs most near the trailing edge, by up to 0.008 in
Cp.

NACA 4412 at Re = 2×10⁵ and α = 4°: Cf drops below zero on the upper
surface from x/c = 0.39 to 0.59. That is the laminar separation bubble;
transition at 0.58c closes it.

<CodeDistributionPlot name="4412_re2e5" quantity="cf" surface="upper" yMax={0.01} />

## The skin-friction floor

By default Scirocco does not let the turbulent skin friction fall below the
laminar value, as XFOIL does and mfoil does not (see
[Trips near the leading edge](../concepts/transition#trips-near-the-leading-edge)).
The comparison above switches that off, to check the implementation against
mfoil. With the default, the same polars converge at the same points and
move slightly towards XFOIL in drag:

| Case | Cd change, mean / max | Cl change, max | Cd vs XFOIL, mean / max |
| --- | --- | --- | --- |
| NACA 0012, Re = 10⁶ | 0.08% / 0.19% | 0.11% | 0.86% / 3.6% |
| NACA 2412, Re = 3×10⁶ | 0.02% / 0.06% | 0.02% | 0.53% / 0.94% |
| NACA 4412, Re = 2×10⁵ | 0.23% / 0.75% | 0.14% | 0.62% / 2.4% |
| NACA 2412, Re = 6×10⁶, M = 0.4 | 0.0005% / 0.004% | 0.00% | 0.48% / 0.77% |

The floor acts at the first few turbulent stations behind transition. Its
effect is largest at the lowest Reynolds number, where a laminar separation
bubble precedes transition.

## Convergence

Scirocco converged at 58 of the 59 viscous points, XFOIL at 57, and mfoil,
from cold starts, at 46. XFOIL's miss at 1° on NACA 0012 is a limit cycle:
its Newton iteration settles on Cl = 0.1075 and stalls there for all 100
iterations. mfoil's misses are all at α ≥ 4°, each solved from a cold
start.

## Reproducing

The scripts are in
[`scirocco/validation/vv`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/vv).
`vv.py codes` runs all three codes, and [flexfoil](./flexfoil), and writes
`docs/src/data/vv_codes.json`, which the figures on this page read. It needs
XFOIL on the `PATH` and the `flexfoil` Python package, and takes
about 80 seconds on an Apple M1 Pro; the mfoil solves dominate.
