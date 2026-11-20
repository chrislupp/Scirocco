---
title: Transition validation
description: How the envelope onset ramp and the captured transition model were checked against stability theory, wind-tunnel drag and complex step.
---

# Transition validation

The two opt-in transition options on the [Transition](../concepts/transition) page
were checked three ways:

| What | Checked against | Result |
| --- | --- | --- |
| Envelope onset ramp | Orr–Sommerfeld linear stability | RMS 0.02–0.06, against 0.13–0.24 for the default ramp |
| Captured transition, drag | Eppler 387 wind-tunnel drag (UIUC) | RMS Cd error 5.3%, against 16% for the default model |
| Captured transition, gradients | Complex step through mfoil | dCd/dα ripple 1.6%, against 25% |
| Adjoint derivatives | Complex step through mfoil | Shape to 1e-8; α, Re and Mach to 2e-11 |

## Onset ramp against stability theory

The default onset ramp, from XFOIL and mfoil, is a numerical smoothing
device. It isn't fitted to stability data. The envelope ramp is fitted to
spatial Orr–Sommerfeld envelope N-factor curves for nine Falkner–Skan
profiles, with H from 2.36 to 3.48. The solver reproduces two Blasius
benchmarks:

- Jordinson's eigenvalue at Re_δ* = 998;
- the critical Reynolds number, 519.4.

For each profile, the growth rate is scaled by its asymptotic value and
compared to both ramps:

| H | envelope ramp RMS | default ramp RMS |
| --- | --- | --- |
| 2.36 | 0.056 | 0.126 |
| 2.59 (Blasius) | 0.020 | 0.184 |
| 2.96 | 0.019 | 0.213 |
| 3.48 | 0.025 | 0.237 |

In the stability results, growth starts 0.12–0.29 decades of Re_θ below the
correlation's intercept, then rises to its asymptotic rate slowly. The
default ramp starts only 0.1 decades below. Both ramps share the asymptote,
which matches the stability results:

- the intercept agrees within 0.06 decades;
- the slope agrees within 10% for H ≥ 2.48.

On NACA 2412 at α = 4° and Re = 10⁶, the envelope ramp changes Cd by
−0.06% and moves transition 0.0007c. It also smooths the shape gradient at
the onset. The test is a Gaussian bump on the upper surface, 0.01c wide,
swept through the onset:

| | peak positive dCd/dA | peak negative dCd/dA | secondary lobe |
| --- | --- | --- | --- |
| default ramp | +0.77 | −0.64 | +0.10 |
| envelope ramp | +0.29 | −0.15 | none |

## Captured transition against measured drag

### Test case

The Eppler 387 is the standard low-Reynolds-number transition test. Every
condition below has a laminar separation bubble. That makes it a hard test
of the transition model, because the bubble's length sets the drag.

- **Reference:** UIUC drag polars from Williamson et al., *Summary of
  Low-Speed Airfoil Data, Vol. 5* (2012), at Re = 1, 2, 3 and 4.6×10⁵.
- **Comparison:** Cd at matched Cl over 0.1 ≤ Cl ≤ 1.0, the attached-flow
  part of the polar, interpolating the measured polar at each computed Cl.
- **Points:** 50 at ncrit = 9.

<E387DragPlot />

### Results

| Re | default (cut cell) | captured, L = 40θ | **L = 60θ** | L = 80θ | L = 120θ | L = 240θ |
| --- | --- | --- | --- | --- | --- | --- |
| 1×10⁵ | −28.4 / 29.0 | −4.9 / 7.1 | −4.0 / 6.3 | −3.9 / 6.2 | −3.9 / 6.2 | −3.8 / 6.2 |
| 2×10⁵ | −10.0 / 10.2 | −3.5 / 3.7 | −2.5 / 3.1 | −2.7 / 3.5 | −2.9 / 4.0 | −3.2 / 4.7 |
| 3×10⁵ | −7.3 / 8.1 | −3.9 / 5.0 | −3.7 / 5.2 | −4.0 / 5.7 | −4.5 / 6.3 | −5.5 / 7.7 |
| 4.6×10⁵ | −6.7 / 7.6 | −4.5 / 5.5 | −4.8 / 6.0 | −5.3 / 6.6 | −6.1 / 7.3 | −7.8 / 9.3 |
| all, ncrit = 9 | −12.8 / 16.0 | −4.2 / 5.4 | **−3.7 / 5.3** | −4.0 / 5.6 | −4.3 / 6.1 | −5.1 / 7.2 |
| all, ncrit = 8 | −14.3 / 17.8 | −7.3 / 9.3 | **−6.7 / 8.5** | −6.9 / 8.6 | −7.5 / 9.1 | |

Each entry is the Cd error in percent, as mean signed / RMS.

- **The default model under-predicts bubble drag.** It makes the flow
  turbulent at a single point, so the bubble reattaches too early. At
  Re = 10⁵ it under-predicts drag by 28%.
- **The captured model is about three times more accurate.** It spreads
  transition over a finite length, so the bubble grows to a realistic size.
- **L = 60θ is the calibrated default** (`capt_L`):
  - it is the best length at ncrit = 9, and within 0.1% of the best at
    ncrit = 8;
  - the error is flat from 40θ to 80θ;
  - longer lengths keep the boundary layer laminar too long at the higher
    Reynolds numbers.
- **Remaining bias:** Cd is 3–4% low at every length. That points at the
  rest of the closure set, or at ncrit, not at the transition length.
- **Convergence:** both models converge at 13–16 of the 23 angles from −3°
  to 8°. The failures are at the ends of the polar, and the captured model
  adds none of its own.

The UIUC tunnel's turbulence intensity is about 0.1%. Mack's
correlation puts ncrit at about 8 for that level, so both ncrit values are
plausible. ncrit = 9 is closer to the measurements with either model.

### Setup and caveats

- **Panels and trailing edge:** 160 panels, with the trailing edge opened to
  0.25% chord, because an exactly closed trailing edge fails (see
  [Troubleshooting](../getting-started/troubleshooting#the-solve-doesnt-converge)).
- **Panel limit:** above 160 panels, the coarse UIUC coordinates (61
  points) break the boundary-layer initialization at the leading edge. This
  happens in mfoil too.
- **Scope:** one airfoil from one tunnel. Measured transition locations,
  such as the E387 oil-flow data in NASA TM-4062, would be an independent
  second check.

## Captured transition against the default model

On NACA 0012, 2412 and 4412 at Re = 2×10⁵, 10⁶ and 3×10⁶, and α = 0, 4
and 8° (160 panels):

- the captured model converges in 26 of 27 cases. The failure (2412,
  Re = 3×10⁶, 8°) also fails with the default model.
- Where the default model converges, the captured Cd is 0.2–3.7% higher at
  Re ≥ 10⁶. At Re = 2×10⁵, where bubbles form, it ranges from 0.5% lower
  to 8.3% higher. The Eppler 387 comparison above shows that the captured
  value is the more accurate one.

## Gradient smoothness

Smoothness is measured by the exact (complex-step) dCd/dα over 21 angles
from 3.9° to 4.3°. Over that range dCd/dα should vary smoothly, so ripple
is the RMS deviation from a quadratic fit, relative to the mean.

<TransitionSweepPlot />

| case | panels | default (cut cell) | captured, L = 60θ |
| --- | --- | --- | --- |
| NACA 2412, Re = 10⁶ | 160 | 24.6% (max 47.5%) | 1.6% (max 2.9%) |
| NACA 0012, Re = 2×10⁵ (bubble) | 160 | 29.6% (max 51.1%) | 0.2% (max 0.5%) |
| NACA 4412, Re = 3×10⁶ | 160 | 21.2% (max 43.1%) | 7.6% (max 15.2%) |
| NACA 4412, Re = 3×10⁶ | 240 | 20.4% (max 35.9%) | 3.7% (max 6.3%) |
| NACA 4412, Re = 3×10⁶ | 320 | 85% (max 359%)¹ | 1.3% (max 2.1%) |

¹ The sweep doesn't check cut-cell convergence, so some of these points may
not have converged.

- **The captured gradient converges as the panels are refined; the cut-cell
  gradient doesn't.**
- **The remaining ripple is a resolution effect.** The transition zone has
  to span several panels. At high Reynolds numbers θ is small, so 60θ
  covers fewer panels. On NACA 4412 at Re = 3×10⁶, use 240 panels or more
  for smooth gradients.
- **Smoothness against accuracy:** a longer transition length smooths the
  gradient on coarse panels, but it is less accurate against the Eppler 387
  drag. At L = 80θ the 4412 ripple is 4.3% at 160 panels.

## Derivatives against complex step

Every derivative is checked against complex step through the mfoil
prototype of the same model. That reference is exact to machine precision.
The tests are in `scirocco/tests/sensitivity_reference_test.cpp`.

| Quantity | Case | Agreement |
| --- | --- | --- |
| Envelope onset: shape | NACA 2412, 160 panels, α = 4°, Re = 10⁶ | 1e-8 |
| Envelope onset: α, Re, Mach | same | 1e-9 |
| Default onset: Re, lower surface laminar to the trailing edge | same | 1e-9 |
| Captured: cl, cd, cm | NACA 2412 at Re = 2×10⁵, 10⁶, 3×10⁶ | 1e-10 |
| Captured: shape, 9 nodes across onset and the blend | NACA 2412, Re = 10⁶ | 1e-8 |
| Captured: α, Re, Mach | same | 2e-11 |

## Reproducing

The scripts and computed data are in
[`scirocco/validation/onset_envelope`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/onset_envelope)
and
[`scirocco/validation/captured_transition`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/captured_transition).

`e387_calibrate.py` downloads the UIUC data when it runs. The measurements
aren't stored in the repository, because UIUC's terms require their notices
to accompany any copy. The figures on this page show only the derived
error, which `e387_calibrate.py docs` writes.
