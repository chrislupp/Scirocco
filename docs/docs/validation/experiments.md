---
title: Wind-tunnel data
description: Scirocco and XFOIL against published NACA 0012 and NACA 4412 measurements, lift, drag and surface pressure.
---

# Comparison with wind-tunnel data

Scirocco and XFOIL agree closely with each other (see
[XFOIL and mfoil](./xfoil-mfoil)), so against measurements they share the
same strengths and errors:

| Quantity | Case | Scirocco | XFOIL | Measured |
| --- | --- | --- | --- | --- |
| Lift slope, α = 0…8° | NACA 0012, Re = 6×10⁶ | 0.112 /deg | 0.113 /deg | 0.108 /deg (Ladson) |
| Cl,max | same | 1.74 at 18° | 1.74 at 18° | 1.60 at 16.4° (Abbott), 1.61 at 17.1° (Ladson) |
| Cd at Cl = 0…1.25 | same | 9–17% below, 14% on average | 9–17% below, 13% on average | Abbott & von Doenhoff |
| Cd at Cl = 0…1.25, tripped at 5% chord | same | 2–4% below on average | 2–4% below on average | Ladson, three grit sizes |
| Upper-surface Cp, RMS error | NACA 0012, α = 0 / 10 / 15° | 0.018 / 0.047 / 0.25 | 0.016 / 0.043 / 0.26 | Gregory & O'Reilly |
| Cp, RMS error | NACA 4412, α = 13.87° | 0.19 | 0.19 | Coles & Wadcock |

- **Attached flow is predicted well.** The lift slope is 3.5% above
  Ladson's, and the pressure at α = 0° is within 0.02 RMS of the data.
- **Stall is not.** Both codes keep gaining lift for 1–2° past the measured
  stall and over-predict Cl,max by about 8%. The NACA 0012 stalls abruptly
  from the leading edge, which an integral boundary layer coupled to a
  panel method doesn't capture.
- **Tripped drag is close.** With the trips of Ladson's test, both codes
  are within 2–4% of his drag on average.
- **Free-transition drag is low.** Both codes are 9–17% below Abbott and
  von Doenhoff's free-transition drag. The codes agree with each other to
  1%, so the gap is in the shared model or in the conditions, not in
  Scirocco's implementation. Possible contributors are ncrit (the tunnel's turbulence
  level isn't known well enough to set it) and the accuracy of the
  digitized data.

## NACA 0012 at Re = 6×10⁶

<ExperimentPolarPlot />

M = 0.15, ncrit = 9, 199 panels, free transition. Scirocco didn't converge
at α = 3°, between two converged points.

- **Lift:** Ladson's data are tripped with 80-grit roughness strips; the
  trip barely changes lift. Abbott and von Doenhoff's are untripped.
- **Drag:** only Abbott and von Doenhoff's data are shown. Ladson's tripped
  drag is 0.0081 at zero lift, 37% above Abbott and von Doenhoff's
  untripped value; the next section compares it with tripped solutions.

### Tripped drag

Ladson fixed transition with roughness strips at 5% chord on both surfaces.
Setting [`xtrip_upper` and `xtrip_lower`](../concepts/transition#forced-transition)
to 0.05 (`XTR 0.05 0.05` in XFOIL) reproduces that:

<TrippedPolarPlot />

- **Zero lift:** Cd is 0.00791 (Scirocco) and 0.00793 (XFOIL) against
  0.0080–0.0081 measured. With free transition the codes give 0.0051.
- **Up to Cl = 1.25:** Scirocco is 4.0%, 2.7% and 1.9% below the 80-, 120-
  and 180-grit data on average (XFOIL 3.7%, 2.4% and 1.6%), against 25%
  below with free transition. Point by point the difference runs from 8.5%
  below to 2.9% above.
- **Code to code:** Scirocco's Cd is within 1.0% of XFOIL's and its Cl
  within 0.4%. Above α = 5° the upper surface transitions freely ahead of
  the trip in both codes, and from 15° the stagnation point has passed the
  lower trip, which both codes then ignore.
- **Near the stagnation point:** at α = 13° and 14° the lower trip sits
  just behind the stagnation point. Scirocco's Cd is 0.9% and 1.3% below
  XFOIL's there. These points rely on the
  [laminar skin-friction floor](../concepts/transition#trips-near-the-leading-edge).

## NACA 0012 surface pressure

Gregory and O'Reilly measured the upper-surface pressure at Re = 2.88×10⁶.
Their boundary layer was tripped, the codes' isn't; the pressure barely
depends on it.

<GregoryCpPlot alpha={0} />

<GregoryCpPlot alpha={10} />

<GregoryCpPlot alpha={15} />

At 0° the codes are within the scatter of the data. At 10° and 15° they
over-predict the suction peak, by 10% at 10° (−Cp of 5.76 against 5.24)
and 8% at 15°. Part of that is the measurement: the peak is narrow, and
pressure taps can miss its top. The data were also digitized from a
photocopy. At 15° the codes again run ahead of the measured lift, which
stalls at about 16°.

## NACA 4412 near maximum lift

Coles and Wadcock measured the pressure on both surfaces at α = 13.87° and
Re = 1.52×10⁶, with the flow separated over the last part of the upper
surface. This is a hard case for any method:

<ColesWadcockCpPlot />

- **Leading edge:** the codes over-predict the suction peak, −Cp of 6.97
  (Scirocco) and 7.05 (XFOIL) against 6.21.
- **Separated region:** the measured upper-surface pressure flattens at
  Cp ≈ −0.40 from about x/c = 0.8 to the trailing edge. The codes separate
  at 0.82c (Scirocco) and 0.84c (XFOIL), but their pressure recovers to
  about −0.07: they predict a smaller separated region.
- **Lower surface:** the codes' Cp is 0.05–0.18 higher than measured.

The experiment tripped both surfaces (at 2.5% chord on the upper surface,
10.3% on the lower) and ran in a small tunnel. The
[NASA Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/naca4412sep_val.html)
considers it a weak validation case for that reason.

## Sources

The measurements are from the
[NASA Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/naca0012_val.html),
stored with their headers in
[`scirocco/validation/vv/experiments`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/vv/experiments):

- Abbott, I. H. and von Doenhoff, A. E., *Theory of Wing Sections*, Dover,
  1959. NACA 0012 at Re = 6×10⁶, untripped. Digitized by the TMR, which
  notes the digitizing is approximate.
- Ladson, C. L., "Effects of Independent Variation of Mach and Reynolds
  Numbers on the Low-Speed Aerodynamic Characteristics of the NACA 0012
  Airfoil Section," NASA TM-4074, 1988. Re = 6×10⁶, M = 0.15, 80-grit trip.
- Gregory, N. and O'Reilly, C. L., "Low-Speed Aerodynamic Characteristics
  of NACA 0012 Aerofoil Sections, including the Effects of Upper-Surface
  Roughness Simulation Hoar Frost," ARC R&M 3726, 1970.
- Coles, D. and Wadcock, A. J., "Flying-Hot-Wire Study of Flow Past an NACA
  4412 Airfoil at Maximum Lift," *AIAA Journal* 17(4), 1979, pp. 321–329.

The Eppler 387 comparison against the UIUC low-speed measurements is on
[Transition validation](./transition#captured-transition-against-measured-drag).

## Reproducing

`vv.py experiments` in
[`scirocco/validation/vv`](https://github.com/chrislupp/Scirocco/tree/develop/scirocco/validation/vv)
runs both codes and writes `docs/src/data/vv_experiments.json`. The
airfoils are the standard NACA sections with their blunt trailing edges,
paneled by mfoil with 199 panels, as on the code comparison page.
