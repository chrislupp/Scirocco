# Captured transition prototype (mfoil)

`mfoil/mfoil_captured.py` prototypes a smooth, "captured" laminar–turbulent
transition model, as an alternative to the cut cell that XFOIL and mfoil use.
Its module docstring documents the formulation and the parameters. This
directory holds the evidence and the scripts that produce it. The C++
library implements the same model, with its adjoint, as
`TransitionModel::Captured`.

## The problem

With the cut cell, the discrete equations switch from laminar to turbulent
at a transition point inside one interval. When that point crosses a node,
the station changes equation type. Outputs stay continuous, but their
derivatives jump. As the transition point moves from panel to panel, the
exact (complex-step) dCd/dα saw-tooths by ±20–50%.

## The model

Both the amplification factor ñ and the shear stress c_τ are carried at every
node, so the state is `[θ, δ*, c_τ, u_e, ñ]`. Every interval uses the same four
equations. Smooth functions of ñ do three things:

- release c_τ into the lag equation just before onset;
- switch the growth of ñ at onset (ñ = ncrit, the cut cell's transition
  point) to at least C/θ;
- blend the closures H*, c_f and c_D from laminar to turbulent over a
  transition length of L momentum thicknesses.

No station ever changes equation type.

## Results

The validation results are on the documentation page
[Transition validation](../../../docs/docs/concepts/transition-validation.md):

- drag against the UIUC Eppler 387 measurements, which sets `capt_L`;
- dCd/dα smoothness against the cut cell;
- the drag test matrix;
- the complex-step derivative checks.

## Reproducing

- `alpha_sweep.py MODE PARAMS A0 A1 N [NPANEL]` writes a complex-step
  dCd/dα sweep (`data/sweep_*.txt`).
- `test_matrix.py` runs the 27-case drag matrix (`data/matrix.txt`).
- `e387_calibrate.py fetch WORKDIR` downloads the UIUC data and writes the
  gapped coordinates. Build `e387_polar.cpp` against the library and run it
  once per Re and L. The polars go in `data/e387/`, with prefix `pol` for
  ncrit 9 and `pol8` for ncrit 8.
- `e387_calibrate.py compare WORKDIR PREFIX L...` prints the error table.
- `e387_calibrate.py docs WORKDIR ../../../docs/src/data/transition.json`
  writes the documentation figures.

The UIUC data are not stored here: their terms require their notices to
accompany any copy.

## Transition length

`transition_rate` implements an optional physical length, selected with
`capt_Lmodel = 'narasimha'`:

- Dhawan–Narasimha intermittency is γ = 1 − exp(−0.41ξ²), with
  ξ = (x − x_t)/λ. The 10–90% length is 1.863λ (Dey et al. 1995, J. Indian
  Inst. Sci. 75, eq. 3).
- λ/θ_t = √(0.41 Re_θt / N), from N = 0.41 Re_θt³ / Re_λ² (same source,
  eq. 10).
- N₀ = 0.7×10⁻³ for turbulence-driven transition. Narasimha (1990, NASA
  CR-187487, §5) notes that N is higher, and facility dependent, when
  residual disturbances drive transition.
- The pressure-gradient increase in N comes from Narasimha (1990, Fig. 3),
  after Dey & Narasimha (1989). For favorable gradients it is 0.25 L², as
  printed on the figure. For adverse gradients, 5.25 L² + 25.2 L³ is my fit
  to the digitized curve through Gostelow's (1989) data, for
  −0.07 ≤ L ≤ 0. Here L = θ²/ν dU/dx is evaluated locally. The report calls
  this figure "a preliminary proposal".

The resulting lengths are 3–6 times longer than 80θ (the default before
calibration), and they lower Cd by 3–11% at Re ≥ 10⁶. The likely reason is double counting. The shear-lag
equation starts c_τ at XFOIL's c_τ,tr ≈ 0.2 c_τ,eq, so the lag already models
much of the transition ramp, and the full intermittency length stacks a
second ramp on top. A consistent physical model would need to recalibrate
c_τ,tr (XFOIL's CtauC and CtauE) together with L. The default is therefore a
constant L, calibrated against measured drag (see Results).

## What did not work

These dead ends shaped the formulation:

1. **Blending whole station residuals.** The shape equation contains
   `ln(H*₂/H*₁)`, and a blend of laminar and turbulent log-differences does not
   telescope. The leftover source term depends on node position, and the
   result was rougher than the cut cell. Blending the closures per station
   fixes this.
2. **Letting the laminar amplification rate carry ñ through the blend.** The
   blend lowers H_k, so that rate collapses and ñ stalls. This produced
   multiple converged solutions: at two adjacent α, transition landed 0.05c
   apart. The fix is the profile-independent floor C/θ.
3. **A hard floor C/θ alone.** It is too slow inside laminar separation
   bubbles, and over-predicted their drag by up to 22%. A smooth
   max(laminar rate, C/θ) fixes this.
4. **Releasing c_τ at onset.** This coincides with the fast growth of ñ, and
   Newton fails. c_τ is released 2.5 units of ñ earlier.
5. **Extrapolating c_τ to the stagnation point.** It can go negative because
   c_τ,eq changes quickly at Re_θ ≈ 10. Nine Re = 3×10⁶ cases failed until
   the first node's value was used instead.
6. **A short transition length (30θ).** The zone then falls inside about one
   panel, which brings back a cell-scale event.

## Open items

- L is calibrated on a single airfoil (E387). Measured transition
  locations, such as the E387 oil-flow data in NASA TM-4062, would tighten
  it. A joint fit of c_τ,tr and the Narasimha length is the more physical
  alternative.
- Only prescribed α is supported (no lift-target mode).
