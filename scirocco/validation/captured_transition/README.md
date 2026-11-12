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

## Smoothness

`alpha_sweep.py` computes complex-step dCd/dα for α = 3.9–4.3° in 21 points;
the data are in `data/sweep_*.txt`. Ripple is the RMS (and maximum)
deviation of dCd/dα from a quadratic fit, relative to its mean. Every point
converged.

| case | panels | cut cell | captured (L = 80θ) |
| --- | --- | --- | --- |
| NACA 2412, Re = 10⁶ | 160 | 24.6% (max 47.5%) | 1.6% (max 3.4%) |
| NACA 0012, Re = 2×10⁵ (bubble) | 160 | 29.6% (max 51.1%) | 0.2% (max 0.4%) |
| NACA 4412, Re = 3×10⁶ | 160 | 21.2% (max 43.1%) | 4.3% (max 9.7%) |
| NACA 4412, Re = 3×10⁶ | 240 | 20.4% (max 35.9%) | 1.5% (max 2.9%) |
| NACA 4412, Re = 3×10⁶ | 320 | 85% (max 359%)¹ | 1.2% (max 2.1%) |

¹ The sweep script does not check cut-cell convergence, so some of these
points may be unconverged.

- The captured gradient converges under panel refinement; the cut-cell
  gradient does not.
- The remaining ripple comes from resolution: the transition zone has to span
  several panels. At high Re, θ is small, so 80θ covers fewer panels. A
  longer L also helps: at 160 panels, L = 120θ gives 2.4% here.
- These sweeps used L = 80θ, the default before the E387 calibration. With
  the calibrated default L = 60θ, the C++ adjoint gives the same 1.6% ripple
  (max 2.9%) for NACA 2412 at Re = 10⁶.

## Drag against the cut cell (L = 80θ)

`test_matrix.py` runs NACA 0012, 2412 and 4412 at Re = 2×10⁵, 10⁶ and 3×10⁶
and α = 0, 4 and 8°, with 160 panels (`data/matrix.txt`). The captured model
converges in 26 of 27 cases. The failure (2412, 3×10⁶, 8°) also fails with the
cut cell.

The table shows the Cd difference from the cut cell, for cases where the cut
cell converged:

| Re | captured, L = 80θ | captured, Narasimha length |
| --- | --- | --- |
| 3×10⁶ | +0.6% to +2.7% | −2.8% to −9% |
| 10⁶ | −0.5% to +2.9% | −0.8% to −11% |
| 2×10⁵ (laminar separation bubbles) | −1.3% to +8.9% | −2.7% to +9.8% |

The cut cell is not ground truth. The E387 calibration below compares both
models against measurements, and the captured model is the more accurate of
the two.

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

The resulting lengths are 3–6 times longer than 80θ, and they lower Cd by
3–11% at Re ≥ 10⁶. The likely reason is double counting. The shear-lag
equation starts c_τ at XFOIL's c_τ,tr ≈ 0.2 c_τ,eq, so the lag already models
much of the transition ramp, and the full intermittency length stacks a
second ramp on top. A consistent physical model would need to recalibrate
c_τ,tr (XFOIL's CtauC and CtauE) together with L. The default is therefore a
constant L, calibrated against measured drag (next section).

## Calibration against experiment (Eppler 387)

The E387 is the standard low-Reynolds-number transition case: laminar
separation bubbles at every condition below. The reference is the UIUC
drag polars at Re = 1, 2, 3 and 4.6×10⁵ (Williamson et al., *Summary of
Low-Speed Airfoil Data, Vol. 5*, 2012). Cl is within 0.1–1.0, which is
attached flow.

`e387_calibrate.py` downloads the data, and `e387_polar.cpp` computes the
polars. The UIUC data are not stored here. Our polars are in `data/e387/`.
Panel count is 160, and the trailing edge is opened to 0.25% chord because
of issue #18.

Each entry is the Cd error at matched Cl, given as mean signed % / RMS %:

| Re | cut cell | L = 40 | **L = 60** | L = 80 | L = 120 | L = 240 |
| --- | --- | --- | --- | --- | --- | --- |
| 1×10⁵ | −28.4 / 29.0 | −4.9 / 7.1 | −4.0 / 6.3 | −3.9 / 6.2 | −3.9 / 6.2 | −3.8 / 6.2 |
| 2×10⁵ | −10.0 / 10.2 | −3.5 / 3.7 | −2.5 / 3.1 | −2.7 / 3.5 | −2.9 / 4.0 | −3.2 / 4.7 |
| 3×10⁵ | −7.3 / 8.1 | −3.9 / 5.0 | −3.7 / 5.2 | −4.0 / 5.7 | −4.5 / 6.3 | −5.5 / 7.7 |
| 4.6×10⁵ | −6.7 / 7.6 | −4.5 / 5.5 | −4.8 / 6.0 | −5.3 / 6.6 | −6.1 / 7.3 | −7.8 / 9.3 |
| all, ncrit = 9 | −12.8 / 16.0 | −4.2 / 5.4 | **−3.7 / 5.3** | −4.0 / 5.6 | −4.3 / 6.1 | −5.1 / 7.2 |
| all, ncrit = 8 | −14.3 / 17.8 | −7.3 / 9.3 | **−6.7 / 8.5** | −6.9 / 8.6 | −7.5 / 9.1 | |

- **The captured model is about three times more accurate than the cut
  cell.** In laminar separation bubbles the cut cell's instantaneous
  transition reattaches the flow too quickly and under-predicts drag, by 28%
  at Re = 10⁵.
- **L = 60θ is best overall.** It is the best length at ncrit = 9, and within
  0.1% of the best at ncrit = 8. Mack's relation puts ncrit at about 8 for
  the ~0.1% turbulence of the UIUC tunnel. It is the default (`capt_L`). The
  result is flat for 40–80θ, and lengths beyond that over-predict the
  laminar portion at the higher Reynolds numbers.
- **What remains** is a 3–4% under-prediction at all lengths. That points to
  the rest of the closure set (or ncrit), not L.
- **Convergence:** both models converge at 13–16 of the 23 angles (−3 to 8°),
  with failures at the polar ends. The captured model adds no failures of
  its own.
- **A separate limit:** above 160 panels the E387's coarse leading-edge
  definition (61 points) breaks the boundary-layer initialization in both
  Scirocco and mfoil.

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

- L is calibrated on a single airfoil family (E387). More cases with
  measured drag and transition, such as the E387 oil-flow transition data in
  NASA TM-4062, would tighten it. A joint fit of c_τ,tr and the Narasimha
  length is the more physical alternative.
- Only prescribed α is supported (no lift-target mode).
- C++ port: done as `TransitionModel::Captured`, with a fully analytic adjoint
  (`CapturedTransition*` tests). Shape sensitivities match mfoil complex step
  to 1e-8, and alpha, Re and Mach to 2e-11. The Re and Mach residual partials
  chain the closures' `*_dRet` / `*_dMa` derivatives through the captured
  residual, as `ComputeStationLinRet` and `ComputeStationLinMa` do for the cut
  cell.
