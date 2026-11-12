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

## Drag against the cut cell

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

The cut cell is not ground truth here. It is XFOIL's calibrated model, and the
low-Re differences sit in separation bubbles, where neither model has been
checked against experiment in this study.

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
c_τ,tr (XFOIL's CtauC and CtauE) together with L. That is why the default is
a constant L = 80θ, chosen to stay close to the validated cut-cell drag at
Re ≥ 10⁶.

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

- Calibrate L against experiment, ideally together with c_τ,tr, or adopt the
  Narasimha length with a recalibrated c_τ,tr.
- Only prescribed α is supported (no lift-target mode).
- C++ port: done as `TransitionModel::Captured`, with a fully analytic adjoint
  (`CapturedTransition*` tests). Shape sensitivities match mfoil complex step
  to 1e-8, and alpha, Re and Mach to 2e-11. The Re and Mach residual partials
  chain the closures' `*_dRet` / `*_dMa` derivatives through the captured
  residual, as `ComputeStationLinRet` and `ComputeStationLinMa` do for the cut
  cell.
