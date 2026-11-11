# Captured transition prototype (mfoil)

`mfoil/mfoil_captured.py` prototypes a smooth, "captured" laminar–turbulent
transition model, as an alternative to the cut cell that XFOIL and mfoil use.
Its module docstring documents the formulation and the parameters. This
directory holds the evidence and the sweep script. The model is not in the
C++ library yet.

## The problem

With the cut cell, the discrete equations switch from laminar to turbulent
at a transition point inside one interval. When that point crosses a node,
the station changes equation type. Outputs stay continuous, but their
derivatives jump. On NACA 2412 at Re = 10⁶ with 160 panels, the exact
(complex-step) dCd/dα saw-tooths between 4.0e-4 and 9.7e-4 for α = 3.9–4.3°.
It jumps five times, each time the upper-surface transition moves one panel.

## The model

Both the amplification factor ñ and the shear stress c_τ are carried at every
node, so the state is `[θ, δ*, c_τ, u_e, ñ]`. Every interval uses the same four
equations. Smooth functions of ñ blend the closures, release c_τ into the lag
equation, and switch the growth of ñ at onset. No station ever changes
equation type.

## Results

These are exact complex-step derivatives, from `alpha_sweep.py`, with data in
`data/`:

| model | dCd/dα, α = 3.9 → 4.3° | ripple | Cd at α = 4° |
| --- | --- | --- | --- |
| cut cell (mfoil default) | 4.0e-4 ↔ 9.7e-4, 5 jumps | ±40% | 0.0069320 |
| captured, b = 0, L = 80θ (default) | 6.4e-4 → 7.7e-4 | ±3% | 0.0069045 (−0.4%) |
| captured, b = 0, L = 120θ | 6.2e-4 → 7.4e-4 | ±1.5% | 0.0067671 (−2.4%) |
| captured, b = 2, L = 120θ | 5.7e-4 → 6.8e-4 | ±0.3% | 0.0065392 (−5.7%) |

- Every captured point converges quadratically, in 7–8 Newton iterations
  from the converted cut-cell solution.
- The captured Jacobian matches complex step to 1e-15 (excluding the
  stagnation-point chain, which only `stagpoint_move` applies).
- The mean slopes agree: 6.6e-4 for the cut cell and 7.0e-4 captured.

## What did not work

These dead ends shaped the formulation:

1. **Blending whole station residuals.** This converges, but the result is
   rougher than the cut cell. The shape equation contains `ln(H*₂/H*₁)`, and a
   blend of laminar and turbulent log-differences does not telescope. The
   leftover source term depends on where the blend sits relative to the
   nodes. Blending the closures per station fixes this.
2. **Letting the laminar amplification rate drive ñ through the blend.**
   Turbulent closures lower H_k, the laminar rate collapses, and ñ stalls.
   The result is multiple converged solutions: at two adjacent α, transition
   landed 0.05c apart. After onset, ñ now grows at the profile-independent
   rate C/θ.
3. **Releasing c_τ at onset.** This coincides with the fast growth of ñ, and
   Newton fails. c_τ is released 2.5 units of ñ earlier. Until then it tracks
   c_τ,tr, which is only meaningful for laminar profiles.
4. **A short transition length (30θ).** This puts the blend inside one panel,
   which brings back a cell-scale event. Transition lengths from
   Dhawan–Narasimha are much longer: λ/θ_t = √(0.412 Re_θt / N), several
   hundred θ in zero pressure gradient and about 70–200θ under airfoil
   adverse gradients. A physically sized zone spans several panels.

## Open items

- L is a constant. It should come from a spot-formation-rate correlation
  (pressure gradient and Tu, for example Gostelow, Blunden and Walker),
  evaluated locally.
- The ±3% ripple at b = 0 is panel-scale residue from the onset switch.
- Only prescribed α is supported (no lift-target mode), and only NACA 2412 at
  Re = 10⁶ has been tested.
- The C++ port: a 5-variable boundary-layer state touches the residuals, the
  Jacobians, the adjoint chains (geometry, α, Re, Mach), the wake, the
  stagnation handling, and the HDF5 reference data.
