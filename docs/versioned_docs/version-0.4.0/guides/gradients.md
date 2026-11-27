---
title: Compute gradients
description: Adjoint derivatives of Cl, Cd, and Cm with respect to coordinates, angle of attack, Reynolds number, and Mach number.
---

# Compute gradients

After a successful solve, `ComputeSensitivities()` returns the derivatives of
$C_l$, $C_d$, and $C_m$ with respect to every panel-node coordinate and the
three operating conditions. One call gives you all of them.

```cpp
scirocco::Scirocco foil("2412", 160);
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
if (!foil.Solve()) return 1;

const scirocco::Sensitivity s = foil.ComputeSensitivities();
std::printf("dCl/dalpha = %.6f per deg\n", s.dcl_dalpha);
std::printf("dCd/dalpha = %.6f per deg\n", s.dcd_dalpha);
std::printf("dCd/dRe    = %.4e\n", s.dcd_dRe);
std::printf("dCd/dx is %ld x %ld\n", (long)s.dcd_dx.rows(), (long)s.dcd_dx.cols());
```

```output
dCl/dalpha = 0.092430 per deg
dCd/dalpha = 0.000748 per deg
dCd/dRe    = -1.3061e-09
dCd/dx is 2 x 161
```

A central finite difference of $C_d$ in $\alpha$, with a step of $10^{-3}$°,
gives `0.000748` as well.

## What you get

| Field | Shape | Derivative of |
| --- | --- | --- |
| `dcl_dx`, `dcd_dx`, `dcm_dx` | `2 × N` | Each output with respect to each node's x (row 0) and z (row 1) |
| `dcl_dalpha`, `dcd_dalpha`, `dcm_dalpha` | scalar | Per degree of angle of attack |
| `dcl_dRe`, `dcd_dRe`, `dcm_dRe` | scalar | Per unit Reynolds number |
| `dcl_dMa`, `dcd_dMa`, `dcm_dMa` | scalar | Per unit Mach number |

Node order matches `GetFoilCoords()`. The same values are available from
getters such as `foil.GetDcdDx()` and `foil.GetDclDalpha()` until the next
solve.

## Chain them to your design variables

The coordinate gradients are with respect to **panel nodes**. If your design
variables are something else, say CST coefficients or a thickness scale,
multiply by the Jacobian of node coordinates with respect to those variables:

$$
\frac{dC_d}{dA} = \frac{\partial C_d}{\partial x}\,\frac{\partial x}{\partial A}
$$

An MDO framework does this for you when the geometry is its own discipline;
[Check the gradients](../python/sensitivities) shows it with OpenMDAO.

## Lift-target mode

With a lift target (`cl ≥ 0` in `SetOper`), the derivatives include the change
in angle of attack needed to keep $C_l$ on target. `dcl_dx` is then zero and
extra fields describe the trim:

| Field | Meaning |
| --- | --- |
| `dcd_dcltarget`, `dcm_dcltarget` | Change per unit change of the target |
| `dalpha_dx`, `dalpha_dRe`, `dalpha_dMa` | How the solved angle moves |

## Cost

The adjoint solves one linear system per output, so the cost does not grow
with the number of design variables. A finite-difference gradient needs one
extra solve per input: `2N + 3` of them.

On an Apple M1 Pro, the full gradient takes about 9 ms at 99 panels and
18 ms at 149 panels, less than the viscous solve before it. Finite
differences over the same inputs would take about 2.5 s and 7 s. See
[Performance](../validation/performance) for the full table.

## Where the shape gradient is rough

On the same case, the drag gradient along the upper surface is smooth except
for a few panels near the leading edge, where it jumps and alternates in sign:

<ShapeGradientPlot />

This is the exact derivative of the discrete model, not an adjoint error.
Complex step through mfoil gives the same values at every node, and the
test suite checks this case. The spike sits where the
[amplification factor](../concepts/transition#limits) starts to grow. The
model switches growth on over a narrow band of $Re_\theta$, so near that
point the growth rate is very sensitive to the local boundary-layer state.
Moving one node there changes how much the disturbance grows downstream,
and so where transition happens at $x/c \approx 0.39$.

The spike comes from how sharply the model switches growth on, not from the
panel count: refining the paneling moves it but does not smooth it away.
Other airfoils and operating points put the onset elsewhere. For
gradient-based shape optimization:

- Parameterize the shape with smooth design variables, such as CST
  coefficients. A smooth parameterization averages the gradient over many
  nodes, so single-node spikes don't dominate it.
- Don't read node-by-node gradients near the onset as physical design
  guidance.

## Rules

- **Solve first.** Gradients are taken at the converged state. A failed
  solve, or a call before any solve, throws.
- **Don't touch the state in between.** A new solve, a new operating point,
  a geometry change, or the non-const `GetState()` clears cached gradients.
- **Positions, not shape.** Moving a node doesn't change the reference
  chord or the moment reference point; they are held fixed.

The derivatives are exact for the discretized equations, to the solver's
tolerance. Every piece of the chain is tested against complex-step reference
data from mfoil. [The adjoint](../concepts/adjoint) explains the method.

<Next items={[
  { to: "/docs/api/sensitivity", title: "Sensitivity reference", body: "Every field, arc-length modes, and near-field drag." },
  { to: "/docs/concepts/adjoint", title: "The adjoint", body: "How one linear solve gives the whole gradient." },
]} />
