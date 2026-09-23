---
sidebar_position: 5
---

# 5. Check analytic sensitivities

`SiroccoDiscipline` exposes analytic partial derivatives of `cl`, `cd`, and
`cm` with respect to coordinates, angle of attack, Reynolds number, and Mach
number. The Philote-MDO Python binding maps those derivatives into OpenMDAO,
so the Python side uses the same differentiation tools as a local component.

This tutorial checks the complete two-server CST model because that is the
derivative chain used by the optimization examples. Run the tested file:

```bash
python examples/tutorials/05_sensitivities.py
```

## Declare the derivative problem

Declare the design variables and an objective before `problem.setup()`:

```python
problem.model.add_design_var("cst_upper")
problem.model.add_design_var("cst_lower")
problem.model.add_design_var("analysis.alpha")
problem.model.add_objective("analysis.cd")
```

Use the same inviscid operating point as the quickstart. After setup and a
successful `run_model()`, compare analytic totals with finite differences:

```python
totals = problem.check_totals(
    of=["analysis.cl", "analysis.cd", "analysis.cm"],
    wrt=["cst_upper", "cst_lower", "analysis.alpha"],
    method="fd",
    step=1.0e-6,
    compact_print=True,
)
```

OpenMDAO assembles the chain rule across both remote components:

$$
\frac{d C_d}{d A}
= \frac{\partial C_d}{\partial x}
  \frac{\partial x}{\partial A},
$$

where $A$ is a CST coefficient. The first factor comes from Scirocco's
adjoint; the second comes from `KulfanDiscipline`'s analytic Bernstein-basis
Jacobian.

## Reading the report

On the tested environment, the largest finite relative error in a nonzero row
is about `3e-6`. The inviscid `Cd` rows have zero analytic magnitude, so judge
them by absolute error; the maximum is about `6e-5` with this forward-difference
step.

Finite differences are only the check. Optimization drivers receive the
analytic values returned by `SiroccoDiscipline::ComputePartials`. Investigate
larger errors in this order:

1. Confirm the primal solve returns finite, repeatable outputs.
2. Move away from stall or a point where transition jumps between panels.
3. Reduce the finite-difference step only after the first two checks.

Do not run an absolute forward-difference check of `Re` at `Re=0`: any positive
perturbation switches the model from the inviscid to the viscous solve path, so
that point is intentionally not differentiable with respect to `Re`. Check a
Reynolds-number derivative at a positive baseline and scale the finite-
difference step to the magnitude of `Re`.

The repository validates the native adjoint against complex-step reference
data; this tutorial verifies the remote Python-to-Philote chain.

## Next: [6. Record an HDF5 history](./hdf5-output)
