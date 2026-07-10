---
sidebar_position: 5
---

# 5. Check analytic sensitivities

`SiroccoDiscipline` exposes analytic partial derivatives of `cl`, `cd`, and
`cm` with respect to coordinates, angle of attack, Reynolds number, and Mach
number. The Philote-MDO Python binding maps those derivatives into OpenMDAO,
so the Python side uses the same differentiation tools as a local component.

Continue with either the two-server CST model from tutorial 3 or the local
Python geometry model from tutorial 4.

## Check the remote component

Run the model once at a smooth, converged operating point, then compare its
analytic Jacobian with finite differences:

```python
problem.run_model()

problem.check_partials(
    includes=["analysis"],
    method="fd",
    step=1.0e-6,
    compact_print=True,
)
```

The report covers `(cl, cd, cm)` with respect to `(x, alpha, Re, Ma)`.
Coordinate columns are flattened in the same order as `x`: every node
x-coordinate, followed by every node z-coordinate. Alpha derivatives are per
degree.

Finite differences are only the check here. Optimization drivers receive the
analytic values returned by `SiroccoDiscipline::ComputePartials`.

## Check the complete CST-to-aerodynamics chain

For the two-server model, declare outputs of interest and the inputs you may
later turn into design variables:

```python
problem.model.add_design_var("cst_upper")
problem.model.add_design_var("cst_lower")
problem.model.add_design_var("analysis.alpha")
problem.model.add_objective("analysis.cd")
```

Make these declarations before `problem.setup()`. After setup and a successful
`run_model()`, check total derivatives:

```python
problem.check_totals(
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

For a smooth converged case, analytic and finite-difference values should
closely agree. Investigate large errors in this order:

1. Confirm the primal solve returns finite, repeatable outputs.
2. Move away from stall or a point where transition jumps between panels.
3. Reduce the finite-difference step only after the first two checks.

Reynolds-number derivatives can be numerically awkward because `Re` is large;
use a step strategy appropriate to its scale when checking that column. The
repository validates the native adjoint against complex-step reference data;
this finite-difference check verifies the remote Python-to-Philote wiring.

## Next: [6. Record an HDF5 history](./hdf5-output)
