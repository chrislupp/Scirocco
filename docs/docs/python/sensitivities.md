---
title: 5. Check the gradients
description: Compare the analytic CST-to-forces gradient against finite differences.
---

# 5. Check the gradients

Both servers return analytic partials, and OpenMDAO chains them into total
derivatives. This tutorial checks the whole chain, from CST coefficients
to forces, against finite differences. It's the same chain the optimization
examples use.

```bash
python examples/tutorials/05_sensitivities.py
```

## Declare the derivatives

Before `problem.setup()`, declare design variables and an objective:

```python
problem.model.add_design_var("cst_upper")
problem.model.add_design_var("cst_lower")
problem.model.add_design_var("analysis.alpha")
problem.model.add_objective("analysis.cd")
```

After a successful `run_model()` at the quickstart's operating point,
compare the analytic totals with finite differences:

```python
totals = problem.check_totals(
    of=["analysis.cl", "analysis.cd", "analysis.cm"],
    wrt=["cst_upper", "cst_lower", "analysis.alpha"],
    method="fd",
    step=1.0e-6,
    compact_print=True,
)
```

For a CST coefficient $A$, OpenMDAO computes

$$
\frac{dC_d}{dA} = \frac{\partial C_d}{\partial x}\,\frac{\partial x}{\partial A},
$$

with the first factor from Scirocco's adjoint and the second from the CST
server's exact Bernstein-basis Jacobian.

## Reading the report

A reference run gives a largest absolute error of about `5e-5` and a
largest relative error of about `4.4e-6`. The script fails if either
exceeds `1e-4`. This run is inviscid, so the `Cd` rows are zero; judge them
by absolute error.

The finite differences are only a check. Optimizers get the analytic values.
If the errors are larger:

1. Make sure the primal solve converges and repeats exactly.
2. Move away from stall, and away from points where transition jumps between
   panels.
3. Only then try a different finite-difference step.

Don't finite-difference `Re` from `Re = 0`: any positive step switches to
the viscous model. Check `Re` derivatives from a positive baseline, with a
step scaled to `Re`.

Inside the repository, the C++ adjoint is tested against complex-step
reference data, which is exact. This check covers the Python-to-server chain
on top of that.

<Next items={[
  { to: "/docs/python/hdf5-history", title: "6. Record a history", body: "Write every evaluation to HDF5." },
]} />
