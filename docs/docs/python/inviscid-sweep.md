---
title: 1. Sweep angle of attack
description: Reuse one model across operating points.
---

# 1. Sweep angle of attack

Start from the [overview model](./index.md#run-one-analysis) and run it at nine
angles of attack. Keep both servers running.

```bash
python examples/tutorials/01_inviscid_sweep.py
```

Everything up to `problem.setup()` stays the same. Replace the rest with a
loop:

```python
problem.set_val("cst_upper", cst_0012)
problem.set_val("cst_lower", cst_0012)
problem.set_val("analysis.Re", 0.0)
problem.set_val("analysis.Ma", 0.0)

print(" alpha       Cl          Cd          Cm")
for alpha in np.arange(-4.0, 12.1, 2.0):
    problem.set_val("analysis.alpha", alpha, units="deg")
    problem.run_model()
    cl = problem.get_val("analysis.cl")[0]
    cd = problem.get_val("analysis.cd")[0]
    cm = problem.get_val("analysis.cm")[0]
    print(f"{alpha:6.1f}  {cl:9.5f}  {cd:10.3e}  {cm:9.5f}")
```

```output
alpha       Cl          Cd          Cm
  -4.0   -0.48270   0.000e+00    0.00569
  -2.0   -0.24150   0.000e+00    0.00285
   0.0   -0.00000   0.000e+00    0.00000
   2.0    0.24150   0.000e+00   -0.00285
   4.0    0.48270   0.000e+00   -0.00569
   6.0    0.72333   0.000e+00   -0.00850
   8.0    0.96310   0.000e+00   -0.01127
  10.0    1.20171   0.000e+00   -0.01399
  12.0    1.43889   0.000e+00   -0.01663
```

Each `run_model()` sends the current inputs to both servers: the geometry
server rebuilds `x`, and the analysis server solves at the new angle.
Reusing the same `Problem` skips the setup cost.

The section is symmetric, so lift is antisymmetric in $\alpha$ and zero at
0°. Lift grows almost linearly, and drag is zero because 2D potential flow
has no drag. Turning on viscosity changes both, in the next tutorial.

## What `x` contains

The two components connect through `x`: all node x-coordinates, then all
node z-coordinates, from the lower trailing edge around to the upper one.
Both components must use the same `npanel`, giving `npanel + 1` nodes.

<Next items={[
  { to: "/docs/python/viscous-analysis", title: "2. Turn on viscosity", body: "Reynolds and Mach number sweeps." },
]} />
