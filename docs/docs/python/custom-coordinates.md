---
title: 4. Your own coordinates
description: Feed coordinates from Python straight into the analysis server.
---

# 4. Your own coordinates

Any component that outputs `x` in the right layout can replace the geometry
server. This tutorial builds NACA 0012 coordinates in NumPy. Only the
analysis server (port 50051) needs to be running.

```bash
python examples/tutorials/04_custom_airfoil.py
```

## Build the coordinates

The analysis server expects `npanel + 1` nodes from the lower trailing edge
around to the upper one, flattened as all x-values, then all z-values:

```python
import numpy as np

def naca0012_coordinates(npanel=100, thickness=0.12):
    if npanel % 2:
        raise ValueError("use an even panel count for this constructor")

    nhalf = npanel // 2
    beta = np.linspace(0.0, np.pi, nhalf + 1)
    x_lower = 0.5 * (1.0 + np.cos(beta))      # 1 -> 0
    x_upper = x_lower[-2::-1]                 # skip the repeated LE; 0 -> 1

    def half_thickness(x):
        return (thickness / 0.2) * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            - 0.1015 * x**4
        )

    z_lower = -half_thickness(x_lower)
    z_upper = +half_thickness(x_upper)
    return np.concatenate([x_lower, x_upper, z_lower, z_upper])
```

Cosine spacing clusters nodes at both edges. The server uses your nodes as
the panels exactly; it doesn't re-panel them.

## Connect it

```python
import grpc
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent

npanel = 100
coordinates = naca0012_coordinates(npanel=npanel)

problem = om.Problem()
ivc = problem.model.add_subsystem(
    "inputs", om.IndepVarComp(), promotes_outputs=["x"]
)
ivc.add_output("x", val=coordinates, units="m")

problem.model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50051"),
        naca="0012",
        npanel=npanel,
    ),
    promotes_inputs=["x"],
)

problem.setup()
problem.set_val("analysis.alpha", 4.0, units="deg")
problem.set_val("analysis.Re", 1.0e6)
problem.set_val("analysis.Ma", 0.0)
problem.run_model()

print("Cl =", problem.get_val("analysis.cl")[0])
print("Cd =", problem.get_val("analysis.cd")[0])
print("Cm =", problem.get_val("analysis.cm")[0])
```

```output
Cl = 0.42840697
Cd = 0.00723560
Cm = 0.00581315
```

## Checklist for coordinates

- Exactly `npanel + 1` nodes, lower trailing edge → leading edge → upper
  trailing edge, with the leading edge listed once.
- Flattened as `[x_0 … x_N, z_0 … z_N]`, not as interleaved pairs.
- Unit chord. The variable carries length units, but Scirocco treats values
  as fractions of chord and doesn't rescale them.
- Clustered at the leading and trailing edges.
- A small trailing-edge gap. A viscous solve on a closed trailing edge
  (thinner than 1e-6 chord) raises an error; see [Troubleshooting](../getting-started/troubleshooting#the-solve-doesnt-converge).
- No repeated neighbours, no self-intersections. Plot it with equal axes
  before blaming the solver.

<Next items={[
  { to: "/docs/python/sensitivities", title: "5. Check the gradients", body: "Analytic totals against finite differences." },
]} />
