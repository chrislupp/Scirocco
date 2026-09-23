---
sidebar_position: 4
---

# 4. Supply coordinates from Python

So far the remote Kulfan discipline produced `x`. Any Python or OpenMDAO
component can replace it, provided it emits the same flattened panel-node
layout. This tutorial generates normalized NACA 0012 coordinates locally and
keeps Scirocco behind the Philote analysis server.

Only the analysis server on port 50051 is needed. Run the complete example:

```bash
python examples/tutorials/04_custom_airfoil.py
```

## Build panel coordinates

The analysis discipline expects `npanel + 1` nodes ordered lower trailing edge
→ leading edge → upper trailing edge, with all x-values followed by all
z-values:

```python
def naca0012_coordinates(npanel=100, thickness=0.12):
    if npanel % 2:
        raise ValueError("use an even panel count for this constructor")

    nhalf = npanel // 2
    beta = np.linspace(0.0, np.pi, nhalf + 1)
    x_lower = 0.5 * (1.0 + np.cos(beta))      # 1 -> 0
    x_upper = x_lower[-2::-1]                 # omit duplicate LE; 0 -> 1

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

The result has `2 * (npanel + 1)` entries. A finite trailing-edge gap is also
valid; the upper and lower trailing-edge z-values need not coincide.

## Replace the geometry subsystem

Create a new problem with only the remote analysis discipline:

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

The tested result is:

```text
Cl = 0.42840697
Cd = 0.00723560
Cm = 0.00581315
```

`naca="0012"` is still required to initialize `SiroccoDiscipline`, but the
connected `x` array replaces that seed geometry before every solve.

## Coordinate checklist

- Supply exactly `npanel + 1` nodes.
- Order nodes lower TE → LE → upper TE.
- Flatten as `[x_0 ... x_N, z_0 ... z_N]`, not interleaved pairs.
- Do not repeat the leading-edge node.
- Normalize the chord so the leading edge is near `x=0` and the trailing edge
  is near `x=1`; Scirocco treats these values as nondimensional chord
  coordinates even though the Philote variable carries length units.
- Resample measured or CAD curves to the requested node count with adequate
  leading- and trailing-edge clustering. The server replaces its panel nodes;
  it does not repanel an arbitrary point cloud.
- Reject duplicate consecutive nodes, self-intersections, and a reversed loop
  before sending the array.
- Plot the loop with equal axes before diagnosing aerodynamic surprises.

The next tutorial uses the exact same Python model and asks Philote for the
analytic derivatives of its outputs.

## Next: [5. Check analytic sensitivities](./sensitivities)
