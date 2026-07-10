---
sidebar_position: 3
---

# Python quickstart

This is the shortest path from a build to an aerodynamic result. Scirocco runs
in two Philote-Cpp gRPC servers and Python uses the Philote-MDO bindings to
present them as OpenMDAO components:

```text
Python/OpenMDAO -> KulfanDiscipline -> panel coordinates -> SiroccoDiscipline
```

Scirocco does not currently install an in-process `import scirocco` module.
In these tutorials, **Python bindings** means
`philote_mdo.openmdao.RemoteExplicitComponent`, the Python binding for a
remote Philote discipline.

## Start the servers

Build the Philote targets as described in [Installation](./installation), then
run these in two terminals from the repository root:

```bash
scirocco/build-philote/src/philote/scirocco_kulfan_server localhost:50052
```

```bash
scirocco/build-philote/src/philote/scirocco_server localhost:50051
```

The first server creates a Kulfan/CST airfoil. The second analyzes the panel
coordinates it receives.

## Run one analysis

Save this as `quickstart.py`:

```python
import grpc
import numpy as np
import openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent

npanel = 100

problem = om.Problem()
problem.model.add_subsystem(
    "geometry",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50052"),
        n_cst_upper=6,
        n_cst_lower=6,
        npanel=npanel,
    ),
    promotes_inputs=["cst_upper", "cst_lower"],
    promotes_outputs=["x"],
)
problem.model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50051"),
        naca="0012",       # initializes the remote solver
        npanel=npanel,
    ),
    promotes_inputs=["x"],
)

problem.setup()

# Symmetric, NACA-0012-like CST coefficients. Positive coefficients produce
# positive thickness; use the same values on the upper and lower surfaces.
cst_0012 = np.array([0.170, 0.160, 0.145, 0.145, 0.160, 0.170])
problem.set_val("cst_upper", cst_0012)
problem.set_val("cst_lower", cst_0012)
problem.set_val("analysis.alpha", 4.0, units="deg")
problem.set_val("analysis.Re", 0.0)  # zero selects inviscid analysis
problem.set_val("analysis.Ma", 0.0)

problem.run_model()

print("Cl =", problem.get_val("analysis.cl")[0])
print("Cd =", problem.get_val("analysis.cd")[0])
print("Cm =", problem.get_val("analysis.cm")[0])
```

Run it:

```bash
python quickstart.py
```

You should get positive lift and essentially zero drag. Inviscid 2-D potential
flow has no profile drag; the [first tutorial](../tutorials/inviscid-naca)
turns this one-off model into an angle-of-attack sweep, and the tutorial after
that enables the viscous solver.

## Why both disciplines?

`SiroccoDiscipline` takes panel coordinates as its `x` input. Its `naca`
option establishes the solver and expected panel count, but an upstream
component still supplies the coordinates on every evaluation. Starting with
`KulfanDiscipline` makes that data flow explicit and gives us a parameterized
geometry with analytic derivatives that later tutorials can optimize.
