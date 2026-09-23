---
sidebar_position: 3
---

# Python quickstart

This is the shortest path from a build to an aerodynamic result. Scirocco runs
in two Philote-Cpp gRPC servers and Python uses the Philote-MDO bindings to
present them as OpenMDAO components:

```text
Cl = 0.48270391
Cd = 0.000e+00
Cm = -0.00569190
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

The first server creates a Kulfan/CST airfoil with a 0.00126-chord trailing-edge gap. The second analyzes the panel
coordinates it receives.

## Run one analysis

The complete, tested file is
[`examples/tutorials/quickstart.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/tutorials/quickstart.py).
Run it from the repository root:

```bash
python examples/tutorials/quickstart.py
```

Its full model is shown below so the wiring is visible:

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
        te_upper=0.00063,
        te_lower=-0.00063,
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

# Symmetric NACA 0012 CST fit (also used by the multipoint example). Positive coefficients produce
# positive thickness; use the same values on the upper and lower surfaces.
cst_0012 = np.array([
    0.1706635358053721, 0.15626369325824663, 0.1546327069833401,
    0.14462202838023833, 0.13318428664336585, 0.15537328876198608,
])
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

With `philote-mdo==0.8.0` and the Philote-Cpp build pinned in the installation
guide, the output is:

```text
Cl = 0.48319963
Cd = 0.000e+00
Cm = -0.00612678
```

The last digits can vary slightly by compiler and platform. You should still
get `Cl` between 0.47 and 0.50 and drag below `1e-12`. Inviscid 2-D potential
flow has no profile drag; the [first tutorial](../tutorials/inviscid-naca)
turns this one-off model into an angle-of-attack sweep.

## Why both disciplines?

`SiroccoDiscipline` takes panel coordinates as its `x` input. Its `naca`
option establishes the solver and expected panel count, but an upstream
component still supplies the coordinates on every evaluation. Starting with
`KulfanDiscipline` makes that data flow explicit and gives us a parameterized
geometry with analytic derivatives that later tutorials can optimize.

## If the quickstart fails

| Error | Check |
| --- | --- |
| `Option 'n_cst_upper' ... has not been declared` | The server is linked to an older Philote-Cpp. Install the pinned 0.8 implementation, which initializes disciplines during service registration, then rebuild Scirocco. |
| `Exception deserializing response` or invalid UTF-8 | The C++ server is using the Philote 0.7 wire schema. Reinstall the pinned Philote-Cpp 0.8 implementation, then rebuild Scirocco. |
| `UNAVAILABLE` or connection refused | Confirm both server processes are running on ports 50051 and 50052. |
| Results differ substantially from the reference | Confirm `philote-mdo==0.8.0`, `npanel=100`, and `Re=0`. |
