---
title: Python and OpenMDAO
sidebar_label: Overview
description: Run Scirocco from Python as OpenMDAO components, through Philote gRPC servers.
---

# Python and OpenMDAO

From Python, Scirocco runs as OpenMDAO components. Two small servers host
the C++ code, and `philote-mdo` connects to them, so each looks like any
other OpenMDAO component, analytic partials included.

| Server | Port | Inputs | Outputs |
| --- | --- | --- | --- |
| `scirocco_server` | 50051 | `x`, `alpha`, `Re`, `Ma` | `cl`, `cd`, `cm` |
| `scirocco_kulfan_server` | 50052 | `cst_upper`, `cst_lower` | `x`, `smoothness_upper`, `smoothness_lower` |

The geometry server turns CST coefficients into panel coordinates `x`; the
analysis server turns `x` into forces. You can swap the geometry server for
any component that produces `x`.

## Set up

1. **Build the servers.** Follow
   [Installation › Philote servers](../getting-started/installation#philote-servers-for-python),
   or run them in [Docker](./docker).

2. **Install the Python packages** in a virtual environment:

   ```bash
   python3.12 -m venv .venv
   source .venv/bin/activate
   python -m pip install -r examples/tutorials/requirements.txt
   ```

   This pins `philote-mdo==0.8.0`, which must match the server's Philote
   version.

3. **Start both servers**, each in its own terminal:

   ```bash
   scirocco/build-philote/src/philote/scirocco_kulfan_server localhost:50052
   ```

   ```bash
   scirocco/build-philote/src/philote/scirocco_server localhost:50051
   ```

## Run one analysis

```bash
python examples/tutorials/quickstart.py
```

The script builds this model. Here its constants, which live in
`tutorial_airfoil.py`, are written inline:

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

# A CST fit to NACA 0012. Use the same values on both surfaces.
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

print(f"Cl = {problem.get_val('analysis.cl')[0]:.8f}")
print(f"Cd = {problem.get_val('analysis.cd')[0]:.3e}")
print(f"Cm = {problem.get_val('analysis.cm')[0]:.8f}")
```

```output
Cl = 0.48270391
Cd = 0.000e+00
Cm = -0.00569190
```

Expect `Cl` between 0.47 and 0.50; the last digits vary with compiler and
platform. Drag is zero because this run is inviscid. If it fails, see
[Troubleshooting](../getting-started/troubleshooting#connecting-from-python).

The `naca` option only initializes the analysis server and fixes its panel
count. The shape it actually analyzes is whatever arrives on `x`, here
from the CST server.

## Tutorials

Each tutorial changes one thing about this model. The complete scripts are
in [`examples/tutorials`](https://github.com/chrislupp/Scirocco/tree/develop/examples/tutorials)
and run in CI.

<Next items={[
  { to: "/docs/python/inviscid-sweep", title: "1. Sweep angle of attack", body: "Reuse one model across operating points." },
  { to: "/docs/python/viscous-analysis", title: "2. Turn on viscosity", body: "Reynolds and Mach number sweeps." },
  { to: "/docs/python/geometry", title: "3. Change the shape", body: "Thickness and camber through CST coefficients." },
  { to: "/docs/python/custom-coordinates", title: "4. Your own coordinates", body: "Replace the geometry server with Python." },
  { to: "/docs/python/sensitivities", title: "5. Check the gradients", body: "Analytic totals against finite differences." },
  { to: "/docs/python/hdf5-history", title: "6. Record a history", body: "Write every evaluation to HDF5." },
]} />
