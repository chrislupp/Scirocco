---
sidebar_position: 0
---

# Python analysis tutorials

These tutorials form one verified learning path from a single analysis to
geometry derivatives and recorded histories. They use `philote-mdo==0.8.0`
and two local Scirocco gRPC servers.

## Before you start

Complete [Installation](../getting-started/installation), then keep these
commands running in separate terminals:

```bash
scirocco/build-philote/src/philote/scirocco_kulfan_server localhost:50052
```

```bash
scirocco/build-philote/src/philote/scirocco_server localhost:50051
```

Verify the complete quickstart before continuing:

```bash
source .venv/bin/activate
python examples/tutorials/quickstart.py
```

Expect approximately `Cl = 0.48270` and zero drag, using the configuration
and acceptance bounds in the [quickstart](../getting-started/quickstart).
That check catches server, protocol, dependency, and port problems in one run.

## Learning path

| Tutorial | Complete example | What changes |
| --- | --- | --- |
| [1. Inviscid sweep](./inviscid-naca) | `01_inviscid_sweep.py` | Reuse one model across angles of attack |
| [2. Viscous analysis](./viscous-analysis) | `02_viscous_analysis.py` | Activate the boundary-layer solve with `Re > 0` |
| [3. Geometry modifications](./geometry-mods) | `03_geometry_modifications.py` | Change CST thickness and camber |
| [4. Custom coordinates](./custom-airfoil) | `04_custom_airfoil.py` | Replace the geometry server with normalized Python coordinates |
| [5. Sensitivities](./sensitivities) | `05_sensitivities.py` | Check the analytic CST-to-aerodynamics derivative chain |
| [6. HDF5 history](./hdf5-output) | `06_hdf5_history.py` | Record server-side evaluations |

All complete files live under
[`examples/tutorials`](https://github.com/chrislupp/Scirocco/tree/develop/examples/tutorials).
Each page explains the important section, but the complete file is the source
of truth when you want to run the example.
