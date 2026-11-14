---
title: Server reference
description: Launch options and inputs and outputs of the two Philote servers.
---

# Server reference

## Launch

Each server takes one argument, the address to listen on:

```bash
scirocco/build-philote/src/philote/scirocco_server        localhost:50051
scirocco/build-philote/src/philote/scirocco_kulfan_server localhost:50052
```

They run in the foreground and log connections to stdout. To have a script
own its server, start it as a subprocess and wait for the channel:

```python
import subprocess
import grpc

proc = subprocess.Popen([
    "scirocco/build-philote/src/philote/scirocco_server", "localhost:50051"
])
try:
    with grpc.insecure_channel("localhost:50051") as channel:
        grpc.channel_ready_future(channel).result(timeout=15)
        # Connect your Philote client to this channel.
finally:
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
```

The optimization examples start and stop their own servers this way.

## Analysis server

`scirocco_server` hosts the analysis discipline (`SiroccoDiscipline` in the
C++ source).

| Options | |
| --- | --- |
| `npanel` | Panel count. `x` has `2 * (npanel + 1)` entries. |
| `naca` or `coordinates` | Exactly one. A NACA string, or `{"x": [...], "z": [...]}`. Only initializes the solver; `x` supplies the shape analyzed. |
| `hdf5_file` | Optional history file, written by the server. Needs `SCIROCCO_HDF5_OUTPUT=ON`. |

| Inputs | |
| --- | --- |
| `x` | Panel nodes, all x then all z, lower trailing edge → upper trailing edge. Replaces the panels directly. |
| `alpha` | Angle of attack in degrees |
| `Re` | Reynolds number; `0` for inviscid |
| `Ma` | Mach number |

| Outputs | |
| --- | --- |
| `cl`, `cd`, `cm` | With analytic partials with respect to every input |

A failed solve raises an error instead of returning numbers. To impose a
lift coefficient, make `alpha` a design variable and constrain `cl`.

## CST geometry server

`scirocco_kulfan_server` hosts [`KulfanDiscipline`](../api/kulfan-discipline).

| Options | |
| --- | --- |
| `n_cst_upper`, `n_cst_lower` | Coefficients per surface, 3–256 |
| `npanel` | Even panel count |
| `te_upper`, `te_lower` | Trailing-edge z offsets |

| Inputs | |
| --- | --- |
| `cst_upper`, `cst_lower` | Coefficient vectors |

| Outputs | |
| --- | --- |
| `x` | Panel nodes in the analysis server's layout |
| `smoothness_upper`, `smoothness_lower` | Second differences of the coefficients, `n − 2` each |

Valid ranges and layouts are in [Conventions](../api/conventions).
