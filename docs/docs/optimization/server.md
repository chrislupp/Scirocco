---
sidebar_position: 1
---

# Philote gRPC server reference

Reference page for the two gRPC servers that wrap Scirocco as MDO
disciplines. The optimization tutorials —
[NACA 4-digit](./naca-design), [CST design](./cst-design), and
[multipoint](./multipoint) — drive these servers; this page documents the
build, launch, and discipline I/O in one place.

Scirocco ships **two** Philote-Cpp gRPC servers that wrap the library as
disciplines for use inside an MDO framework like OpenMDAO:

| Binary | Discipline | What it does |
| --- | --- | --- |
| `scirocco_server` | `SiroccoDiscipline` | Takes panel coordinates `x`, $\alpha$, $Re$, $M_\infty$ → returns `cl`, `cd`, `cm` with analytic gradients |
| `scirocco_kulfan_server` | [`KulfanDiscipline`](../api/kulfan-discipline) | Takes CST coefficients `cst_upper`, `cst_lower` → returns panel coordinates `x` (and smoothness diagnostics) |

You wire them together in OpenMDAO with `philote_mdo.openmdao.RemoteExplicitComponent`.
For a NACA-parameterized design, you only need the analysis server (the
geometry comes from your own OpenMDAO component). For CST-parameterized
designs, you start both.

This page covers how to build and launch them; the OpenMDAO wiring is in
the [NACA 4-digit tutorial](./naca-design).

## Build

The Philote interface is on by default (`SCIROCCO_BUILD_PHILOTE=ON`), but
it needs Protobuf, gRPC, and the Philote-Cpp/Philote-MDO 0.8 pairing pinned in
the [installation guide](../getting-started/installation):

```bash
# After building Philote-Cpp as described in Installation:
cmake -S scirocco -B scirocco/build-philote -DCMAKE_BUILD_TYPE=Release \
      -DPhiloteCpp_DIR="$PWD/Philote-Cpp/build"
cmake --build scirocco/build-philote -j
```

You'll find the two server binaries at:

```
scirocco/build-philote/src/philote/scirocco_server
scirocco/build-philote/src/philote/scirocco_kulfan_server
```

## Launch

Each server takes one command-line argument: the address to bind. Use a
different port for each server.

```bash
# Analysis server on 50051
scirocco/build-philote/src/philote/scirocco_server localhost:50051

# Kulfan geometry server on 50052
scirocco/build-philote/src/philote/scirocco_kulfan_server localhost:50052
```

The processes stay in the foreground and log accepted connections to
stdout. In practice you start them in the background from your driver
script — every example under `/examples/` does this:

```python
import subprocess, time

proc = subprocess.Popen(
    ["scirocco/build-philote/src/philote/scirocco_server", "localhost:50051"],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
)
time.sleep(2.0)   # give gRPC time to start listening

# ... run your OpenMDAO problem against the channel ...

proc.terminate()
proc.wait()
```

## Discipline interfaces

### Analysis (`SiroccoDiscipline`)

Configured at instantiation with two options: `naca` (digit string for the
initial geometry) and `npanel` (panel count). The CST workflow overrides
`x` on every evaluation so the `naca` option only seeds the panel-count
metadata.

| Direction | Variable | Shape |
| --- | --- | --- |
| Input | `x` | `2 * (npanel + 1)` — TE-low → LE → TE-high, x-block then z-block |
| Input | `alpha` | scalar — degrees |
| Input | `Re` | scalar — 0 means inviscid |
| Input | `Ma` | scalar — freestream Mach |
| Output | `cl`, `cd`, `cm` | scalars with partials w.r.t. every input |

Optionally the discipline can write its iteration history to an HDF5 file —
pass `hdf5_file=...` to the `RemoteExplicitComponent`.

### Geometry (`KulfanDiscipline`)

| Direction | Variable | Shape |
| --- | --- | --- |
| Option | `n_cst_upper`, `n_cst_lower` | int — coefficient counts |
| Option | `npanel` | int — even; `npoint = npanel + 1` |
| Option | `te_upper`, `te_lower` | TE y-offset |
| Input | `cst_upper`, `cst_lower` | CST coefficient vectors |
| Output | `x` | flat `2 * npoint` panel coordinates |
| Output | `smoothness_upper`, `smoothness_lower` | length `n_cst_* - 2` second differences |

Full details on the [`KulfanDiscipline` API page](../api/kulfan-discipline).

## Common pitfalls

| Symptom | Likely cause |
| --- | --- |
| `Address already in use` | A previous server didn't terminate cleanly. Find and kill it: `lsof -i :50051`. |
| Client connects but every call returns `UNAVAILABLE` | gRPC bound to `localhost` but client tried `0.0.0.0` (or vice versa). Match exactly. |
| `RuntimeError: server failed to start` in driver script | Server binary path wrong, or Philote-Cpp not installed. Try running the binary by hand to see the error. |
| Optimizer iterates produce NaN `cl` | A bad CST iterate produced self-intersecting coordinates. Confirm CST coefficient bounds are correct (see [CST design](./cst-design)). |

## Next: [NACA 4-digit optimization](./naca-design)
