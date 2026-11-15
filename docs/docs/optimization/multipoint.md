---
sidebar_position: 5
---

# Multipoint design

One CST shape, three operating points: low, mid, and high lift. The
objective is a weighted average of their drag. This builds on
[CST design](./cst-design); what's new is one analysis server per operating
point, all fed by one geometry server.

This is
[`kulfan_multipoint_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py).

## Problem

| Element | Choice |
| --- | --- |
| Geometry | One shared Kulfan CST (8 upper + 8 lower) |
| Operating points | 3 design points (varying target $C_l$, $Re$) |
| Per-point DV | One $\alpha_i$ per design point |
| Shared DVs | 16 CST coefficients |
| Objective | $\sum_i w_i\, C_{d,i}$ (weighted drag) |
| Constraints | $C_{l,i} = C_{l,\text{tgt}_i}$ for each point; area floor; CST smoothness |

All points share one shape, so there's one geometry server. Each point keeps
its own solver state, so each gets its own analysis server.

## How it fits together

One geometry feeds an analysis per design point; analyses 1, 2, and 3 are the
low, mid, and high lift points. A last component combines their drag into the
objective:

<N2 name="multipoint" alt="N2 diagram: IPOPT drives one CST geometry, AirfoilArea, three Analysis components, and a weighted-drag component." />

The geometry server listens on 50060 and the analysis servers on 50061
upward, one per entry in `DESIGN_POINTS`. The script starts and stops all
of them.

## Design-point configuration

The example declares three points as a list of tuples
`(name, target_cl, Re, Ma, weight)`:

```python
DESIGN_POINTS = [
    ("low_cl",  0.30, 1.0e6, 0.0, 1.0),
    ("mid_cl",  0.50, 1.0e6, 0.0, 1.5),
    ("high_cl", 0.80, 1.0e6, 0.0, 1.0),
]
```

The script normalizes the raw weights, giving
`[0.285714, 0.428571, 0.285714]` for this configuration.

Add a point to the list and you get another server, another $\alpha$
variable, another lift constraint, and another drag term.

Before optimizing, the script uses each server's analytic $dC_l/d\alpha$ to
Newton-iterate the starting angles onto their lift targets: about
`[2.8988°, 4.5334°, 6.7619°]` for the default shape. That only improves the
starting point.

## The model

```python
geometry_channel = grpc.insecure_channel("localhost:50060")
analysis_channels = [
    grpc.insecure_channel(f"localhost:{50061 + i}")
    for i in range(len(DESIGN_POINTS))
]

model.add_subsystem("geometry",
    RemoteExplicitComponent(channel=geometry_channel, **geom_kwargs),
    promotes_outputs=["x", "smoothness_upper", "smoothness_lower"],
    promotes_inputs=["cst_upper", "cst_lower"])

for i, (name, _cl, Re, Ma, _w) in enumerate(DESIGN_POINTS):
    model.add_subsystem(f"analysis_{i}",
        RemoteExplicitComponent(channel=analysis_channels[i],
                                naca="0012", npanel=NPANEL,
                                hdf5_file=str(hdf5_for(name))),
        promotes_inputs=["x"])

# Per-point alpha and Cl-equality constraint
for i, (name, tgt_cl, *_) in enumerate(DESIGN_POINTS):
    model.add_design_var(f"analysis_{i}.alpha",
                         lower=-2, upper=10, ref=5)
    model.add_constraint(f"analysis_{i}.cl", equals=tgt_cl, ref=tgt_cl)

# Weighted-drag aggregator
weights = np.array([w for *_, w in DESIGN_POINTS])
weights = weights / weights.sum()
cd_agg = om.ExecComp("cd_total = " + " + ".join(
    f"w_{i}*cd_{i}" for i in range(N_POINTS)),
    **{f"w_{i}": weights[i] for i in range(N_POINTS)},
    **{f"cd_{i}": 0.0 for i in range(N_POINTS)})
model.add_subsystem("cd_aggregator", cd_agg)
for i in range(N_POINTS):
    model.connect(f"analysis_{i}.cd", f"cd_aggregator.cd_{i}")

model.add_objective("cd_aggregator.cd_total", ref=1e-2)
```

The full script adds setup and cleanup.

## Per-point HDF5

The example writes one HDF5 file per design point:

```
kulfan_multipoint_optimization_history_low_cl.h5
kulfan_multipoint_optimization_history_mid_cl.h5
kulfan_multipoint_optimization_history_high_cl.h5
```

Each file holds the solves of one point's server. Derivative calls aren't
recorded, but line-search trials are, so record counts don't match IPOPT's
iteration count. All three files are written in step, so you can plot them
together:

```python
import h5py, matplotlib.pyplot as plt

for name in ("low_cl", "mid_cl", "high_cl"):
    with h5py.File(f"...history_{name}.h5") as f:
        keys = sorted(
            (k for k in f if k.startswith("iteration_")),
            key=lambda k: int(k.rsplit("_", 1)[1]),
        )
        cd = [f[f"{k}/coefficients/cd"][()] for k in keys]
    plt.plot(cd, label=name)
plt.xlabel("model evaluation"); plt.ylabel("Cd"); plt.legend(); plt.show()
```

This shows the multipoint trade: drag can't fall at every point at once.
The mid-lift point has the largest weight, so it gets the most attention.

## Running it

Build the servers with `SCIROCCO_HDF5_OUTPUT=ON` to produce the histories
described above; see the [HDF5 tutorial](../python/hdf5-history).

This example uses IPOPT through ModOpt. Install IPOPT first:

```bash
# macOS (Homebrew)
brew install ipopt pkg-config

# Ubuntu
sudo apt install coinor-libipopt-dev pkg-config
```

Then its own pinned Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r examples/kulfan_multipoint_optimization/requirements.txt
```

Then launch the optimization:

```bash
python examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py
```

Each point adds a viscous solve and an adjoint per iteration, so the run
time grows with the number of points. On an Apple M1 Pro the reference run
takes 168 IPOPT iterations and about 45 seconds, including server start-up:

```output
Return status: Solve_Succeeded
low_cl:   alpha=-1.007°  Cl=0.3000  Cd=0.00434
mid_cl:   alpha=+0.362°  Cl=0.5000  Cd=0.00446
high_cl:  alpha=+3.098°  Cl=0.8000  Cd=0.00508
Σ w·Cd = 0.004602
constraint violation = 3.1e-12
```

A few trial steps early in the run produce shapes the viscous solve can't
converge. The server reports the failure, the log shows a traceback, and
IPOPT treats it as a rejected step and shortens it. Tracebacks of the form
`Scirocco analysis did not converge` are expected and don't stop the run.

Results from before Philote-Cpp 0.6.0 differ: older versions silently
truncated large Jacobians, so they aren't a valid comparison.

## Exploring the weights

To map the trade-off, rerun with different weights:

```python
for w in [(1, 0, 0), (0.5, 0.5, 0), (0.33, 0.33, 0.33), ...]:
    DESIGN_POINTS = [(n, cl, re, ma, wi) for (n, cl, re, ma, _), wi in zip(BASE, w)]
    # ... rerun optimization, save (cd_low, cd_mid, cd_high) ...
```

Each weight vector is a separate optimization of about a minute, so a
sweep of a few dozen weight vectors finishes in well under an hour.
