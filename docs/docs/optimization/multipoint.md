---
sidebar_position: 5
---

# Multipoint design

A single CST geometry optimized against three operating points
(low-, mid-, high-Cl) with a weighted-drag objective. This is the
showcase example for "real" airfoil design: the same wing has to perform
across a flight envelope, not just at one point.

This builds on the [CST design tutorial](./cst-design) — the geometry
discipline and constraints are identical. What's new is running **N
analysis servers**, one per operating point, against a single shared
geometry server, and aggregating their drag into one weighted objective.

Mirrors
[`examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py).

## Problem

| Element | Choice |
| --- | --- |
| Geometry | One shared Kulfan CST (8 upper + 8 lower) |
| Operating points | 3 design points (varying target $C_l$, $Re$) |
| Per-point DV | One $\alpha_i$ per design point |
| Shared DVs | 16 CST coefficients |
| Objective | $\sum_i w_i\, C_{d,i}$ (weighted drag) |
| Constraints | $C_{l,i} = C_{l,\text{tgt}_i}$ for each point; area floor; CST smoothness |

The shape is single-valued — same panel coordinates feeding every analysis
— so the geometry component drives one Kulfan server. Each operating point
needs its own analysis state, so the script launches **N independent
analysis servers**, one per design point.

## Architecture

```
                       ┌──→ [Analysis @ 50061] ──cl_0, cd_0──→
[Kulfan geometry] ──x──┼──→ [Analysis @ 50062] ──cl_1, cd_1──→  [IPOPT]
   (50060)             └──→ [Analysis @ 50063] ──cl_2, cd_2──→
```

One geometry server (50060), N analysis servers (50061 … 50060+N by
default — N is whatever `len(DESIGN_POINTS)` evaluates to). The script
starts all of them, runs the optimization, and tears them down.

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

Editing this list changes the entire optimization shape — add a point and
you get another analysis server, another $\alpha$ DV, another lift-equality
constraint, and an extra term in the weighted drag objective.

Before constructing the ModOpt wrapper, the script uses the analytic
$dC_l/d\alpha$ from each analysis discipline to Newton-match the initial
angles of attack to their lift targets. For the default shape this gives
approximately `[2.8988°, 4.5334°, 6.7619°]`. This changes only the starting
point, not the optimization problem.

## OpenMDAO wiring (sketch)

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

This excerpt shows the component connections; see the complete script for
setup and cleanup.

## Per-point HDF5

The example writes one HDF5 file per design point:

```
kulfan_multipoint_optimization_history_low_cl.h5
kulfan_multipoint_optimization_history_mid_cl.h5
kulfan_multipoint_optimization_history_high_cl.h5
```

Each file contains the evaluation history seen by **that point's** analysis
server. The records contain successful primal evaluations; derivative-only calls do
not append records. Line searches can evaluate several trial points, so these
records do not correspond one-to-one with IPOPT's iteration number. The
three files receive calls in the same order and can be stacked to compare
per-point drag across model evaluations:

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

The classical multipoint trade is visible here: the optimizer can't pull
$C_d$ down at every point simultaneously. With normalized weights of about
$[0.286, 0.429, 0.286]$, the mid-Cl point has the largest influence; the
low-Cl point might even see slightly worse drag than a single-point design.

## Running it

Build the servers with `SCIROCCO_HDF5_OUTPUT=ON` to produce the histories
described above; see the [HDF5 tutorial](../tutorials/hdf5-output).

This example adds IPOPT through the exact ModOpt revision used for release
validation. Install the system IPOPT library first:

```bash
# macOS (Homebrew)
brew install ipopt pkg-config

# Ubuntu
sudo apt install coinor-libipopt-dev pkg-config
```

Then install its separate pinned environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r examples/kulfan_multipoint_optimization/requirements.txt
```

Then launch the optimization:

```bash
python examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py
```

This is intentionally much larger than the single-point examples. The cost
grows approximately linearly in `N_POINTS` because each point adds an
independent viscous primal and adjoint evaluation.

The Philote 0.8 release-validation run completed in 222 IPOPT iterations
(about 13.5 minutes on the validation machine):

```text
Return status: Solve_Succeeded
low_cl:   alpha=-0.988°  Cl=0.3000  Cd=0.00432
mid_cl:   alpha=+0.363°  Cl=0.5000  Cd=0.00447
high_cl:  alpha=+3.098°  Cl=0.8000  Cd=0.00508
Σ w·Cd = 0.004601
constraint violation = 2.2e-11
```

This differs from results produced through Philote-Cpp v0.5.0. That release
silently truncated streamed arrays after each complete 1,000-value chunk; a
`202 × 8` Kulfan geometry Jacobian therefore lost its final 616 values.
Philote-Cpp 0.8 includes the chunking fix and sends the complete Jacobian, so
the old objective is not a valid regression reference.

## Sensitivity to the weight vector

If you want to see the Pareto front, sweep `weights` along a simplex:

```python
for w in [(1, 0, 0), (0.5, 0.5, 0), (0.33, 0.33, 0.33), ...]:
    DESIGN_POINTS = [(n, cl, re, ma, wi) for (n, cl, re, ma, _), wi in zip(BASE, w)]
    # ... rerun optimization, save (cd_low, cd_mid, cd_high) ...
```

Each weight choice is a separate optimization, so plan for an overnight run
if you want a smooth surface.

## Full source

[`examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py).
