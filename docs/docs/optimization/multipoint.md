---
sidebar_position: 5
---

# Multipoint design

A single CST geometry optimized against three operating points
(low-, mid-, high-Cl) with a weighted-drag objective. This is the
showcase example for "real" airfoil design: the same wing has to perform
across a flight envelope, not just at one point.

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
                       ┌──→ [Analysis @ 50051] ──cl_0, cd_0──→
[Kulfan geometry] ──x──┼──→ [Analysis @ 50053] ──cl_1, cd_1──→  [SLSQP]
   (50052)             └──→ [Analysis @ 50054] ──cl_2, cd_2──→
```

One geometry server (50052), N analysis servers (50051, 50053, 50054 by
default — N is whatever `len(DESIGN_POINTS)` evaluates to). The script
starts all of them, runs the optimization, and tears them down.

## Design-point configuration

The example declares three points as a list of tuples
`(name, target_cl, Re, Ma, weight)`:

```python
DESIGN_POINTS = [
    ("low_cl",  0.30, 1.0e6, 0.0, 0.30),
    ("mid_cl",  0.50, 1.0e6, 0.0, 0.40),
    ("high_cl", 0.80, 1.0e6, 0.0, 0.30),
]
```

Editing this list changes the entire optimization shape — add a point and
you get another analysis server, another $\alpha$ DV, another lift-equality
constraint, and an extra term in the weighted drag objective.

## OpenMDAO wiring (sketch)

```python
geometry_channel = grpc.insecure_channel("localhost:50052")
analysis_channels = [
    grpc.insecure_channel(f"localhost:{50051 + 2*i}")
    for i in range(len(DESIGN_POINTS))
]

model.add_subsystem("geometry",
    RemoteExplicitComponent(channel=geometry_channel, **geom_kwargs),
    promotes_outputs=["x", "smoothness_upper", "smoothness_lower"],
    promotes_inputs=["cst_upper", "cst_lower"])

for i, (name, _cl, Re, Ma, _w) in enumerate(DESIGN_POINTS):
    model.add_subsystem(f"analysis_{name}",
        RemoteExplicitComponent(channel=analysis_channels[i],
                                naca="0012", npanel=NPANEL,
                                hdf5_file=str(hdf5_for(name))),
        promotes_inputs=["x"])

# Per-point alpha and Cl-equality constraint
model.add_design_var("alpha", lower=-5, upper=15, shape=N_POINTS, ref=10)
for i, (name, tgt_cl, *_) in enumerate(DESIGN_POINTS):
    model.connect(f"alpha[{i}]", f"analysis_{name}.alpha")
    model.add_constraint(f"analysis_{name}.cl", equals=tgt_cl, ref=tgt_cl)

# Weighted-drag aggregator
weights = np.array([w for *_, w in DESIGN_POINTS])
cd_agg = om.ExecComp("cd_avg = " + " + ".join(
    f"{w}*cd_{i}" for i, w in enumerate(weights)),
    **{f"cd_{i}": 0.0 for i in range(N_POINTS)})
model.add_subsystem("cd_aggregator", cd_agg)
for i, (name, *_) in enumerate(DESIGN_POINTS):
    model.connect(f"analysis_{name}.cd", f"cd_aggregator.cd_{i}")

model.add_objective("cd_aggregator.cd_avg", ref=1e-2)
```

(The actual script is somewhat tidier — see the source — but this is the
shape.)

## Per-point HDF5

The example writes one HDF5 file per design point:

```
kulfan_multipoint_optimization_history_low_cl.h5
kulfan_multipoint_optimization_history_mid_cl.h5
kulfan_multipoint_optimization_history_high_cl.h5
```

Each file contains the iteration history seen by **that point's** analysis
server. They're indexed by the same iteration count, so you can stack them
to plot per-point drag vs. SLSQP iteration:

```python
import h5py, matplotlib.pyplot as plt

for name in ("low_cl", "mid_cl", "high_cl"):
    with h5py.File(f"...history_{name}.h5") as f:
        cd = [f[f"iterations/{k}/post/cd"][()] for k in sorted(f["iterations"])]
    plt.plot(cd, label=name)
plt.xlabel("SLSQP iteration"); plt.ylabel("Cd"); plt.legend(); plt.show()
```

The classical multipoint trade is visible here: the optimizer can't pull
$C_d$ down at every point simultaneously. With $w = [0.3, 0.4, 0.3]$, the
mid-Cl point benefits the most; the low-Cl point might even see slightly
worse drag than a single-point design.

## Running it

```bash
python examples/kulfan_multipoint_optimization/kulfan_multipoint_optimization.py
```

Larger than the single-point examples — expect ~5–15 minutes depending on
how SLSQP polishes. The total cost grows linearly in `N_POINTS` since each
adds one independent viscous primal + adjoint per iteration. The shared
adjoint trick from the single-point case still applies within each point.

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
