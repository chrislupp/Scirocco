---
sidebar_position: 3
---

# NACA 4-digit optimization

This tutorial turns the analysis model from the
[six-part Python/Philote sequence](../tutorials/inviscid-naca) into an
optimization: minimize $C_d$ at a target $C_l$ with a section-area floor,
using a NACA 4-digit parameterization (3 design variables) plus $\alpha$.
The model reuses the remote analysis, total-derivative, and HDF5 patterns from
those tutorials. The [CST](./cst-design) and [multipoint](./multipoint)
tutorials then build on this optimization.

This mirrors
[`examples/naca_optimization/naca_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/naca_optimization/naca_optimization.py).

## How the pieces fit together

Scirocco plugs into [OpenMDAO](https://openmdao.org/) as a discipline over
gRPC. The library is wrapped in a Philote-Cpp server (`scirocco_server`);
OpenMDAO talks to it through `philote_mdo.openmdao.RemoteExplicitComponent`,
which presents the discipline as a standard `ExplicitComponent` with
analytic partials. For a NACA-parameterized design you only need this one
analysis server — the geometry is a small local Python component. (The
CST tutorials add a second server for the geometry.)

The full server build and discipline interface are documented on the
[Philote gRPC server reference](./server) page; the essentials are below.

## Prerequisites

**Build Scirocco with Philote on** (the default). This needs Protobuf,
gRPC, and the tested Philote-Cpp 0.8 build tree — see the
[server reference](./server#build) for the one-time Philote-Cpp build:

```bash
cmake -S scirocco -B scirocco/build-philote \
  -DCMAKE_BUILD_TYPE=Release \
  -DPhiloteCpp_DIR=/path/to/Philote-Cpp/build
cmake --build scirocco/build-philote -j
```

**Install the Python dependencies:**

```bash
python -m pip install -r examples/tutorials/requirements.txt
```

You do not need to launch the server by hand — every example script starts
`scirocco_server` in the background and tears it down on exit.

## Problem

| Element | Choice |
| --- | --- |
| Geometry | NACA 4-digit: camber $m$ (%), location $p$ (10ths), thickness $t$ (%) |
| Analysis | Scirocco viscous, $Re = 10^6$, $M_\infty = 0$ |
| Design variables | $m \in [0, 9]$, $p \in [1, 9]$, $t \in [1, 40]$, $\alpha \in [-5°, 15°]$ |
| Objective | $\min C_d$ |
| Constraints | $C_l = 0.5$ (equality), area $\ge 0.8 \times$ (initial area) |
| Optimizer | SciPy SLSQP, `tol = 1e-6`, `maxiter = 200` |

The geometry component (`NACA4Geometry`) is implemented locally in Python
because it's only 3 design variables and the analytic partials of the
4-digit thickness+camber formulas are straightforward.

## Architecture

```
[NACA4Geometry] --x--> [Scirocco (Philote)] --cl, cd, cm--> [SLSQP]
   analytical    |        analytical adjoint           ^
   gradients     |          sensitivities              |
                 +--> [AirfoilArea] --area-------------+
```

One server (`scirocco_server` on `localhost:50051`); no Kulfan server
needed because geometry is built in Python.

## Running it

```bash
# Build with Philote on (default), if you haven't already
cmake -S scirocco -B scirocco/build-philote \
  -DCMAKE_BUILD_TYPE=Release \
  -DPhiloteCpp_DIR="$PWD/Philote-Cpp/build"
cmake --build scirocco/build-philote -j

# Run — the script starts and stops the server itself
python examples/naca_optimization/naca_optimization.py
```

You should see roughly:

```text
Optimization: minimize Cd, target Cl = 0.5
  ...
Optimization Results
  camber     = 2.5659   (m = 0.0257)
  camber_loc = 5.4208   (p = 0.5421)
  thickness  = 9.5881   (t = 0.0959)
  alpha      = 1.3204 deg
  Cl = 0.500000
  Cd = 0.004446
  Cm = -0.085310
  Area = 0.065725  (min = 0.065725)
```

The driver prefers thin, well-cambered airfoils at small $\alpha$ — exactly
what you'd expect for "minimum drag at fixed lift". The area floor binds
(area equals its lower bound) because thinner is better for friction-
dominated drag at this $Re$.

## OpenMDAO wiring (annotated)

The script does three passes — a single-point evaluation, a derivative
check, and the actual optimization. The optimization function is the
load-bearing part:

```python
model.add_subsystem(
    "geometry", NACA4Geometry(n_nodes=nn),
    promotes_outputs=["x"],
)
model.add_subsystem(
    "area_comp", AirfoilArea(n_nodes=nn),
    promotes_inputs=["x"],
)
model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(
        channel=channel, naca="0012", npanel=npanel,
        hdf5_file=hdf5_path,         # log every call to opt_history.h5
    ),
    promotes_inputs=["x"],
)

# Design variables
model.add_design_var("geometry.camber",     lower=0,   upper=9,  ref=9)
model.add_design_var("geometry.camber_loc", lower=1,   upper=9,  ref=9)
model.add_design_var("geometry.thickness",  lower=1,   upper=40, ref=40)
model.add_design_var("analysis.alpha",      lower=-5,  upper=15, ref=10)

# Objective and constraints
model.add_objective ("analysis.cd",       ref=1e-2)
model.add_constraint("analysis.cl",       equals=target_cl, ref=target_cl)
model.add_constraint("area_comp.area",    lower=0.8 * init_area, ref=init_area)
```

A handful of things worth understanding:

| Snippet | Why |
| --- | --- |
| `promotes_outputs=["x"]`, `promotes_inputs=["x"]` | Lets `x` flow between geometry, area, and analysis without naming collisions. |
| `RemoteExplicitComponent(..., naca="0012", npanel=npanel)` | The `naca` option seeds the analysis with a 0012 baseline at start-up. Once the optimizer is running, `x` is overridden every iteration, so the seed never matters for results — but `npanel` does, because it sets the discipline's expected input shape. |
| `hdf5_file=hdf5_path` | Records every primal + sensitivity call to disk for post-analysis. |
| `ref=...` scaling | SLSQP works better in $O(1)$ — these scale the internal representation to roughly unit magnitude without changing the math. |

## Verifying derivatives

Before kicking off a long optimization, sanity-check the gradients that
flow through the discipline. OpenMDAO's `check_totals` runs
finite-difference perturbations on the inputs and compares them to the
chained analytic Jacobian:

```python
prob.run_model()
prob.check_totals(
    of=["analysis.cl", "analysis.cd", "analysis.cm"],
    wrt=["geometry.camber", "geometry.camber_loc",
         "geometry.thickness", "analysis.alpha"],
    compact_print=True,
)
```

Expect relative errors around $10^{-5}$–$10^{-7}$ for a smooth converged
primal. Larger errors usually mean either the primal isn't converged at
the linearization point, or the design is sitting on a discontinuity
(transition jumping panels, near-stall). The `naca_optimization.py` script
runs this check automatically before the optimization pass.

## What's actually being differentiated

Every SLSQP iteration triggers one residual evaluation plus one total-
derivative computation. The chain that gets differentiated:

$$
\frac{\partial C_d}{\partial \{m, p, t, \alpha\}}
   = \underbrace{\frac{\partial C_d}{\partial x}}_{\text{Scirocco adjoint}}
     \cdot
     \underbrace{\frac{\partial x}{\partial \{m, p, t\}}}_{\text{NACA4Geometry analytic Jac}}
   + \frac{\partial C_d}{\partial \alpha}
$$

OpenMDAO assembles the two pieces. Scirocco delivers the first one in
constant cost (the adjoint) independent of how many geometry DVs you have.
Adding more DVs (e.g. trying the [CST workflow](./cst-design) with 16 DVs)
costs the same single adjoint solve.

## Full source

[`examples/naca_optimization/naca_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/naca_optimization/naca_optimization.py).

## Next: [CST design](./cst-design) — same problem structure, richer parameterization.
