---
sidebar_position: 3
---

# NACA 4-digit optimization

Minimize drag at a fixed lift coefficient, with an area floor, by changing a
NACA 4-digit section's camber, camber position, and thickness, plus the
angle of attack. The geometry is a small Python component; the analysis
runs on `scirocco_server`.

This is
[`examples/naca_optimization/naca_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/naca_optimization/naca_optimization.py).
It starts and stops its own server. You need the servers built with HDF5
([Installation](../getting-started/installation#philote-servers-for-python))
and the [Python environment](../python#set-up).

## Problem

| Element | Choice |
| --- | --- |
| Geometry | NACA 4-digit: camber $m$ (%), location $p$ (10ths), thickness $t$ (%) |
| Analysis | Scirocco viscous, $Re = 10^6$, $M_\infty = 0$ |
| Design variables | $m \in [0, 9]$, $p \in [1, 9]$, $t \in [1, 40]$, $\alpha \in [-5°, 15°]$ |
| Objective | $\min C_d$ |
| Constraints | $C_l = 0.5$ (equality), area $\ge 0.8 \times$ (initial area) |
| Optimizer | SciPy SLSQP, `tol = 1e-6`, `maxiter = 200` |

The geometry component, `NACA4Geometry`, is plain Python with analytic
partials of the 4-digit formulas.

## How it fits together

```text
[NACA4Geometry] --x--> [Scirocco (Philote)] --cl, cd, cm--> [SLSQP]
   analytical    |        analytical adjoint           ^
   gradients     |          sensitivities              |
                 +--> [AirfoilArea] --area-------------+
```

Only the analysis server is needed, on `localhost:50051`.

## Run it

```bash
python examples/naca_optimization/naca_optimization.py
```

The result:

```output
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

The optimizer goes for a thinner, more cambered section at a small angle of
attack. The area floor is active: at this Reynolds number drag is mostly
friction, and a thinner section has less of it.

## The model

The script runs one evaluation, checks the derivatives, then optimizes. The
optimization sets up the model like this:

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

| Line | Why |
| --- | --- |
| `promotes_outputs=["x"]`, `promotes_inputs=["x"]` | Lets `x` flow between geometry, area, and analysis without naming collisions. |
| `naca="0012", npanel=npanel` | `naca` only initializes the server; `x` replaces the shape on every evaluation. `npanel` fixes the size of `x`. |
| `hdf5_file=hdf5_path` | Records successful primal evaluations for post-analysis; derivative calls do not append records. |
| `ref=...` | Scales each variable to order one, which SLSQP needs to behave well |

## Check the derivatives first

Before a long optimization, compare the chained analytic derivatives with
finite differences:

```python
prob.run_model()
prob.check_totals(
    of=["analysis.cl", "analysis.cd", "analysis.cm"],
    wrt=["geometry.camber", "geometry.camber_loc",
         "geometry.thickness", "analysis.alpha"],
    compact_print=True,
)
```

Expect relative errors of $10^{-5}$ to $10^{-7}$. Larger errors usually mean
the solve hasn't converged, or the design sits where transition jumps
between panels or near stall. The script runs this check before it
optimizes.

## What gets differentiated

Each SLSQP iteration needs one solve and one total derivative:

$$
\frac{\partial C_d}{\partial \{m, p, t, \alpha\}}
   = \underbrace{\frac{\partial C_d}{\partial x}}_{\text{Scirocco adjoint}}
     \cdot
     \underbrace{\frac{\partial x}{\partial \{m, p, t\}}}_{\text{NACA4Geometry analytic Jac}}
   + \frac{\partial C_d}{\partial \alpha}
$$

The framework multiplies the two pieces. Scirocco's adjoint costs the same no
matter how many geometry variables feed `x`, so moving to the 16 variables
of the [CST example](./cst-design) doesn't make the gradient more expensive.

<Next items={[
  { to: "/docs/optimization/cst-design", title: "CST design", body: "The same problem with 16 shape variables." },
]} />
