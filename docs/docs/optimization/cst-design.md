---
sidebar_position: 4
---

# CST design

The Kulfan / class-shape transformation (CST) parameterization replaces
NACA's 3 design variables with 16 (8 upper + 8 lower Bernstein coefficients
by default) — enough freedom to do real design, while still cheap enough
that the adjoint is the only sensible way to compute gradients.

This builds directly on the [NACA 4-digit tutorial](./naca-design) — read
that first for the server build, Python dependencies, OpenMDAO wiring, and
derivative-check basics. The new ingredient here is a **second** Philote
server: instead of a local Python geometry component, the CST geometry is
its own remote discipline (`scirocco_kulfan_server`), so the optimizer
chains gradients through two analytic Jacobians.

This page covers both the **inviscid** ($Re = 0$) and **viscous** ($Re > 0$)
variants. They share architecture; the viscous one adds a smoothness
constraint to keep the BL solver happy on intermediate iterates.

Sources:
- [`examples/kulfan_inviscid_optimization/kulfan_inviscid_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_inviscid_optimization/kulfan_inviscid_optimization.py)
- [`examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py)

## Architecture

Two Philote-Cpp servers, both started by the example script:

```
[KulfanDiscipline] --x--> [SiroccoDiscipline] --cl, cd, cm--> [SLSQP]
   (Philote gRPC          (Philote gRPC                ^
    localhost:50052)        localhost:50051)           |
                       +--> [AirfoilArea] -------------+
```

The CST coefficients (`cst_upper`, `cst_lower`) are SLSQP design variables;
they go into the Kulfan server, which emits panel coordinates `x`, which go
into the analysis server. OpenMDAO chains the analytic Jacobians.

## Inviscid variant

| Element | Choice |
| --- | --- |
| Geometry | Kulfan CST, 8 upper + 8 lower Bernstein coefficients (round LE, closed TE) |
| Analysis | Scirocco inviscid ($Re = 0$, $M_\infty = 0$) |
| Design variables | 16 CST coefficients $\in [0, 0.40]$, $\alpha \in [-2°, 10°]$ |
| Objective | $\min C_m^2$ |
| Constraints | $C_l = 0.5$, area $\ge 0.8 \times$ (initial area) |

**Why $C_m^2$ and not $C_d$?** In 2-D potential flow $C_d \equiv 0$, so it's
useless as an objective. Minimizing $C_m^2$ at a target $C_l$ produces a
naturally pitch-trimmed airfoil — a classical inviscid design problem with
a real smooth minimum.

```bash
python examples/kulfan_inviscid_optimization/kulfan_inviscid_optimization.py
```

The script:

1. Fits a NACA-0012-equivalent initial guess to the CST basis via
   least squares.
2. Starts both Philote servers.
3. Runs SLSQP with the design variables and constraints above.
4. Saves `kulfan_inviscid_optimization_history.h5` and a result plot.

With the pinned Philote 0.8 tutorial environment, a representative run
converges in 4 SLSQP iterations to approximately `alpha = 4.3204°`,
`Cl = 0.500000`, `Cm = -1e-6`, and `area = 0.082185`.

### Bernstein non-negativity

A nice property of the CST basis: setting `cst_upper[i] >= 0` and
`cst_lower[i] >= 0` for all $i$ is **sufficient** to guarantee that the
airfoil is non-self-intersecting at every intermediate optimizer
iterate — no per-evaluation thickness constraint needed. This is what
the `lower=0.0` design-variable bounds are doing.


## Viscous variant

| Element | Choice |
| --- | --- |
| Geometry | Same CST: 8 upper + 8 lower |
| Analysis | Scirocco viscous, $Re = 10^6$, $M_\infty = 0$ |
| Design variables | 16 CST coefficients, $\alpha$ |
| Objective | $\min C_d$ |
| Constraints | $C_l = 0.5$, area $\ge 0.8 \times$ initial, **smoothness $\in [-0.05, 0.05]$ on each surface** |

```bash
python examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py
```

A representative pinned-environment run converges in 86 SLSQP iterations
to approximately `alpha = 0.9910°`, `Cl = 0.500000`, `Cd = 0.004025`,
`Cm = -0.097065`, and `area = 0.065720`.

### The smoothness constraint

The `KulfanDiscipline` emits two extra outputs that are linear in the CST
coefficients:

```
smoothness_upper[i] = cst_upper[i+2] - 2*cst_upper[i+1] + cst_upper[i]
smoothness_lower[i] = cst_lower[i+2] - 2*cst_lower[i+1] + cst_lower[i]
```

These are the second differences of the coefficient sequences. Because the
shape function's second derivative

$$
S''(\psi) = n(n-1) \sum_i D_i\, B^{n-2}_i(\psi),
$$

bounding $|D_i| \le 0.05$ directly caps the curvature of the surface. In
the viscous optimization that's the constraint:

```python
prob.model.add_constraint("smoothness_upper",
                          lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX)
prob.model.add_constraint("smoothness_lower",
                          lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX)
```

Why this constraint matters: without it, SLSQP loves to push individual
coefficients to extreme values to chase tiny drag reductions, producing
geometry with sharp bumps that destabilize the viscous BL Newton solve
(or move transition by several panels between iterations, with the same
effect). The smoothness cap fixes both at once.

For a fully inviscid run the smoothness constraint isn't needed — the
inviscid solver doesn't care about local curvature beyond the basic
non-self-intersecting requirement.

## Optimizer hygiene

Both variants use these settings:

```python
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"]       = 1e-6
prob.driver.options["maxiter"]   = 200
```

With 16 DVs SLSQP runs through quickly — under a minute for the inviscid
case, a few minutes for the viscous one (limited by primal solve cost,
not optimizer overhead).

If you see SLSQP make a few rapid improvements then plateau, try
increasing `maxiter` to 400 — the last bit of polishing on the optimal
$C_m^2$ or $C_d$ often takes more iterations than you expect.

## Reading the output

Each run writes an HDF5 history file and a PNG plot. The HDF5 contains
the per-iteration primal state you can post-process; the PNG
overlays the initial and optimized airfoils plus the converged CST
coefficient values.

```python
import h5py
with h5py.File("kulfan_viscous_optimization_history.h5") as f:
    keys = sorted(
        (k for k in f if k.startswith("iteration_")),
        key=lambda k: int(k.rsplit("_", 1)[1]),
    )
    cds = [f[f"{k}/coefficients/cd"][()] for k in keys]
import matplotlib.pyplot as plt
plt.semilogy(cds); plt.xlabel("iter"); plt.ylabel("Cd"); plt.show()
```

The history records remote model evaluations, including evaluations used for
derivatives, rather than mapping one-to-one to the optimizer's displayed
iteration number.

## Next: [Multipoint design](./multipoint) — minimize a weighted-average drag over multiple operating points.
