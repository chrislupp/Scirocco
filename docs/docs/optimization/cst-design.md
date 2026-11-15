---
sidebar_position: 4
---

# CST design

The same kind of problem as the [NACA example](./naca-design), with a CST
parameterization: 8 coefficients per surface instead of 3 NACA parameters.
The geometry now runs on its own server, `scirocco_kulfan_server`, so the
framework (OpenMDAO in these scripts) chains two analytic Jacobians. The
adjoint's cost doesn't depend on how many coefficients there are.

There are two variants. The inviscid one trims the pitching moment; the
viscous one minimizes drag and adds a smoothness constraint.

- [`kulfan_inviscid_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_inviscid_optimization/kulfan_inviscid_optimization.py)
- [`kulfan_viscous_optimization.py`](https://github.com/chrislupp/Scirocco/blob/develop/examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py)

## How it fits together

Both scripts start both servers themselves: the geometry server
(`KulfanDiscipline`) on `localhost:50052` and the analysis server
(`SiroccoDiscipline`) on `localhost:50051`. The coefficients go to the
geometry server, its `x` goes to the analysis server, and the framework
multiplies the two Jacobians. Each variant below shows its model.

## Inviscid variant

| Element | Choice |
| --- | --- |
| Geometry | Kulfan CST, 8 upper + 8 lower Bernstein coefficients (round LE, closed TE) |
| Analysis | Scirocco inviscid ($Re = 0$, $M_\infty = 0$) |
| Design variables | 16 CST coefficients $\in [0, 0.40]$, $\alpha \in [-2°, 10°]$ |
| Objective | $\min C_m^2$ |
| Constraints | $C_l = 0.5$, area $\ge 0.8 \times$ (initial area) |

Inviscid drag is zero, so it can't be the objective. Minimizing $C_m^2$ at
a target $C_l$ gives a section with no pitching moment, a classic inviscid
design problem with a smooth minimum.

The model adds a small component that squares `cm` for the objective:

<N2 name="cst-inviscid" alt="N2 diagram: SLSQP drives CST geometry, AirfoilArea, Analysis, and a cm-squared component." />

```bash
python examples/kulfan_inviscid_optimization/kulfan_inviscid_optimization.py
```

It starts from a least-squares CST fit to NACA 0012 and writes
`kulfan_inviscid_optimization_history.h5` and a plot. A reference run
converges in 4 SLSQP iterations to `alpha = 4.3204°`, `Cl = 0.500000`,
`Cm = -1e-6`, and `area = 0.082185`.

The `lower=0.0` bounds on the coefficients keep the upper surface above
the lower one at every iterate. The area constraint stops the section from
getting too thin.


## Viscous variant

| Element | Choice |
| --- | --- |
| Geometry | Same CST: 8 upper + 8 lower |
| Analysis | Scirocco viscous, $Re = 10^6$, $M_\infty = 0$ |
| Design variables | 16 CST coefficients, $\alpha$ |
| Objective | $\min C_d$ |
| Constraints | $C_l = 0.5$, area $\ge 0.8 \times$ initial, **smoothness $\in [-0.05, 0.05]$ on each surface** |

The model drops the `cm` component, and the geometry server also returns its
smoothness outputs to SLSQP:

<N2 name="cst-viscous" alt="N2 diagram: SLSQP drives CST geometry, AirfoilArea, and Analysis." />

```bash
python examples/kulfan_viscous_optimization/kulfan_viscous_optimization.py
```

A reference run converges in 86 SLSQP iterations to `alpha = 0.9910°`, `Cl = 0.500000`, `Cd = 0.004025`,
`Cm = -0.097065`, and `area = 0.065720`.

### The smoothness constraint

The `KulfanDiscipline` emits two extra outputs that are linear in the CST
coefficients:

```
smoothness_upper[i] = cst_upper[i+2] - 2*cst_upper[i+1] + cst_upper[i]
smoothness_lower[i] = cst_lower[i+2] - 2*cst_lower[i+1] + cst_lower[i]
```

These are second differences of the coefficients. Because the shape
function's second derivative is

$$
S''(\psi) = n(n-1) \sum_i D_i\, B^{n-2}_i(\psi),
$$

bounding $|D_i| \le 0.05$ keeps the shape function from wiggling. The
viscous optimization adds it as a linear constraint:

```python
prob.model.add_constraint("smoothness_upper",
                          lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX)
prob.model.add_constraint("smoothness_lower",
                          lower=-SMOOTHNESS_MAX, upper=SMOOTHNESS_MAX)
```

Without it, SLSQP pushes single coefficients to extremes to chase tiny drag
gains. The bumps that result upset the boundary-layer solve and make
transition jump between panels. The inviscid variant doesn't need the
constraint.

## Optimizer settings

Both variants use:

```python
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"]       = 1e-6
prob.driver.options["maxiter"]   = 200
```

On an Apple M1 Pro each variant finishes in under 20 seconds, including
server start-up. The viscous one takes about 80 iterations, almost all of
the time in the solves. If SLSQP improves quickly and then
stalls short of convergence, raise `maxiter` to 400.

## Reading the output

Each run writes an HDF5 history and a PNG comparing the initial and final
sections and coefficients. To plot drag over the run:

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

The history has one record per solve, which includes line-search trials, so
its count won't match SLSQP's iteration count.

<Next items={[
  { to: "/docs/optimization/multipoint", title: "Multipoint design", body: "One shape, three operating points, weighted drag." },
]} />
