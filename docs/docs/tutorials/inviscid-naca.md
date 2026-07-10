---
sidebar_position: 1
---

# 1. Sweep an inviscid airfoil

This tutorial extends the [Python quickstart](../getting-started/quickstart).
Keep both Philote servers running and reuse its `quickstart.py` model. The
geometry is a symmetric CST approximation to NACA 0012; the analysis is
inviscid because `Re` is zero.

By the end you will know how Philote variables map into OpenMDAO, how to reuse
one remote model for several operating points, and why inviscid drag is zero.

## Replace the single evaluation with a sweep

Everything through `problem.setup()` stays the same. Replace the input and
output section with:

```python
cst_0012 = np.array([0.170, 0.160, 0.145, 0.145, 0.160, 0.170])
problem.set_val("cst_upper", cst_0012)
problem.set_val("cst_lower", cst_0012)
problem.set_val("analysis.Re", 0.0)
problem.set_val("analysis.Ma", 0.0)

print(" alpha       Cl          Cd          Cm")
for alpha in np.arange(-4.0, 12.1, 2.0):
    problem.set_val("analysis.alpha", alpha, units="deg")
    problem.run_model()
    cl = problem.get_val("analysis.cl")[0]
    cd = problem.get_val("analysis.cd")[0]
    cm = problem.get_val("analysis.cm")[0]
    print(f"{alpha:6.1f}  {cl:9.5f}  {cd:10.3e}  {cm:9.5f}")
```

`problem.run_model()` sends the current inputs to both remote components. The
Kulfan server evaluates the same CST geometry, then the analysis server solves
it at the new angle of attack. Reusing the channels and `Problem` avoids
repeating component discovery and setup.

## The interface you are driving

`RemoteExplicitComponent` discovers the variables declared by each C++
discipline:

| Component | Inputs | Outputs |
| --- | --- | --- |
| `KulfanDiscipline` | `cst_upper`, `cst_lower` | `x`, `smoothness_upper`, `smoothness_lower` |
| `SiroccoDiscipline` | `x`, `alpha`, `Re`, `Ma` | `cl`, `cd`, `cm` |

The promoted `x` variable connects them. It is a flat array containing all
node x-coordinates followed by all node z-coordinates. The nodes run from the
lower trailing edge, around the leading edge, to the upper trailing edge.

The `naca="0012"` analysis option initializes the remote solver, while the
connected `x` input supplies the evaluated geometry. `npanel` must agree on
both components (`npanel + 1` nodes).

## Interpret the result

Lift should be nearly linear with angle of attack. Drag should be zero to
solver precision: a two-dimensional body in inviscid, irrotational flow has no
profile drag (D'Alembert's paradox). The next tutorial changes only one input
to activate boundary-layer coupling.

## Next: [2. Enable viscous analysis](./viscous-analysis)
