---
sidebar_position: 2
---

# 2. Enable viscous analysis

Continue with the model from [the inviscid sweep](./inviscid-naca). A positive
Reynolds number activates Scirocco's coupled integral boundary-layer solver;
the Python and Philote wiring does not change.

Run the complete example:

```bash
python examples/tutorials/02_viscous_analysis.py
```

## Change the operating point

Replace the inviscid sweep with:

```python
problem.set_val("analysis.alpha", 2.0, units="deg")
problem.set_val("analysis.Ma", 0.0)

print("       Re       Cl        Cd        Cm")
for reynolds in (2.5e5, 5.0e5, 1.0e6, 2.0e6):
    problem.set_val("analysis.Re", reynolds)
    problem.run_model()
    cl = problem.get_val("analysis.cl")[0]
    cd = problem.get_val("analysis.cd")[0]
    cm = problem.get_val("analysis.cm")[0]
    print(f"{reynolds:9.2e}  {cl:8.5f}  {cd:8.5f}  {cm:8.5f}")
```

The analysis server calls the same `SiroccoDiscipline::Compute` method as in
the first tutorial. The value of `Re` selects the solve path:

- `Re == 0`: linear inviscid panel solve.
- `Re > 0`: inviscid solve coupled to laminar/turbulent boundary-layer and
  wake equations.

The tested result is:

```text
       Re       Cl        Cd        Cm
 2.50e+05   0.21382   0.01123   0.00728
 5.00e+05   0.15099   0.00795   0.01536
 1.00e+06   0.13889   0.00554   0.02043
 2.00e+06   0.07474   0.00484   0.03382
```

Expect non-zero drag; for this case it decreases as Reynolds number rises.
Lift and moment also move because the boundary-layer displacement thickness
feeds back into the outer-flow solution. Use the table as a smoke-test, not as
validation data: small numerical changes can result from compiler and solver
differences.

## Add compressibility

Mach number is another ordinary remote input:

```python
problem.set_val("analysis.Re", 1.0e6)
for mach in (0.0, 0.2, 0.4, 0.6):
    problem.set_val("analysis.Ma", mach)
    problem.run_model()
    print(mach, problem.get_val("analysis.cl")[0])
```

Scirocco applies its subsonic compressibility correction internally. Keep the
tutorial cases below roughly Mach 0.7; this is not a transonic solver.

## Keep the model for the next tutorial

The next page changes the upstream CST inputs instead of the operating-point
inputs. The two remote components and their `x` connection remain exactly the
same.

## Next: [3. Modify parameterized geometry](./geometry-mods)
