---
sidebar_position: 3
---

# 3. Modify parameterized geometry

The first two tutorials held the CST coefficients fixed. Now use the Python
bindings to drive `KulfanDiscipline` as a parameterized geometry generator and
send each new shape directly into `SiroccoDiscipline`.

Run the complete example:

```bash
python examples/tutorials/03_geometry_modifications.py
```

## Thickness and camber

Start from the same symmetric coefficient vector:

```python
base = np.array([
    0.1706635358053721, 0.15626369325824663, 0.1546327069833401,
    0.14462202838023833, 0.13318428664336585, 0.15537328876198608,
])

problem.set_val("analysis.alpha", 2.0, units="deg")
problem.set_val("analysis.Re", 1.0e6)
problem.set_val("analysis.Ma", 0.0)

cases = {
    "thin":      (0.80 * base, 0.80 * base),
    "baseline":  (base, base),
    # Upper coefficients larger than lower coefficients add positive camber.
    "cambered":  (1.12 * base, 0.88 * base),
}

for name, (upper, lower) in cases.items():
    problem.set_val("cst_upper", upper)
    problem.set_val("cst_lower", lower)
    problem.run_model()
    print(
        f"{name:9s}",
        f"Cl={problem.get_val('analysis.cl')[0]:.5f}",
        f"Cd={problem.get_val('analysis.cd')[0]:.5f}",
    )
```

The tested result is:

```text
thin      Cl=0.21032 Cd=0.00533
baseline  Cl=0.20489 Cd=0.00565
cambered  Cl=0.27794 Cd=0.00574
```

Each `run_model()` follows the same chain:

```text
(upper CST, lower CST) -> Kulfan server -> x -> analysis server -> (Cl, Cd, Cm)
```

Because Bernstein basis functions are non-negative on the chord, non-negative
coefficients keep each half-thickness non-negative. That is a useful starting
bound for design work, though practical optimizations also need thickness,
area, or smoothness constraints.

## Inspect the coordinates in Python

The geometry output is already available in the Python process:

```python
x_flat = problem.get_val("x")
nnode = npanel + 1
x_coord = x_flat[:nnode]
z_coord = x_flat[nnode:]

import matplotlib.pyplot as plt
plt.plot(x_coord, z_coord, "o-")
plt.axis("equal")
plt.xlabel("x/c")
plt.ylabel("z/c")
plt.show()
```

Install the tutorial requirements before using this optional plot; they
include Matplotlib.

This flattened `[x; z]` layout is the contract between the two Philote
disciplines. It is also the layout to use when a local Python geometry model
replaces `KulfanDiscipline`, which is the subject of the next tutorial.

## What about `GeomFlap` and `GeomAddCamber`?

Those are direct C++ library operations, not inputs on the current Philote
discipline. In a Python/Philote workflow, represent shape changes through CST
coefficients or create a local OpenMDAO geometry component that outputs `x`.
The [C++ geometry API](../api/geometry) documents the in-process alternatives.

## Next: [4. Supply coordinates from Python](./custom-airfoil)
