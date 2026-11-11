---
title: 3. Change the shape
description: Change thickness and camber through CST coefficients.
---

# 3. Change the shape

The CST coefficients are ordinary inputs, so changing the shape is a
`set_val` away.

```bash
python examples/tutorials/03_geometry_modifications.py
```

## Thickness and camber

Scale both surfaces to thin the section. Make the upper coefficients larger
than the lower ones to add camber:

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

```output
thin      Cl=0.21032 Cd=0.00533
baseline  Cl=0.20489 Cd=0.00565
cambered  Cl=0.27794 Cd=0.00574
```

Keep coefficients non-negative: then the upper surface can't cross the
lower one. They don't stop the section getting arbitrarily thin, which is
why the optimization examples add an area constraint.

## Plot the shape

The coordinates are an output, so you can look at them directly:

```python
import matplotlib.pyplot as plt

x_flat = problem.get_val("x")
nnode = npanel + 1
plt.plot(x_flat[:nnode], x_flat[nnode:], "o-")
plt.axis("equal")
plt.xlabel("x/c")
plt.ylabel("z/c")
plt.show()
```

## Flaps and camber lines

The C++ methods `GeomFlap` and `GeomAddCamber` aren't available through
the server. In Python, express shape changes through CST coefficients, or
write a component that produces `x` yourself, as the next tutorial does.

<Next items={[
  { to: "/docs/python/custom-coordinates", title: "4. Your own coordinates", body: "Replace the geometry server with Python." },
]} />
