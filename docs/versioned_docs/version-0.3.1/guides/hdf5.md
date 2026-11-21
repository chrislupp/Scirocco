---
title: Save results to HDF5
description: Record every solve, with distributions and gradients, to an HDF5 file.
---

# Save results to HDF5

`Hdf5Writer` appends one group per solve to an HDF5 file: operating point,
coefficients, panel geometry, surface and boundary-layer distributions, and,
if you pass them, the gradients. It's built when you configure with
`-DSCIROCCO_HDF5_OUTPUT=ON`; link against `scirocco_hdf5`.

```cpp
#include <stdexcept>
#include <scirocco/io/hdf5_writer.h>
#include <scirocco/scirocco.h>

int main() {
    scirocco::Scirocco foil("2412", 160);
    scirocco::Hdf5Writer writer("alpha_sweep.h5");   // truncates an existing file

    for (double alpha = -2.0; alpha <= 10.0; alpha += 1.0) {
        foil.SetOper(alpha, -1.0, 1.0e6, 0.0);
        if (!foil.Solve()) {
            continue;                                // skip failed points
        }
        const auto grad = foil.ComputeSensitivities();
        writer.WriteIteration(foil, &grad);          // or WriteIteration(foil)
    }
}
```

The file is closed when `writer` goes out of scope. Call `writer.Close()` to
close it earlier, for example before another process reads it.

## What's in the file

```text
alpha_sweep.h5
├── attributes: scirocco_version, created, num_iterations
├── iteration_000/
│   ├── operating_conditions   alpha, Re, Ma, Vinf, rho, viscous
│   ├── geometry               chord, hTE, source points
│   ├── panels/foil            x [2, N], s [N], t [2, N]
│   ├── panels/wake            x [2, Nw], s [Nw], t [2, Nw]
│   ├── coefficients           cl, cd, cdf, cdp, cdpi, cm, converged
│   ├── distributions          cp, cf, ue, th, ds, ...
│   └── sensitivities          only when you pass a Sensitivity
├── iteration_001/
└── ...
```

## Read it in Python

```python
import h5py

with h5py.File("alpha_sweep.h5", "r") as f:
    keys = sorted((k for k in f if k.startswith("iteration_")),
                  key=lambda k: int(k.split("_")[1]))
    for key in keys:
        g = f[key]
        alpha = g["operating_conditions/alpha"][()]
        cl = g["coefficients/cl"][()]
        cp = g["distributions/cp"][:]
        print(f"alpha = {alpha:5.1f}  Cl = {cl:.4f}  ({cp.size} Cp values)")
```

Sort the groups numerically, as above. Names are padded to three digits, so
plain string order breaks after `iteration_999`.

A complete version is
[`scirocco/examples/hdf5_output_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/hdf5_output_example.cpp).
For servers driven from Python, see [Record an HDF5 history](../python/hdf5-history).

<Next items={[
  { to: "/docs/api/hdf5-writer", title: "Hdf5Writer reference", body: "Every method and the full file layout." },
]} />
