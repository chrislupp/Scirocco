---
sidebar_position: 5
---

# HDF5 output

Scirocco can persist every solved operating point — state, distributions,
and optionally sensitivities — to a single HDF5 file. The
[`Hdf5Writer`](../api/hdf5-writer) is provided by an opt-in extension
library that's not built by default. This tutorial mirrors
[`scirocco/examples/hdf5_output_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/hdf5_output_example.cpp).

## Build with HDF5 enabled

The output writer is gated behind a CMake option that pulls in HDF5 and an
extra `scirocco_hdf5` library:

```bash
cmake -S scirocco -B scirocco/build-hdf5 -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-hdf5 -j
```

This also builds the `hdf5_output_example` executable.

## The program

```cpp
#include <cstdio>

#include "scirocco/io/hdf5_writer.h"
#include "scirocco/scirocco.h"

int main() {
    scirocco::Scirocco  foil("2412", 199);
    scirocco::Hdf5Writer writer("alpha_sweep.h5");

    const double Re = 1.0e6;
    const double Ma = 0.0;

    for (double alpha = -2.0; alpha <= 10.0; alpha += 1.0) {
        foil.SetOper(alpha, -1.0, Re, Ma);
        bool converged = foil.Solve();

        auto sens = foil.ComputeSensitivities();
        writer.WriteIteration(foil, &sens);

        std::printf("alpha=%5.1f  cl=%7.4f  cd=%8.5f  cm=%7.4f  %s\n",
                    alpha, foil.GetCl(), foil.GetCd(), foil.GetCm(),
                    converged ? "OK" : "NOT CONVERGED");
    }

    std::printf("Wrote %d iterations to alpha_sweep.h5\n",
                writer.NumIterations());
}
```

## The API

```cpp
class Hdf5Writer {
public:
    explicit Hdf5Writer(const std::string& filepath);
    ~Hdf5Writer();                              // closes the file

    void WriteIteration(const Scirocco& solver,
                        const Sensitivity* sensitivity = nullptr);

    void WriteIteration(const FoilState& state,
                        const Sensitivity* sensitivity = nullptr);

    void Close();
    int  NumIterations() const;
};
```

| Method | Use it when |
| --- | --- |
| `Hdf5Writer(path)` | Opening / creating the output file. Truncates if it exists. |
| `WriteIteration(solver, sens?)` | Most common — pass the `Scirocco` object after a `Solve()`. |
| `WriteIteration(state, sens?)` | Lower-level — pass a raw `FoilState`. Use when you've been driving the solver below the wrapper. |
| `Close()` | Optional; the destructor closes too. Call it if you need the file flushed before the writer goes out of scope. |
| `NumIterations()` | How many `WriteIteration` calls have been made so far. |

The writer is move-only (`= delete`'d copy semantics) — pass it by reference
or by `std::move`.

## File layout

Each call to `WriteIteration` appends a new group under `/iterations/`:

```
alpha_sweep.h5
├── /iterations/
│   ├── 000/
│   │   ├── oper           # alpha, cl_target, Re, Ma
│   │   ├── post           # cl, cd, cdf, cdp, cm, cp[…], cf[…], theta[…], ds[…]
│   │   ├── foil           # x[2,N]
│   │   ├── wake           # x[2,Nw]
│   │   └── sensitivity    # (optional) dcl_dx[2,N], dcd_dx[2,N], …
│   ├── 001/
│   └── …
└── /metadata               # solver version, panel count
```

In Python:

```python
import h5py

with h5py.File("alpha_sweep.h5", "r") as f:
    for key in sorted(f["iterations"]):
        g = f[f"iterations/{key}"]
        alpha = g["oper/alpha"][()]
        cl    = g["post/cl"][()]
        cd    = g["post/cd"][()]
        print(f"{alpha:5.1f}  {cl:7.4f}  {cd:8.5f}")
```

## Passing sensitivities

The `sensitivity` argument is optional. Pass it when you want the gradient
information persisted alongside the primal — useful for offline gradient
postprocessing, finite-difference verification, or feeding into an
optimization audit trail:

```cpp
auto sens = foil.ComputeSensitivities();
writer.WriteIteration(foil, &sens);
```

Pass `nullptr` (the default) to skip the sensitivity group entirely. See
the [Sensitivities tutorial](./sensitivities) for what
`ComputeSensitivities` populates.

## Full source

[`scirocco/examples/hdf5_output_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/hdf5_output_example.cpp).

## Next: [Adjoint sensitivities](./sensitivities) — the gradients you just persisted, explained.
