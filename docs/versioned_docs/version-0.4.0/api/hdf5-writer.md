---
sidebar_position: 6
---

# `Hdf5Writer`

Optional HDF5 iteration output. The writer lives in a separate
`scirocco_hdf5` CMake target gated behind `-DSCIROCCO_HDF5_OUTPUT=ON`,
which also pulls in the HDF5 C++ dependency. Header:
[`scirocco/io/hdf5_writer.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/io/hdf5_writer.h).

```cpp
#include <scirocco/io/hdf5_writer.h>
```

## Interface

```cpp
class Hdf5Writer {
public:
    explicit Hdf5Writer(const std::string& filepath);
    ~Hdf5Writer();

    Hdf5Writer(const Hdf5Writer&)            = delete;
    Hdf5Writer& operator=(const Hdf5Writer&) = delete;

    void WriteIteration(const Scirocco&  solver,
                        const Sensitivity* sensitivity = nullptr);

    void WriteIteration(const FoilState&   state,
                        const Sensitivity* sensitivity = nullptr);

    void Close();
    [[nodiscard]] int NumIterations() const;
};
```

| Method | Use it when |
| --- | --- |
| `Hdf5Writer(path)` | Opening / creating the file. Truncates if it exists. |
| `WriteIteration(solver, sens?)` | Most common — pass the `Scirocco` after `Solve()`. |
| `WriteIteration(state, sens?)` | Lower-level — useful when driving the solver below the wrapper. |
| `Close()` | Optional — destructor closes too. Call explicitly to flush. |
| `NumIterations()` | Number of `WriteIteration` calls so far. |

The writer is neither copyable nor movable. Pass it by reference, or keep it
in a `std::unique_ptr` when ownership must be transferred.

## File layout

Every call to `WriteIteration` appends a new top-level group
`/iteration_<NNN>/`:

```
file.h5
├── attributes: scirocco_version, created, num_iterations
├── /iteration_000/
│   ├── operating_conditions # alpha, Re, Ma, Vinf, rho, viscous
│   ├── geometry             # chord, hTE, source points
│   ├── panels/foil          # x[2, N], s[N], t[2, N]
│   ├── panels/wake          # x[2, Nw], s[Nw], t[2, Nw]
│   ├── coefficients         # cl, cd, cdf, cdp, cdpi, cm, converged
│   ├── distributions        # cp, cf, ue, th, ds, ...
│   └── sensitivities        # optional shape and scalar derivatives
└── /iteration_001/ …
```

Iteration indices are 0-padded and monotonically increasing. The
optional `sensitivities` group is written only when you pass a non-null
`Sensitivity*`.

## Reading the file

Python:

```python
import h5py
with h5py.File("alpha_sweep.h5", "r") as f:
    for key in sorted(k for k in f if k.startswith("iteration_")):
        g     = f[key]
        alpha = g["operating_conditions/alpha"][()]
        cl    = g["coefficients/cl"][()]
        cp    = g["distributions/cp"][:]
        x     = g["panels/foil/x"][:]   # shape (2, N)
        print(f"α={alpha:5.1f}  Cl={cl:7.4f}")
```

The HDF5 datatypes are the natural ones for Eigen: scalar doubles,
1-D `VectorXd` as `dtype=float64` arrays, 2-D `MatrixXd` as `(rows, cols)`.

## Typical usage

```cpp
#include <stdexcept>

scirocco::Scirocco foil("2412", 199);
scirocco::Hdf5Writer writer("alpha_sweep.h5");

for (double a = -2.0; a <= 10.0; a += 1.0) {
    foil.SetOper(a, -1.0, 1.0e6, 0.0);
    if (!foil.Solve()) {
        throw std::runtime_error("Airfoil analysis did not converge");
    }

    auto sens = foil.ComputeSensitivities();
    writer.WriteIteration(foil, &sens);
}
```

See also [Save results to HDF5](../guides/hdf5) and the example
[`scirocco/examples/hdf5_output_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/hdf5_output_example.cpp).

## Building with HDF5 enabled

```bash
cmake -S scirocco -B scirocco/build-hdf5 -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-hdf5 -j
```

Consumers of the library link against the `scirocco_hdf5` target instead
of (or in addition to) `scirocco`.
