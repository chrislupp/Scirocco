---
title: 6. Record a history
description: Have the analysis server write every evaluation to an HDF5 file.
---

# 6. Record a history

Give the analysis component an `hdf5_file` option and the server appends a
record for every evaluation: operating point, coefficients, distributions,
and coordinates.

```bash
python examples/tutorials/06_hdf5_history.py
# stop scirocco_server with Ctrl-C, then:
python examples/tutorials/read_hdf5_history.py analysis_history.h5
```

## Build the server with HDF5

```bash
cmake -S scirocco -B scirocco/build-philote -DCMAKE_BUILD_TYPE=Release \
      -DPhiloteCpp_DIR="$PWD/Philote-Cpp/build" \
      -DSCIROCCO_BUILD_PHILOTE=ON \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-philote -j
```

Restart `scirocco_server` after rebuilding.

## Turn it on

Add one option to the analysis component:

```python
problem.model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(
        channel=grpc.insecure_channel("localhost:50051"),
        naca="0012",
        npanel=npanel,
        hdf5_file="analysis_history.h5",
    ),
    promotes_inputs=["x"],
)

for alpha in np.linspace(-2.0, 8.0, 11):
    problem.set_val("analysis.alpha", alpha, units="deg")
    problem.set_val("analysis.Re", 1.0e6)
    problem.set_val("analysis.Ma", 0.0)
    problem.run_model()
```

The **server** writes the file, so the path is relative to the server's
working directory. That matters when the server runs in Docker or on another
machine.

Stop the server before reading the file, so everything is flushed.

## Read it

```python
import h5py

with h5py.File("analysis_history.h5", "r") as history:
    keys = sorted(
        (key for key in history if key.startswith("iteration_")),
        key=lambda key: int(key.rsplit("_", 1)[1]),
    )
    for key in keys:
        iteration = history[key]
        alpha = iteration["operating_conditions/alpha"][()]
        cl = iteration["coefficients/cl"][()]
        cd = iteration["coefficients/cd"][()]
        print(f"{alpha:6.2f}  {cl:9.5f}  {cd:9.6f}")
```

The layout is the same as for the [C++ writer](../guides/hdf5#whats-in-the-file).
The server records solves, not gradients; the C++ `Hdf5Writer` can also
store a `Sensitivity`.

:::warning Each new component truncates the file
Creating another analysis component with the same `hdf5_file` starts a new
file. Use a separate path for each run, and create the component once per
optimization.
:::

That covers everything the optimization examples use: parameterized
geometry, remote analysis, analytic totals, and a history.

<Next items={[
  { to: "/docs/optimization/naca-design", title: "NACA 4-digit optimization", body: "Minimize drag at fixed lift with three shape variables." },
]} />
