---
sidebar_position: 6
---

# 6. Record an HDF5 history

The final analysis tutorial adds persistent output to the same
Python/Philote model. HDF5 is configured when the remote analysis component is
created; every call to that component then appends an iteration.

## Rebuild with HDF5 support

The server must link the optional writer:

```bash
cmake -S scirocco -B scirocco/build-philote -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_BUILD_PHILOTE=ON \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-philote -j
```

Restart `scirocco_server` after rebuilding it.

## Add the component option

In the model from the previous tutorials, change only the analysis component
construction:

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
```

The file is opened by the C++ server process, so the path is relative to the
server's working directory—not necessarily the Python client's directory when
the server runs on another machine or in a container.

Now run a sweep exactly as in tutorial 1:

```python
for alpha in np.linspace(-2.0, 8.0, 11):
    problem.set_val("analysis.alpha", alpha, units="deg")
    problem.set_val("analysis.Re", 1.0e6)
    problem.set_val("analysis.Ma", 0.0)
    problem.run_model()
```

Stop the analysis server (or otherwise close its writer) before opening the
file in another process, so all buffers are flushed.

## Read the history in Python

```python
import h5py

with h5py.File("analysis_history.h5", "r") as history:
    for key in sorted(k for k in history if k.startswith("iteration_")):
        iteration = history[key]
        alpha = iteration["operating_conditions/alpha"][()]
        cl = iteration["coefficients/cl"][()]
        cd = iteration["coefficients/cd"][()]
        print(f"{alpha:6.2f}  {cl:9.5f}  {cd:9.6f}")
```

Each evaluation records operating conditions, force coefficients,
distributions, foil coordinates, and wake coordinates under a numbered group.
Through the Philote option, the server currently records the primal state. The
native C++ `Hdf5Writer` API can additionally receive a `Sensitivity` object;
see the [writer reference](../api/hdf5-writer) for that lower-level path.

## Avoid accidental truncation

Creating a new analysis component with the same `hdf5_file` constructs a new
writer and truncates the existing file. Use a unique path for each run. For an
optimization, construct one component once and let all driver evaluations
append to it.

You now have the pieces used by the optimization tutorials: parameterized
geometry, remote analysis, analytic total derivatives, and iteration history.

## Next: [NACA 4-digit optimization](../optimization/naca-design)
