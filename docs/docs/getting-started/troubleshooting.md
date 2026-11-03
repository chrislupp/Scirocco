---
sidebar_position: 4
---

# Troubleshooting

First reproduce the [quickstart](./quickstart) with its specified geometry,
panel count, operating conditions, and Philote versions. This separates setup
problems from failures introduced by a new shape or operating point.

## Build and connection failures

| Symptom | Action |
| --- | --- |
| CMake cannot find Eigen, Protobuf, gRPC, or HDF5 | Install the development packages in [Installation](./installation); check `CMAKE_PREFIX_PATH`. |
| CMake cannot find PhiloteCpp | Build the pinned Philote-Cpp revision and point `PhiloteCpp_DIR` at its build directory. |
| Unknown component options, deserialization errors, or invalid UTF-8 | Check the pinned C++ revision and `philote-mdo==0.8.0`, then rebuild and restart the servers. |
| `UNAVAILABLE` or connection refused | Check the server log, address, port, and process; start the binary by hand to see startup errors. |
| Address already in use | On macOS/Linux, `lsof -i :50051` identifies listeners. Stop the process you own or choose a free port in both server and client. |
| HDF5 history is absent | Rebuild with `SCIROCCO_HDF5_OUTPUT=ON`, restart the analysis server, and check its working directory. |

`0.0.0.0` is a server bind address for all IPv4 interfaces. Clients connect to
`localhost` or a reachable hostname/IP, not to the wildcard address. For Docker,
also check that the container port is published.

## Handle native solver failures

`Solve()` returns `false` for nonconvergence, including failure to reach a lift
target. Invalid inputs and nonfinite results can throw. A failed solve must not
be treated as a valid force evaluation or a source of optimizer gradients.

```cpp
#include <exception>
#include <iostream>
#include <scirocco/scirocco.h>

int main() {
    try {
        scirocco::Scirocco foil("0012", 100);
        foil.GetState().param.verb = 1;
        foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
        if (!foil.Solve()) {
            std::cerr << "Analysis did not converge\n";
            return 1;
        }
        const auto sens = foil.ComputeSensitivities();
        std::cout << "Cl = " << foil.GetCl()
                  << ", dCl/dalpha = " << sens.dcl_dalpha << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Analysis failed: " << error.what() << '\n';
        return 1;
    }
}
```

Set verbosity before solving. Calling the mutable `GetState()` accessor after
computing sensitivities invalidates their cache, even if you only wanted to
inspect a field. Use a const solver reference for read-only inspection.

## Handle remote failures

The analysis discipline raises an error when its solve fails; it does not
return a convergence output alongside `cl`, `cd`, and `cm`. Let a failed
`problem.run_model()` stop the current evaluation. Preserve the exception and
server log, and record the shape and operating conditions needed to reproduce
it. Do not substitute zeros for failed results or gradients.

In a sweep, a caller may record a failed case and continue to the next case,
but should read outputs only after a successful `run_model()`. Optimizer
recovery depends on the client and driver; there is no documented guarantee
that an arbitrary remote error will be treated as a recoverable trial point.

## Diagnose convergence

1. Check [input conventions](../api/conventions): finite values, node order,
   unit chord, positive thickness, no adjacent duplicates, and correct layout.
2. Try the same shape at a less demanding angle of attack. A failure limited
   to a narrow operating range may reflect separation or a transition change.
3. Check panel clustering and repeat with a different panel count. Convergence
   alone does not establish that a result is insensitive to discretization.
4. In C++, enable `param.verb = 1` and inspect the residual history. The
   [solver parameters](../concepts/solver-params) and [PTC guide](../concepts/ptc)
   describe iteration limits and damping. These fields are not options on the
   current remote discipline.
5. If the residual decreases but reaches the iteration limit, consider a larger
   `niglob`. If it oscillates or grows, investigate geometry and damping first.
   Relaxing tolerance changes the accuracy criterion and can degrade gradients.

The solver already enables continuation retries by default. More iterations
cannot make an invalid shape or an unsuitable physical model reliable.

## Diagnose gradients

Run the primal first at exactly the differentiation inputs. `ComputePartials`
requires a preceding `Compute` at those same inputs. Native cached getters
require a successful `ComputeSensitivities()` after the latest state change.

For a finite-difference check, choose a converged, repeatable operating point
away from transition changes and compare several suitably scaled steps. Judge
zero derivatives by absolute error. Do not perturb `Re=0` into positive values
to check an inviscid derivative: that changes the solve mode.

See the [sensitivity tutorial](../tutorials/sensitivities) for a runnable check.
When reporting a failure, include the revision, compiler/client versions,
build options, coordinates, operating point, and diagnostic log.
