---
sidebar_position: 2
---

# OpenMDAO setup

OpenMDAO does not call Scirocco directly — it talks to the [Philote gRPC
server](./server) through `philote_mdo.openmdao.RemoteExplicitComponent`,
which presents the discipline as a standard `ExplicitComponent` with
analytic partials. The pattern is the same for the NACA, CST, and
multipoint workflows; what changes is how you stitch the components
together.

## Required Python dependencies

```bash
pip install openmdao philote-mdo grpcio numpy matplotlib h5py
```

`philote-mdo` provides `RemoteExplicitComponent`. Scirocco's examples are
also tested against the versions pinned in [Philote-MDO's repo](https://github.com/chrislupp/Philote-MDO).

## Minimal connection

```python
import grpc, openmdao.api as om
from philote_mdo.openmdao import RemoteExplicitComponent

# Server already running at localhost:50051

channel = grpc.insecure_channel("localhost:50051")

prob = om.Problem()
prob.model.add_subsystem(
    "analysis",
    RemoteExplicitComponent(channel=channel, naca="0012", npanel=100),
)

prob.setup()
prob.set_val("analysis.alpha", 5.0)
prob.set_val("analysis.Re",    1.0e6)
prob.set_val("analysis.Ma",    0.0)

prob.run_model()

print("Cl =", prob.get_val("analysis.cl")[0])

channel.close()
```

The `RemoteExplicitComponent` constructor takes the gRPC channel plus any
discipline-specific options (`naca`, `npanel`, `hdf5_file`, …) and
exposes the discipline's variables on the OpenMDAO graph.

## Connecting geometry → analysis

For NACA-4-digit-parameterized optimization you build the geometry as a
local OpenMDAO component (`NACA4Geometry`) and connect its `x` output to
the analysis discipline's `x` input. For CST, the geometry is a **second**
remote discipline (`KulfanDiscipline`) and you connect the two remote
components together.

```python
model.add_subsystem("geometry",
    RemoteExplicitComponent(channel=geom_channel,
                            n_cst_upper=8, n_cst_lower=8, npanel=100),
    promotes_outputs=["x"],
    promotes_inputs=["cst_upper", "cst_lower"])

model.add_subsystem("analysis",
    RemoteExplicitComponent(channel=analysis_channel, naca="0012", npanel=100),
    promotes_inputs=["x"])
```

`promotes_outputs=["x"]` on geometry and `promotes_inputs=["x"]` on
analysis is the easiest way to make the data flow visible at the top
level. OpenMDAO will compute the chain rule through both disciplines'
analytic partials when SLSQP asks for total derivatives.

## Driving the optimizer

```python
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options["optimizer"] = "SLSQP"
prob.driver.options["tol"]       = 1e-6
prob.driver.options["maxiter"]   = 200
prob.driver.options["disp"]      = True

prob.model.add_design_var("cst_upper", lower=0.0, upper=0.4)
prob.model.add_design_var("cst_lower", lower=0.0, upper=0.4)
prob.model.add_design_var("analysis.alpha", lower=-2.0, upper=10.0)

prob.model.add_objective("analysis.cd", ref=0.01)
prob.model.add_constraint("analysis.cl", equals=0.5)
```

The `ref` arguments are scaling — they're applied to the optimizer's
internal representation but not to your printed results.

SLSQP is the default-recommended optimizer because it handles the typical
"minimize $C_d$ at a Cl-equality and area floor" structure cleanly with
small DV counts (≤ ~25). For larger CST counts or many points consider
SNOPT or IPOPT.

## Verifying derivatives

Before kicking off a long optimization, sanity-check the gradients
through both disciplines:

```python
prob.run_model()
totals = prob.check_totals(
    of=["analysis.cl", "analysis.cd", "analysis.cm"],
    wrt=["cst_upper", "cst_lower", "analysis.alpha"],
    compact_print=True,
)
```

`check_totals` runs finite-difference perturbations on the inputs and
compares to OpenMDAO's chained analytic Jacobian. Expect relative errors
around $10^{-5}$ – $10^{-7}$ for a smooth converged primal. Larger errors
usually mean either a bad initial guess (primal not converged at the
linearization point) or a configuration where one of the disciplines is
hitting a discontinuity (transition jumping panels, near-stall).

## Recording iteration history

The analysis discipline accepts an `hdf5_file` option that writes every
primal + sensitivity to disk through Scirocco's
[`Hdf5Writer`](../tutorials/hdf5-output):

```python
RemoteExplicitComponent(
    channel=analysis_channel,
    naca="0012", npanel=100,
    hdf5_file="opt_history.h5",
)
```

The file is appended once per primal call, so for a SLSQP run with $N$
total residual evaluations you get $N$ iterations. Useful for offline
plotting, gradient debugging, or post-hoc analysis of what the optimizer
visited.

## Cleaning up

End every run by closing the channel and terminating the server
processes:

```python
channel.close()
server_proc.terminate()
server_proc.wait()
```

The examples wrap this in a `try/finally` block so the servers don't
linger on `Ctrl+C`.

## Next

Pick a parameterization and go end-to-end:

- [NACA-4-digit optimization](./naca-design)
- [CST inviscid and viscous design](./cst-design)
- [Multipoint CST design](./multipoint)
