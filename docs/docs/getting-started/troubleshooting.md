---
title: Troubleshooting
description: Fixes for build errors, server connections, convergence failures, and suspicious gradients.
---

# Troubleshooting

Find your symptom below. If it isn't here, run the
[quickstart](./quickstart) first: if that works, the problem is in your
geometry or operating point, not your setup.

## Building

| Symptom | Fix |
| --- | --- |
| CMake can't find Eigen, gRPC, Protobuf, or HDF5 | Install the development packages from [Installation](./installation). On macOS add `-DCMAKE_PREFIX_PATH="$(brew --prefix)"`. |
| CMake can't find PhiloteCpp | Build Philote-Cpp at the pinned commit and pass `-DPhiloteCpp_DIR=<path>/Philote-Cpp/build`. |
| Tests you expected aren't running | Reference suites need HDF5 and the test suites need Google Test. Check with `ctest -N`. |
| Odd CMake behaviour after changing options | Delete the build directory and configure again. |

## Connecting from Python

| Symptom | Fix |
| --- | --- |
| `UNAVAILABLE` or connection refused | Start both servers, and check the ports match your script (50051 and 50052 by default). |
| Address already in use | Find the old process with `lsof -i :50051`, stop it, or pick another port on both sides. |
| `Option 'n_cst_upper' ... has not been declared` | The server is linked against an old Philote-Cpp. Rebuild with the pinned commit. |
| `Exception deserializing response` or invalid UTF-8 | Client and server use different Philote versions. Use `philote-mdo==0.8.0` with the pinned Philote-Cpp. |
| HDF5 history file is missing | Rebuild the server with `SCIROCCO_HDF5_OUTPUT=ON`. The file is written in the server's working directory. |

Servers bind to an address such as `0.0.0.0:50051`, but clients connect to
`localhost` or the host's name, never to `0.0.0.0`.

## The solve doesn't converge

`Solve()` returns `false` when the viscous Newton iteration fails or a lift
target can't be met. Don't use the forces or gradients from a failed solve.
Through Philote, a failed solve comes back as an error from the server; in
the OpenMDAO examples it raises an exception from `run_model()`.

Work through these in order:

1. **Check the geometry.** Plot it with equal axes. Nodes must run from the
   lower trailing edge, around the leading edge, to the upper trailing edge,
   with no repeated points. See [Conventions](../api/conventions).
2. **Open a closed trailing edge.** A viscous solve on a closed trailing
   edge (thinner than 1e-6 chord) throws `std::invalid_argument` instead of
   running, because the wake model can't handle it. Give the trailing edge a
   small gap, such as the 0.25% of chord that the built-in NACA generator
   uses, or solve inviscid. This limitation comes from mfoil; see
   [why a closed trailing edge fails](../guides/geometry#why-a-closed-trailing-edge-fails).
3. **Back off the angle of attack.** If a lower angle converges, you're
   likely near stall, where the integral boundary layer is least reliable.
4. **Change the panel count.** Try 120 or 160 panels instead of 199.
5. **Watch the residual.** Set `GetState().param.verb = 1` before solving to
   print it each iteration. If it is still falling when the iteration limit
   hits, raise `param.niglob`. If it oscillates or grows, the problem is
   the geometry or the operating point, not the limit.

[Solver parameters](../api/solver-params) and
[pseudo-transient continuation](../concepts/ptc) describe the knobs.

## Gradients look wrong

- Gradients need a converged solve at exactly the inputs you're
  differentiating. Call `ComputeSensitivities()` after a successful
  `Solve()`.
- Calling the non-const `GetState()` throws away cached gradients. Use a
  `const` reference to inspect state after computing them.
- When you compare against finite differences, use central differences at a
  point away from stall, and try several step sizes. A step that moves
  transition to another panel gives a meaningless difference.
- Never finite-difference `Re` starting from `Re = 0`: any positive step
  switches from the inviscid to the viscous model.

## Reporting a bug

Open an [issue](https://github.com/chrislupp/Scirocco/issues) with the
commit, compiler, build options, the coordinates, the operating point, and
the output with `verb = 1`.
