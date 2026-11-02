---
sidebar_position: 1
---

# Installation

Scirocco is a CMake-based C++17 library. The tutorials use its Philote-Cpp
servers through the Philote-MDO 0.8.0 Python bindings, so install both the
native server dependencies and the Python client environment.

## Prerequisites

| Tool | Version |
| --- | --- |
| CMake | 3.23+ (required by Philote-Cpp 0.8) |
| C++ compiler | C++17 capable; GCC 14 and Clang 18 are tested in CI |
| Eigen | 3.4+ |

Optional, depending on the features you want:

| Component | Needed for |
| --- | --- |
| HDF5 (CXX) | `SCIROCCO_HDF5_OUTPUT`, reference validation tests |
| Google Test | The unit and reference test suites |
| Protobuf, gRPC, [Philote-Cpp](https://github.com/MDO-Standards/Philote-Cpp) | The MDO discipline server (`SCIROCCO_BUILD_PHILOTE`, off by default) |
| LLVM tools (`llvm-cov`, `llvm-profdata`) | Coverage measurement |

## Clone

```bash
git clone https://github.com/chrislupp/Scirocco.git
cd Scirocco
```

## Install Philote-Cpp

Philote-Cpp provides the discipline and gRPC server infrastructure. Use the
tested 0.8 commit below; it pins the matching `v0.8.0` protocol submodule and
is interoperable with `philote-mdo==0.8.0`.

```bash
git clone --recurse-submodules https://github.com/MDO-Standards/Philote-Cpp.git
git -C Philote-Cpp checkout cd19666a23ddca009ee8b27ecdd40c0a0a2ea4b4
git -C Philote-Cpp submodule update --init --recursive
cmake -S Philote-Cpp -B Philote-Cpp/build -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF
cmake --build Philote-Cpp/build -j
```

The explicit commit keeps the C++ implementation reproducible until the 0.8
release is tagged. Its submodule already points at the exact Python 0.8
message schema. The tutorial consumes Philote-Cpp's build-tree CMake package,
so no system-wide installation or `sudo` step is needed.

## Configure and build the servers

A lean Release build that produces the library and the example programs:

```bash
cmake -S scirocco -B scirocco/build-philote -DCMAKE_BUILD_TYPE=Release \
      -DPhiloteCpp_DIR="$PWD/Philote-Cpp/build" \
      -DSCIROCCO_BUILD_PHILOTE=ON
cmake --build scirocco/build-philote -j
```

This produces the two binaries used throughout the tutorials:

```text
scirocco/build-philote/src/philote/scirocco_server
scirocco/build-philote/src/philote/scirocco_kulfan_server
```

## Install the Python bindings

Use a virtual environment so OpenMDAO and the Philote client do not affect
your system Python:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r examples/tutorials/requirements.txt
```

There is no separate `scirocco` Python wheel. Python communicates with the
compiled Scirocco servers through `philote-mdo==0.8.0` and gRPC. Confirm the
critical client version with:

```bash
python -c "from importlib.metadata import version; print(version('philote-mdo'))"
```

The command must print `0.8.0`.

Local build directories are git-ignored and disposable; if CMake or test
discovery state ever gets stale, delete the build dir and re-configure.

### Build options

| Option | Default | Purpose |
| --- | --- | --- |
| `SCIROCCO_BUILD_TESTS` | `ON` | Default test suite (five gtest executables — primal + Jacobians + end-to-end adjoint vs. complex-step) |
| `SCIROCCO_BUILD_EXAMPLES` | `ON` | `naca_example`, `viscous_example` |
| `SCIROCCO_HDF5_OUTPUT` | `OFF` | Builds the `scirocco_hdf5` library and the `hdf5_output_example` |
| `SCIROCCO_BUILD_PHILOTE` | `OFF` | Philote-Cpp gRPC discipline server (requires Protobuf, gRPC, PhiloteCpp) |
| `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS` | `OFF` | Chain-piece sensitivity tests for adjoint debugging |
| `SCIROCCO_BUILD_PROBES` | `OFF` | 14 `probe_*` diagnostic binaries (no gtest, no `ctest`) |
| `SCIROCCO_NATIVE_ARCH` | `ON` | `-march=native` on Release builds |

Example — enable HDF5 output so the `hdf5_output_example` is built:

```bash
cmake -S scirocco -B scirocco/build-hdf5 -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-hdf5 -j
```

## Run the tests

The default suite is fast (normally well under a minute) and covers the
primal solver, geometry, Jacobians, Philote disciplines, and end-to-end
adjoint checks. The discovered test count varies with optional HDF5 and
debug suites:

```bash
ctest --test-dir scirocco/build-philote --output-on-failure
```

If you don't have HDF5 installed, the reference test executables (which load
`.h5` fixtures) silently drop out of the build; the smaller spline/geometry
tests still run.

## Consuming Scirocco from your own CMake project

The currently supported downstream path is an in-tree CMake dependency:

```cmake
add_subdirectory(path/to/Scirocco/scirocco)
target_link_libraries(my_app PRIVATE scirocco)
```

The install tree exports targets but does not yet provide the package config
needed by `find_package(scirocco)`.

## CI matrix

Every push and PR runs against the following on Ubuntu 24.04:

| Compiler | Build | Suite |
| --- | --- | --- |
| GCC 14 | Release | Default tests |
| Clang 18 | Release | Default tests |
| GCC 14 | Release | Default + sensitivity-debug + probes |
| Clang 18 | Release | Default + LLVM coverage instrumentation (uploaded to Codecov) |
| Clang 18 | — | clang-tidy on PRs |

If your environment matches one of these, you have a known-good configuration.

## Next: [Python quickstart](./quickstart) — run your first analysis through Philote.
