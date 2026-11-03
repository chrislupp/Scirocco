---
sidebar_position: 1
---

# Installation

Choose the interface you need:

- **C++ applications:** build and install the core library. Eigen is the only
  required library dependency.
- **Python analysis and optimization:** build the Philote servers and install
  the Python client environment below.
- **Containers:** follow [Docker](./docker) to run the servers without a native
  server build. Python clients still need the tutorial environment.

Commands below run from the repository root unless stated otherwise.

## Prerequisites

| Tool | Version |
| --- | --- |
| CMake | 3.18+ for the core; 3.23+ for Philote-Cpp 0.8 |
| C++ compiler | C++17 capable; GCC 14 and Clang 18 are tested in CI |
| Eigen | 3.4+ |
| Python (tutorials) | 3.12, used by release validation |

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

## Install native dependencies

Ubuntu 24.04 is the CI platform. For the core library:

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git libeigen3-dev
```

For the servers, Python tutorials, HDF5 output, and reference tests, also install:

```bash
sudo apt-get install -y libgrpc++-dev libprotobuf-dev protobuf-compiler \
    protobuf-compiler-grpc libhdf5-dev libgtest-dev pkg-config python3.12-venv
```

On macOS, with the Xcode command-line tools and Homebrew installed:

```bash
brew install cmake eigen
# Additional dependencies for servers, tutorials, and tests:
brew install grpc protobuf hdf5 googletest pkg-config python@3.12
```

The macOS commands use the packages available in your Homebrew installation;
Ubuntu CI is the reference configuration. If CMake cannot locate Homebrew
packages, add `-DCMAKE_PREFIX_PATH="$(brew --prefix)"` to configuration commands.

## Build the core library

```bash
cmake -S scirocco -B scirocco/build-default -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build-default -j
```

For a C++ application, continue to [Consuming Scirocco](#consuming-scirocco-from-your-own-cmake-project).
The remaining server and Python steps are needed for the Python tutorials.

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

The explicit commit makes the C++ implementation reproducible. Its submodule
points at the matching Python 0.8
message schema. The tutorial consumes Philote-Cpp's build-tree CMake package,
so no system-wide installation or `sudo` step is needed.

## Configure and build the servers

Build the core library, both servers, and the example programs:

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
python3.12 -m venv .venv
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
| `SCIROCCO_BUILD_TESTS` | `ON` | Solver, geometry, Jacobian, and adjoint tests; reference suites require HDF5 |
| `SCIROCCO_BUILD_EXAMPLES` | `ON` | `naca_example`, `viscous_example` |
| `SCIROCCO_HDF5_OUTPUT` | `OFF` | Builds the `scirocco_hdf5` library and the `hdf5_output_example` |
| `SCIROCCO_BUILD_PHILOTE` | `OFF` | Philote-Cpp gRPC discipline server (requires Protobuf, gRPC, PhiloteCpp) |
| `SCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS` | `OFF` | Chain-piece sensitivity tests for adjoint debugging |
| `SCIROCCO_BUILD_PROBES` | `OFF` | 14 `probe_*` diagnostic binaries (no gtest, no `ctest`) |
| `SCIROCCO_SANITIZERS` | `OFF` | AddressSanitizer and UndefinedBehaviorSanitizer with GCC or Clang |
| `SCIROCCO_NATIVE_ARCH` | `ON` | `-march=native` on Release builds |

Example — enable HDF5 output so the `hdf5_output_example` is built:

```bash
cmake -S scirocco -B scirocco/build-hdf5 -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build-hdf5 -j
```

## Run the tests

The suite covers the solver, geometry, Jacobians, and adjoints. Philote builds
also register discipline and RPC tests. Counts and runtime depend on features
and dependencies; inspect the actual suite before running it:

```bash
ctest --test-dir scirocco/build-philote -N
ctest --test-dir scirocco/build-philote --output-on-failure
```

Without HDF5, reference suites are omitted. Without Google Test, only the
installed-consumer check is registered. A passing small suite does not imply
that the reference or RPC tests ran.

## Consuming Scirocco from your own CMake project

Install the core library into a prefix you control:

```bash
cmake --install scirocco/build-default --prefix "$PWD/install"
```

In your application's `CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.18)
project(my_app LANGUAGES CXX)
find_package(scirocco CONFIG REQUIRED)
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE scirocco::scirocco)
```

Use the [minimal C++ example](../api/scirocco#minimal-example) as `main.cpp`,
then configure your application with the absolute install prefix:

```bash
cmake -S /path/to/my_app -B /path/to/my_app/build \
    -DCMAKE_PREFIX_PATH=/absolute/path/to/Scirocco/install
cmake --build /path/to/my_app/build
```

Eigen must also be discoverable by CMake. Alternatively, include the source
as an in-tree dependency:

```cmake
add_subdirectory(path/to/Scirocco/scirocco scirocco-build)
target_link_libraries(my_app PRIVATE scirocco)
```

The installed package exports the core library. The optional HDF5 writer is
currently consumed from a source build using the `scirocco_hdf5` target.

## CI and validation

The repository tests GCC 14 and Clang 18 Release builds on Ubuntu 24.04,
including Philote and HDF5 paths. Additional jobs cover sensitivity diagnostics,
installed consumers, sanitizers, and static analysis. The workflow definitions
in `.github/workflows/` specify the enabled jobs and their triggers.

Coverage CI is disabled. For a local measurement, follow the coverage commands
in the repository README; historical coverage percentages are not a guarantee
for the current revision.

## Next steps

Run the [Python quickstart](./quickstart), or use the
[C++ API reference](../api/scirocco). If setup or analysis fails, see
[Troubleshooting](./troubleshooting).
