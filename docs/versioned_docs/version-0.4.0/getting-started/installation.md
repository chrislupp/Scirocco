---
title: Installation
description: Build Scirocco, with optional HDF5 output and Philote servers.
---

# Installation

:::tip Using Scirocco from Python or an MDO framework?
You don't have to build anything. Both Philote servers come as Docker
images:

<DockerRun />

Then continue with [Python and MDO frameworks](../python). Build from source
when you want the C++ library or HDF5 history files. See [Docker](../python/docker) for the image tags.
:::

The core library needs only Eigen. Everything else is optional: turn on what
you use. All commands run from the repository root.

| You want | You need | Build option |
| --- | --- | --- |
| The C++ library | CMake 3.18+, a C++17 compiler, Eigen 3.4+ | (default) |
| The C++ library, without building | Homebrew, on macOS or Linux | none: use the [formula](#homebrew) |
| HDF5 history files | HDF5 with C++ bindings | `SCIROCCO_HDF5_OUTPUT=ON` |
| Philote servers, for Python and MDO frameworks | Protobuf, gRPC, Philote-Cpp 0.6.0 | `SCIROCCO_BUILD_PHILOTE=ON` |
| Philote servers, without building | Docker | none: use the [images](../python/docker) |
| The test suites | Google Test, and HDF5 for the reference tests | (on by default) |

CI builds with GCC 14 and Clang 18 on Ubuntu 24.04. macOS with Homebrew works
too.

## Homebrew

On macOS or Linux, Homebrew installs the C++ library from
[chrislupp/tools](https://github.com/chrislupp/homebrew-tools), a tap that
follows each release:

```bash
brew install chrislupp/tools/scirocco
```

This builds the core library from the release's source and installs its
headers and CMake package. On macOS the dense solves use Accelerate; on
Linux, OpenBLAS. CMake finds the package under Homebrew's prefix, so link it
as shown in [the quickstart](./quickstart#use-it-in-your-own-project) and add
`-DCMAKE_PREFIX_PATH="$(brew --prefix)"` if CMake can't find it.

The formula has the core library only. Build from source for HDF5 history
files, the Philote servers, the tests or the examples. Upgrade with
`brew upgrade scirocco`.

[fauXFOIL](../guides/fauxfoil), the XFOIL-style command line, has its own
formula, so library users don't get it or Qt:

```bash
brew install chrislupp/tools/fauxfoil
```

## The C++ library from source

<Tabs groupId="os">
<TabItem value="mac" label="macOS">

```bash
brew install cmake eigen
```

</TabItem>
<TabItem value="ubuntu" label="Ubuntu">

```bash
sudo apt-get install -y build-essential cmake git libeigen3-dev
```

</TabItem>
</Tabs>

```bash
git clone https://github.com/chrislupp/Scirocco.git
cd Scirocco
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j
```

This builds the library, the tests (if Google Test is installed), and the
example programs. To use Scirocco from your own CMake project, see
[the quickstart](./quickstart#use-it-in-your-own-project). You can also add it
as a subdirectory:

```cmake
add_subdirectory(path/to/Scirocco/scirocco scirocco-build)
target_link_libraries(my_app PRIVATE scirocco)
```

:::tip If CMake can't find Homebrew packages
Add `-DCMAKE_PREFIX_PATH="$(brew --prefix)"` to the configure command.
:::

## HDF5 output

<Tabs groupId="os">
<TabItem value="mac" label="macOS">

```bash
brew install hdf5
```

</TabItem>
<TabItem value="ubuntu" label="Ubuntu">

```bash
sudo apt-get install -y libhdf5-dev
```

</TabItem>
</Tabs>

```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release \
      -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build scirocco/build -j
```

Link against the `scirocco_hdf5` target. The writer is available from a source
build only; the installed CMake package exports the core library. See
[Save results to HDF5](../guides/hdf5).

## Philote servers for Python

Python talks to Scirocco through two gRPC servers built on
[Philote-Cpp](https://github.com/MDO-Standards/Philote-Cpp). Install gRPC and
Protobuf:

<Tabs groupId="os">
<TabItem value="mac" label="macOS">

```bash
brew install grpc protobuf pkg-config
```

</TabItem>
<TabItem value="ubuntu" label="Ubuntu">

```bash
sudo apt-get install -y libgrpc++-dev libprotobuf-dev protobuf-compiler \
    protobuf-compiler-grpc pkg-config
```

</TabItem>
</Tabs>

Build Philote-Cpp v0.6.0, the tested release. It matches `philote-mdo==0.8.0` on
the Python side; other versions use a different wire format.

```bash
git clone --recurse-submodules https://github.com/MDO-Standards/Philote-Cpp.git
git -C Philote-Cpp checkout 49e1b42c3e46a4b9b70619c29e705985a6e6555d
git -C Philote-Cpp submodule update --init --recursive
cmake -S Philote-Cpp -B Philote-Cpp/build -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=OFF
cmake --build Philote-Cpp/build -j
```

Then build Scirocco with the servers. No install step is needed; Scirocco
uses Philote-Cpp's build tree directly.

```bash
cmake -S scirocco -B scirocco/build-philote -DCMAKE_BUILD_TYPE=Release \
      -DPhiloteCpp_DIR="$PWD/Philote-Cpp/build" \
      -DSCIROCCO_BUILD_PHILOTE=ON
cmake --build scirocco/build-philote -j
```

This produces `scirocco_server` and `scirocco_kulfan_server` in
`scirocco/build-philote/src/philote/`. Continue with
[Python and MDO frameworks](../python) to set up a client. If you'd rather not
build the servers, use [Docker](../python/docker).

## Run the tests

```bash
ctest --test-dir scirocco/build --output-on-failure
```

Which suites run depends on what you installed. Without HDF5 the reference
suites are skipped; without Google Test only the install check runs. List
what's registered with `ctest --test-dir scirocco/build -N`.

All build options are listed under [Build options](../api/build-options).
If something goes wrong, see [Troubleshooting](./troubleshooting).
