---
title: Quickstart
description: Build Scirocco and solve a viscous airfoil in C++.
---

# Quickstart

Build the library and run a viscous analysis of a NACA 2412. You need CMake
3.18+, a C++17 compiler, and Eigen 3.4+.

:::tip Working in Python instead?
Skip the build: run the servers from Docker images and follow
[Python and MDO frameworks](../python#set-up).
:::

## 1. Install the dependencies

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

## 2. Build

```bash
git clone https://github.com/chrislupp/Scirocco.git
cd Scirocco
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j --target quickstart
```

## 3. Run

The `quickstart` target is this program:

<QuickstartSource />

```bash
scirocco/build/quickstart
```

```output
Cl = 0.720156
Cd = 0.00693341
Cm = -0.0585026
```

That's a full viscous solve, with boundary layer, transition, and wake, in
about a quarter of a second. Your last digits may differ slightly with a
different compiler.

`SetOper(alpha, cl, Re, Ma)` sets the operating point:

| Argument | Meaning |
| --- | --- |
| `alpha` | Angle of attack in degrees |
| `cl` | Target lift coefficient. Pass `-1.0` to solve at the given `alpha` instead. |
| `Re` | Chord Reynolds number. `0` gives an inviscid solve. |
| `Ma` | Freestream Mach number |

## Use it in your own project

Install the library, then link it with CMake:

```bash
cmake --install scirocco/build --prefix "$PWD/install"
```

```cmake title="CMakeLists.txt"
cmake_minimum_required(VERSION 3.18)
project(my_app LANGUAGES CXX)
find_package(scirocco CONFIG REQUIRED)
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE scirocco::scirocco)
```

Point CMake at the install prefix when you configure:

```bash
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/Scirocco/install
```

If you installed Scirocco with [Homebrew](./installation#homebrew), skip the
install step and use `-DCMAKE_PREFIX_PATH="$(brew --prefix)"`.

[Installation](./installation) covers the other ways to build, including
the optional HDF5 and Philote components.

<Next items={[
  { to: "/docs/guides/analysis", title: "Run an analysis", body: "Sweeps, lift targets, Mach number, and surface distributions." },
  { to: "/docs/guides/gradients", title: "Compute gradients", body: "Derivatives of Cl, Cd, and Cm for optimization." },
  { to: "/docs/python", title: "Use it from Python", body: "Connect an MDO framework through Philote." },
]} />
