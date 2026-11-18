---
sidebar_position: 3
---

# `FoilState`

Everything the solver knows. [`Scirocco`](./scirocco) holds one and
exposes the common fields through getters; read it directly for the rest.
Header
[`scirocco/data/state.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/state.h).

```cpp
struct FoilState {
    std::string version = "1.0.0";
    Geom   geom;     // input geometry (xpoint, NACA digits)
    Panel  foil;     // airfoil panel discretization
    Panel  wake;     // wake panel discretization
    Oper   oper;     // operating conditions
    Isol   isol;     // inviscid linear-system data
    Vsol   vsol;     // viscous-system bookkeeping
    Glob   glob;     // global Newton state and Jacobian
    Post   post;     // forces and distributions
    Param  param;    // solver parameters
};
```

Access via:

```cpp
FoilState& state = foil.GetState();
state.param.verb = 2;
double cl = state.post.cl;
```

## Sub-structs in detail

### `Panel`

The panel discretization. Used for both `foil` and `wake`.
[`scirocco/data/panel.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/panel.h).

| Field | Shape | Meaning |
| --- | --- | --- |
| `N`  | int | Number of nodes |
| `x`  | `2 × N` | Node coordinates |
| `s`  | `N` | Arc-length at each node |
| `t`  | `2 × N` | Unit tangent $(dx/ds, dz/ds)$ at each node |

### `Oper`

Operating conditions, mostly set by `SetOper`.
[`scirocco/data/oper.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/oper.h).

| Field | Default | Meaning |
| --- | --- | --- |
| `Vinf` | `1.0` | Freestream speed magnitude |
| `alpha` | `0.0` | Angle of attack, degrees |
| `rho` | `1.0` | Density |
| `cltgt` | `0.0` | Target lift coefficient (active when `givencl = true`) |
| `givencl` | `false` | True if the solver should iterate $\alpha$ to hit `cltgt` |
| `initbl` | `true` | Initialize the boundary layer at the next solve |
| `viscous` | `false` | Viscous (BL-coupled) analysis |
| `redowake` | `false` | Rebuild wake after an $\alpha$ change |
| `Re` | `1e5` | Reynolds number |
| `Ma` | `0.0` | Mach number |

### `Post`

Output quantities, populated by the solver and `CalcForce`.
[`scirocco/data/post.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/post.h).

Scalars:

| Field | Meaning |
| --- | --- |
| `cl`, `cd`, `cm` | Force coefficients |
| `cdf`, `cdp` | Friction / pressure drag components |
| `cdpi` | Near-field inviscid pressure drag |
| `cl_alpha`, `cdpi_alpha` | Linearizations w.r.t. $\alpha$ |

Distributions (`Eigen::VectorXd`):

| Field | Meaning |
| --- | --- |
| `cp`, `cpi` | Compressible and inviscid $C_p$ |
| `cf` | Skin friction $C_f$ |
| `th`, `ds`, `sa` | $\theta$, $\delta^*$, amplification or $C_\tau$ |
| `ue`, `uei` | Edge velocity, compressible and inviscid |
| `Ret`, `Hk` | $Re_\theta$ and kinematic shape parameter $H_k$ |

`Post` also holds partial Jacobians for the adjoint, such as
`dcl_dx_partial`, `cl_V`, and `cd_V`. They're internal.

### `Param`

Solver parameters.
[`scirocco/data/param.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/param.h).

### `Geom`, `Isol`, `Vsol`, `Glob`

Internal bookkeeping structs:

| Struct | Role | Header |
| --- | --- | --- |
| `Geom` | Raw input point cloud (`xpoint`) plus NACA-digit string | [`data/geom.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/geom.h) |
| `Isol` | Inviscid AIC matrix, $\gamma$ solutions, edge-velocity sensitivities | [`data/isol.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/isol.h) |
| `Vsol` | Viscous-station bookkeeping (transition indices, station-type masks) | [`data/vsol.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/vsol.h) |
| `Glob` | Global Newton state vector $U$, sparse Jacobian $R_U$, residual, convergence flag | [`data/glob.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/glob.h) |

These belong to the solver kernels and tests.

## When to reach below the wrapper

Common patterns:

```cpp
// Solver verbosity
foil.GetState().param.verb = 1;

// Non-default transition threshold
foil.GetState().param.ncrit = 11.0;

// Turn off PTC for an A/B comparison
foil.GetState().param.use_ptc = false;

// After solving, inspect through a const reference to preserve cached gradients.
const scirocco::Scirocco& solved = foil;
const Eigen::VectorXd& Hk = solved.GetState().post.Hk;

// Inspect the stored BL residual (not the complete coupled residual).
double rnorm = solved.GetState().glob.R.norm();
```

Writing to `glob.U` or the Jacobian directly is outside the supported API;
the unit tests show how it's done.

### Watch the Newton iteration

`glob.iteration_observer` is an optional callback that the coupled viscous
solve calls at every Newton iteration, after the residual and forces are
evaluated and before the convergence test. It receives the state and the
0-based iteration number. From iteration 1 on, `glob.dU` holds the step
that produced the current `glob.U` and `glob.relax` the under-relaxation
factor applied to it. The callback reads the state; it must not modify it,
and it does not change the solve. [`fauXfoil`](../guides/fauxfoil.md) uses it to
print XFOIL's iteration lines.

```cpp
scirocco::Scirocco foil("2412", 160);
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.GetState().glob.iteration_observer =
    [](const scirocco::FoilState& s, int iteration) {
        std::printf("%3d  CL = %.4f  CD = %.5f\n", iteration, s.post.cl,
                    s.post.cd);
    };
foil.Solve();
```

`Glob::relax` and `Glob::iteration_observer` are declared in
[`data/glob.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/glob.h).
