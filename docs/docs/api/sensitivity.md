---
sidebar_position: 5
---

# Sensitivities

Adjoint-computed gradients of $C_l$, $C_d$, $C_m$, and the inviscid
pressure-drag component with respect to panel-node coordinates and
operating-condition inputs. Headers:
[`scirocco/data/sensitivity.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/sensitivity.h)
and
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h).

## Calling

```cpp
foil.SetOper(α, -1.0, Re, Ma);
foil.Solve();                          // must converge first
Sensitivity sens = foil.ComputeSensitivities();
```

`ComputeSensitivities` runs the reverse-mode adjoint against the converged
primal in `foil.GetState()` and populates the returned struct. It also
caches the result internally so the direct accessors on `Scirocco` work
without re-running the adjoint.

## The `Sensitivity` struct

```cpp
struct Sensitivity {
    // Coordinate gradients — 2 × N matrices matching foil.x
    Eigen::MatrixXd dcl_dx;
    Eigen::MatrixXd dcd_dx;
    Eigen::MatrixXd dcm_dx;
    Eigen::MatrixXd dcdpi_dx;   // d(cdpi)/d(x) — inviscid pressure-drag piece

    // Operating-condition gradients
    double dcl_dalpha = 0.0;    double dcd_dalpha = 0.0;    double dcm_dalpha = 0.0;
    double dcl_dMa    = 0.0;    double dcd_dMa    = 0.0;    double dcm_dMa    = 0.0;
    double dcl_dRe    = 0.0;    double dcd_dRe    = 0.0;    double dcm_dRe    = 0.0;
};
```

Layout of the `2 × N` coordinate matrices:

```
sens.dcl_dx(0, j)  =  d(cl) / d(x_j)
sens.dcl_dx(1, j)  =  d(cl) / d(z_j)
```

Node ordering matches `Scirocco::GetFoilCoords()` (TE-lower → LE → TE-upper).

`dalpha` derivatives are **per degree** to match the units of
`Oper::alpha`. `dMa` and `dRe` are per unit of $M_\infty$ and $Re$
respectively.

## Direct accessors on `Scirocco`

These read from the cached `Sensitivity` after `ComputeSensitivities`:

| Method | Returns |
| --- | --- |
| `GetDclDx()`, `GetDcdDx()`, `GetDcmDx()` | `const Eigen::MatrixXd&` (2 × N) |
| `GetDclDalpha()`, `GetDcdDalpha()`, `GetDcmDalpha()` | per-degree scalars |
| `GetDclDMa()`, `GetDcdDMa()`, `GetDcmDMa()` | per-Mach scalars |
| `GetDclDRe()`, `GetDcdDRe()`, `GetDcmDRe()` | per-Re scalars |

```cpp
double dCl_dα   = foil.GetDclDalpha();
const auto& g  = foil.GetDcdDx();         // 2 × N
double dCd_dx_j = g(0, j);
```

## Computational cost

The adjoint solve reuses the SparseLU factorization stored by
`StoreAdjointFactorization` during the primal Newton iteration, so the
adjoint cost is one back-substitution per output (three outputs:
$C_l$, $C_d$, $C_m$ — actually four with $C_{d,\pi}$). This is
**independent of $N$**, which is what makes panel-node gradients
practical.

## Adjoint ridge fallback

If the SparseLU factorization is numerically singular (typically at a
converged state with low-damping transition rows), the adjoint solve
adds a tiny diagonal ridge `Param::adjoint_ridge` (default `1e-10`) to
$R_V$ and retries. The bias on the resulting gradient is
$O(\varepsilon \cdot \|R_V\|)$ — small enough to be invisible in normal
work, large enough to disambiguate when it triggers.

If you suspect a ridge fired, set `state.param.verb = 1` before
`ComputeSensitivities` and watch for the warning line.

## See also

- Tutorial: [Adjoint sensitivities](../tutorials/sensitivities) — including the
  finite-difference verification recipe

- Theory: [adjoint.md](https://github.com/chrislupp/Scirocco/blob/develop/theory/adjoint.md) — full derivation
