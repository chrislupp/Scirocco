---
sidebar_position: 5
---

# Sensitivities

Adjoint-computed gradients of $C_l$, $C_d$, $C_m$, and near-field
pressure drag with respect to panel-node coordinates and
operating-condition inputs. Headers:
[`scirocco/data/sensitivity.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/sensitivity.h)
and
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h).

## Calling

```cpp
#include <stdexcept>
#include <scirocco/core/paneling.h>

// Optional: reconstruct spline arc lengths when panel nodes are the inputs.
scirocco::RefreshPanelGeometry(foil.GetState().foil);
foil.SetOper(α, -1.0, Re, Ma);
if (!foil.Solve()) {
    throw std::runtime_error("Airfoil analysis did not converge");
}
auto sens = foil.ComputeSensitivities();
```

`ComputeSensitivities` runs the reverse-mode adjoint against the converged
primal in `foil.GetState()` and populates the returned struct. It also
caches the result internally so the direct accessors on `Scirocco` work
without re-running the adjoint.

## Common `Sensitivity` fields

This excerpt shows the primary force gradients; the header also declares
near-field drag, lift-target, and solved-angle derivatives discussed below.

```cpp
struct Sensitivity {
    // Coordinate gradients — 2 × N matrices matching foil.x
    Eigen::MatrixXd dcl_dx;
    Eigen::MatrixXd dcd_dx;
    Eigen::MatrixXd dcm_dx;
    Eigen::MatrixXd dcdpi_dx;   // d(cdpi)/d(x) — near-field pressure drag

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
adjoint cost is one back-substitution per output (four outputs:
$C_l$, $C_d$, $C_m$, $C_{d,\pi}$). This is
a constant number of adjoint solves as $N$ grows. Factorization, geometric
linearization, and back-substitution still grow with the number of nodes.

## Geometry and validity

Arc lengths have an explicit `ArcLengthMode`. The constructor preserves
Python's original paneling and selects `Fixed`: node-coordinate derivatives
hold those arc-length labels fixed, exactly as changing `mfoil.foil.x` alone
does. `RefreshPanelGeometry` selects `Spline` and reconstructs arc lengths
through the current nodes; its derivatives include all ten reintegration
passes. Philote uses this spline path. `RecomputeArcLength` selects
`Cumulative` straight-panel lengths. Each mode now uses its own matching
foil/wake derivative chain. These are panel-node derivatives; they do not
include the earlier curvature-based repaneling from source geometry.

Reference chord, moment reference point, wake extent and other model
parameters are held fixed when differentiating the exposed inputs. In
particular, moving an endpoint does not automatically redefine reference chord.

`ComputeSensitivities` requires a converged solve. Failed adjoint solves raise
an exception instead of substituting regularized, stale, or zero gradients.
The cached getters require a successful `ComputeSensitivities` after the last
solve, operating-condition/geometry change, or mutable `GetState()` access.
When retaining a mutable state reference, solve again after modifying it.

Lift-target derivatives include the angle adjustment needed to keep Cl equal
to the target. `dcl_dx`, `dcl_dMa`, and `dcl_dRe` are consequently zero in this
mode. New `d*_dcltarget` fields differentiate with respect to the target;
`dalpha_dx`, `dalpha_dMa`, `dalpha_dRe`, `dalpha_dcltarget`, and
`dalpha_dalpha` describe the resulting angle. The constraint uses the Schur
complement of the fixed-angle linearization and fails explicitly when the
lift/angle slope is singular.

For viscous trim with Python's default `redowake=false`, the input angle
sets the frozen wake geometry. Thus `d*_dalpha` is the derivative with respect
to that **initial angle**, including the resulting trim adjustment. With
`redowake=true`, the wake follows the solved angle and the converged solution
is independent of its initial-angle guess. Inviscid trim has the same
initial-guess independence.

`dcdpi_dx`, `dcdpi_dalpha`, `dcdpi_dMa`, and `dcdpi_dRe` differentiate the
actual near-field pressure-drag output in both inviscid and viscous modes.
The Philote discipline continues to expose `cl`, `cd`, and `cm` at prescribed
alpha; its existing input/output interface is unchanged.

Inviscid `cd` and its derivatives are zero, matching Python mfoil. Near-field
pressure drag `cdpi` is a separate output. The existing C++ force conventions are unchanged: pressure-integrated
Cl/Cdpi use reference chord, Cm uses its square, and total Cd retains the
Python/XFOIL unit-reference convention. Use unit-chord geometry when
conventional chord-normalized coefficients are required.

## Numerical precision

Distant linear-vortex streamfunction coefficients use a stable quadrature of
mfoil's log-distance integral to avoid cancellation for short panels. Python
and C++ use the same mathematical integral, and its coordinate derivatives
are checked against independent 60-digit references. The end-to-end Python
comparisons use `2e-11 + 2e-8*abs(reference)` for both unit and nonunit chord.

## See also

- Tutorial: [Check analytic sensitivities](../tutorials/sensitivities) — use
  the Philote Python binding and verify component and total derivatives

- Theory: [Discrete adjoint](../theory) — residual formulation and derivative chain
