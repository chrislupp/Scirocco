---
sidebar_position: 6
---

# Adjoint sensitivities

Scirocco exposes derivatives of $C_l$, $C_d$, and $C_m$ with respect to every
panel-node coordinate and to all four operating-condition inputs ($\alpha$,
$Re$, $M_\infty$, ... target lift via $\alpha$). Internally it uses a
reverse-mode adjoint so the cost of all three gradients is roughly the cost
of **one** primal solve, independent of the number of design variables.

This tutorial shows how to call it and how to verify the gradients against
finite differences.

## Computing sensitivities

```cpp
scirocco::Scirocco foil("2412", 129);
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.Solve();

auto sens = foil.ComputeSensitivities();   // returns a Sensitivity struct
```

`ComputeSensitivities` must be called **after** `Solve()` — it consumes the
converged state.

The returned `Sensitivity` struct:

```cpp
struct Sensitivity {
    Eigen::MatrixXd dcl_dx;     // 2 × N — w.r.t. foil node (x, z)
    Eigen::MatrixXd dcd_dx;     // 2 × N
    Eigen::MatrixXd dcm_dx;     // 2 × N
    Eigen::MatrixXd dcdpi_dx;   // 2 × N — inviscid pressure drag (often useful in viscous decomposition)

    double dcl_dalpha, dcd_dalpha, dcm_dalpha;   // per degree
    double dcl_dMa,    dcd_dMa,    dcm_dMa;
    double dcl_dRe,    dcd_dRe,    dcm_dRe;
};
```

The same data is also exposed as direct member accessors on the `Scirocco`
class:

```cpp
const Eigen::MatrixXd& dcl_dx = foil.GetDclDx();
double dcl_dalpha             = foil.GetDclDalpha();
double dcd_dRe                = foil.GetDcdDRe();
```

These accessors read from a cached `Sensitivity` that
`ComputeSensitivities` populates internally. The free-standing struct
returned by `ComputeSensitivities` is the same data.

## Layout of `dcl_dx`

`dcl_dx` is `2 × N` where $N$ is the panel-node count:

```
dcl_dx(0, j)  =  d(cl) / d(x_j)   — sensitivity to the x-coordinate of node j
dcl_dx(1, j)  =  d(cl) / d(z_j)   — sensitivity to the z-coordinate of node j
```

Nodes are ordered the same way as `GetFoilCoords()`. To get a sensitivity
with respect to a parametric design variable $p$, chain through the
geometry Jacobian $\partial x / \partial p$ — that's exactly what the
[Philote MDO discipline](../optimization/server) does for NACA digits and
CST coefficients.

## Verifying with finite differences

A quick sanity check: bump $\alpha$ by a small amount and compare the
analytic gradient to a centered difference.

```cpp
scirocco::Scirocco foil("0012", 99);
foil.SetOper(4.0, -1.0, 1.0e6, 0.0);
foil.Solve();
double analytic = foil.GetDclDalpha();   // [per degree]

const double h = 1.0e-4;
double cl_p, cl_m;

foil.SetOper(4.0 + h, -1.0, 1.0e6, 0.0);
foil.Solve();
cl_p = foil.GetCl();

foil.SetOper(4.0 - h, -1.0, 1.0e6, 0.0);
foil.Solve();
cl_m = foil.GetCl();

double fd = (cl_p - cl_m) / (2.0 * h);

std::cout << "dCl/dalpha  analytic = " << analytic << "\n";
std::cout << "            FD       = " << fd       << "\n";
std::cout << "            rel err  = "
          << std::abs(analytic - fd) / std::abs(fd) << "\n";
```

You should see relative errors at the level of the finite-difference
truncation (≈ $10^{-7}$ for centered differences in `double`). For
coordinate gradients, the same pattern works — perturb one entry of
`coords` by $h$, rebuild the airfoil, re-solve.

If you need higher accuracy than finite differences afford, the
[`tests/sensitivity_reference_test.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/tests/sensitivity_reference_test.cpp)
unit test compares the adjoint against complex-step references; it serves
as a worked template for the more rigorous verification.

## Reusing the converged solve

`ComputeSensitivities` does not re-solve the primal — it reuses the state
already in the solver. So in an optimizer loop you call:

```cpp
for (each design iteration) {
    foil.SetOper(...);
    if (!foil.Solve()) break;
    auto sens = foil.ComputeSensitivities();
    feed (foil.GetCl(), foil.GetCd(), sens) into optimizer
}
```

The two calls together cost one primal + one adjoint, which is the right
shape for gradient-based optimization.

## Full source

API in
[`scirocco/scirocco.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/scirocco.h);
struct in
[`scirocco/data/sensitivity.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/data/sensitivity.h);
reference test in
[`tests/sensitivity_reference_test.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/tests/sensitivity_reference_test.cpp).

## Next

You now have the full primal + adjoint workflow. Continue to:


- [Optimization → Server](../optimization/server) — using these gradients
  inside OpenMDAO via Philote
