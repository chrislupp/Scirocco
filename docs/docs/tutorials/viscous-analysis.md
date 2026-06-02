---
sidebar_position: 2
---

# Viscous analysis

Pass a non-zero Reynolds number and Scirocco switches from a linear inviscid
solve to a coupled inviscid + integral boundary-layer Newton iteration. This
tutorial mirrors
[`scirocco/examples/viscous_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/viscous_example.cpp)
and adds context on what the new outputs mean.

## The program

```cpp
#include <iomanip>
#include <iostream>

#include "scirocco/scirocco.h"

int main() {
    scirocco::Scirocco foil("0012", 129);   // 129 panels is plenty for inviscid + BL

    const double Re = 1.0e6;
    const double Ma = 0.0;                  // incompressible

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Alpha      CL        CD        CDf       CDp       CM\n";

    for (double alpha = 0.0; alpha <= 10.0; alpha += 2.0) {
        foil.SetOper(alpha, -1.0, Re, Ma);
        bool converged = foil.Solve();

        std::cout << std::setw(7) << alpha << "  "
                  << std::setw(8) << foil.GetCl()  << "  "
                  << std::setw(8) << foil.GetCd()  << "  "
                  << std::setw(8) << foil.GetCdf() << "  "
                  << std::setw(8) << foil.GetCdp() << "  "
                  << std::setw(8) << foil.GetCm();
        if (!converged) std::cout << "  (not converged)";
        std::cout << "\n";
    }
}
```

## What changed from inviscid

Just one slot in `SetOper`:

```cpp
foil.SetOper(alpha, -1.0, 1.0e6, 0.0);   // Re = 1e6 turns viscosity on
```

That single change cascades through the solve. After the linear inviscid
solution is built, Scirocco:

1. Marches the integral BL equations from the stagnation point downstream on
   each surface (and along the wake), using Thwaites for laminar starts.
2. Tracks the amplification factor $\tilde n(s)$ along each laminar segment
   and triggers transition when $\tilde n$ crosses $n_\text{crit}$ (default
   9 — see [Concepts → Transition](../concepts/transition)).
3. Couples the BL displacement effect back into the inviscid system through
   the edge velocity $u_e$ and runs a **Newton iteration** over the joint
   state.

The full nonlinear solve is the reason `Solve()` now actually returns
`true`/`false`: it can fail to converge. See
[Concepts → Pseudo-transient continuation](../concepts/ptc) for what to turn
on when it does.

## Reading the output

```text
  Alpha      CL        CD        CDf       CDp       CM
  -----   -------   -------   -------   -------   -------
 0.0000    0.0000    0.0055    0.0042    0.0013    0.0000
 2.0000    0.2412    0.0058    0.0042    0.0016   -0.0001
 4.0000    0.4837    0.0066    0.0041    0.0025   -0.0005
 ...
```

Total drag is now decomposed:

| Quantity | Getter | Meaning |
| --- | --- | --- |
| $C_d$ | `GetCd()` | Total drag, integrated to far-wake momentum deficit |
| $C_{d,f}$ | `GetCdf()` | Friction component (integral of $C_f \cdot \text{tangent}$) |
| $C_{d,p}$ | `GetCdp()` | Pressure component (form drag from BL displacement) |

In attached flow they're roughly comparable; as you approach stall,
$C_{d,p}$ grows faster than $C_{d,f}$.

The four distribution getters return $\theta(s)$, $\delta^*(s)$, $C_p(s)$,
$C_f(s)$:

```cpp
const auto& theta = foil.GetTheta();      // momentum thickness
const auto& ds    = foil.GetDeltaStar();  // displacement thickness
const auto& cp    = foil.GetCp();
const auto& cf    = foil.GetCf();
```

These are the standard outputs you'd plot in an integral BL paper.

## Convergence

`Solve()` returns `false` if the Newton iteration runs out of steps without
reducing the residual under `state.param.rtol` (default `1e-11`). When that
happens, three knobs help in order of preference:

1. Turn on PTC: `foil.GetState().param.use_ptc = true;` (see
   [PTC](../concepts/ptc)).
2. Refine panels — try 199 instead of 129, especially near the leading edge.
3. Raise `n_crit` toward 11–13 if the airfoil is at a high $Re$ in a clean
   wind tunnel; lower toward 5–7 for a noisy free flight.

You can also bump verbosity to see what the inner solver is doing:

```cpp
foil.GetState().param.verb = 1;   // 0=silent (default), 1=summary, 2=verbose
```

## Full source

[`scirocco/examples/viscous_example.cpp`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/examples/viscous_example.cpp).

## Next: [Custom airfoil coordinates](./custom-airfoil) — bring your own geometry.
