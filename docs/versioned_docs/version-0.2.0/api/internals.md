---
sidebar_position: 8
---

# Internals

The solver kernels under `scirocco/core/`, `scirocco/util/`, and
`scirocco/properties/`. You can include them, but they aren't a stable API
and may change between minor versions. Use [`Scirocco`](./scirocco) and
[`FoilState`](./foil-state) first.

## Inviscid

| Header | One-liner |
| --- | --- |
| [`core/inviscid.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/inviscid.h) | `BuildGamma`, `GetUeInv`, `SolveInviscid`, `InitThermo`, `ClearSolution` |
| [`core/panel_influence.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/panel_influence.h) | Closed-form linear-vortex influence-coefficient kernels |
| [`core/panel_influence_deriv.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/panel_influence_deriv.h) | Analytic derivatives of the above for the adjoint chain |
| [`core/inviscid_deriv.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/inviscid_deriv.h) | $R_U$, $R_x$ partials of the inviscid Newton system |
| [`core/paneling.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/paneling.h) | Curvature-graded paneler, `MakePanels`, `RefreshPanelGeometry` |
| [`core/spline.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/spline.h) | Cubic-spline interpolation used by the paneler |

## Viscous boundary layer

| Header | One-liner |
| --- | --- |
| [`core/boundary_layer.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/boundary_layer.h) | Integral BL residuals; Thwaites laminar initialization; station-by-station march |
| [`core/coupled_solver.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/coupled_solver.h) | `SolveCoupled`, `BuildGlobSys`, `UpdateState`, `CalcUeM`, `StoreAdjointFactorization` — the top-level viscous Newton loop |
| [`core/transition.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/transition.h) | `MarchAmplification`, `UpdateTransition` — $e^N$ method machinery |
| [`core/viscous_geom_deriv.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/viscous_geom_deriv.h) | Geometric chain in the viscous adjoint (foil-x → wake-s → BL stations) |
| [`core/viscous_residual_deriv.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/viscous_residual_deriv.h) | Residual partials of the BL system |
| [`core/wake.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/wake.h) | Wake geometry generation and march |
| [`properties/bl_properties.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/properties/bl_properties.h) | $H^*$, $C_f$, $C_D$ closure relations; Sutherland viscosity; Kármán–Tsien constants |

## Forces and adjoint

| Header | One-liner |
| --- | --- |
| [`core/forces.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/forces.h) | `CalcForce` (pressure + friction integration), `GetDistributions` |
| [`core/adjoint.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/core/adjoint.h) | Full reverse-mode adjoint kernel |

## Utilities

| Header | One-liner |
| --- | --- |
| [`util/math_utils.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/util/math_utils.h) | Sign, atan2 helpers, vector-norm utilities |
| [`util/linearization.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/util/linearization.h) | Linearization-recording helpers for the adjoint chain |
| [`config.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/config.h) | Eigen type aliases (`VectorXd`, `MatrixXd`, `Vector2d`, …) and global constants |

## Philote integration

| Header | One-liner |
| --- | --- |
| [`philote_discipline.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/philote_discipline.h) | `SiroccoDiscipline` — the analysis discipline that the `scirocco_server` binary serves |
| [`kulfan_discipline.h`](https://github.com/chrislupp/Scirocco/blob/develop/scirocco/include/scirocco/kulfan_discipline.h) | [`KulfanDiscipline`](./kulfan-discipline) — has its own page |

## Missing something?

If you find yourself depending on one of these, open an issue describing
what you need. Patterns that come up often get added to `Scirocco`.
