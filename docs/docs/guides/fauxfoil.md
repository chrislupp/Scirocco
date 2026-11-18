---
title: Use the XFOIL-style command line
description: fauXfoil runs Scirocco behind XFOIL's menus, prompts, output and polar files.
---

# Use the XFOIL-style command line

`fauXfoil` is an interactive program that looks and behaves like XFOIL 6.99:
the same menus, prompts, command abbreviations, iteration printout and
polar file layout. Scirocco does the analysis. If you already drive XFOIL
by hand or from a script, you can use `fauXfoil` the same way.

It's built by default (`SCIROCCO_BUILD_CLI=ON`) and installed next to the
library:

```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j --target fauXfoil
scirocco/build/fauXfoil
```

## A session

```text
 fauXFOIL   c>  naca 2412

 fauXFOIL   c>  oper

.OPERi   c>  visc 1e6

.OPERv   c>  pacc
...
.OPERva   c>  aseq 0 2 1
```

Each Newton iteration of a viscous point prints XFOIL's block as it runs:

```text
 Side 1  free  transition at x/c =  0.6524   57
 Side 2  free  transition at x/c =  0.6788   58

   6   rms: 0.1015E-08   max: -.1467E-07   D at   59  1
       a =  0.000      CL =  0.2353
      Cm = -0.0516     CD =  0.00563   =>   CDf =  0.00442    CDp =  0.00121
```

The numbers come from Scirocco's own Newton iteration, read through
`Glob::iteration_observer` (see [FoilState](../api/foil-state.md)). `rms`
and `max` are the size of the step just taken, normalized as XFOIL
normalizes them, and `RLX` is Scirocco's under-relaxation factor when it is
below 1. Scirocco converges on its own residual tolerance, so it often takes
an iteration or two more than XFOIL's `rms` criterion would. The transition
location solves the same amplification equation as the solver's transition
residual.

## Scripts

`fauXfoil` reads commands from standard input, so XFOIL batch scripts run
unchanged, including the usual `PLOP`, `G`, `<return>` opening that turns
off graphics:

```bash
fauXfoil < polar_script.txt
```

As with XFOIL, a file name on the command line is loaded at startup
(`fauXfoil myfoil.dat`).

## Commands

| Menu | Commands |
|---|---|
| Top level | `QUIT`, `OPER`, `GDES`, `NACA`, `LOAD`, `SAVE`, `PSAV`, `REVE`, `NORM`, `XYCM`, `PCOP`, `PANE`, `PPAR`, `PLOP`, `NAME`, `NINC` |
| `OPER` | `VISC`, `RE`, `MACH`, `TYPE` (type 1), `ITER`, `INIT`, `ALFA`, `CL`, `CLI`, `ASEQ`, `CSEQ`, `!`, `VPAR`, `PACC`, `PGET`, `PWRT`, `PSUM`, `PLIS`, `PDEL`, `CPWR`, `DUMP`, `CPMN`, `NAME`, `NINC` |
| `VPAR` | `SHOW`, `N`, `NT`, `NB`, `WAKE`, `INIT`, `LAG`, `GB`, `CTR`, `REST` |
| `GDES` | `GSET`, `EXEC`, `ADEG`, `ARAD`, `TRAN`, `SCAL`, `UNIT`, `DERO`, `FLAP`, `NAME`, `NINC` |

Abbreviations work as in XFOIL (`a` for `ALFA`, `as` for `ASEQ`, `v` for
`VISC`, and so on). The menus list every XFOIL command; one that `fauXfoil`
lacks prints `command not available in Scirocco`.

`VPAR` settings map to Scirocco's [`Param`](../api/solver-params.md):
`N` sets `ncrit`, `LAG` sets `SlagK` and `Cuq`, `GB` sets `GA` and `GB`,
`CTR` sets `CtauC` and `CtauE`, and `WAKE` sets the wake length. `ITER`
sets `niglob`.

## Files

- `SAVE` and `PSAV` write coordinates in XFOIL's order (counterclockwise
  from the upper trailing edge) and number format. `LOAD` reads labeled,
  plain, ISES and MSES single-element files in either direction.
- Polar save files (`PACC`, `PWRT`) use XFOIL's layout column for column,
  so tools that read XFOIL polars read these. The header line names
  fauXfoil and the Scirocco version instead of XFOIL. `PACC` appends to an existing
  polar file after checking its Mach, Reynolds number and Ncrit, as XFOIL
  does.
- `CPWR` writes x and Cp; `DUMP` writes s, x, y, Ue, δ*, θ, Cf, H, H*, P,
  m and K on the surface and wake.

## How it differs from XFOIL

- **No graphics.** Plot commands are listed but unavailable. Inviscid
  results, which XFOIL shows only on the Cp plot, are printed instead.
- **Paneling.** `PANE` uses Scirocco's curvature-based spacing. Only the
  node count (`PPAR N`, default 160) can be changed. As in XFOIL, `LOAD`
  uses the file's points directly as panel nodes, so run `PANE` after
  loading a coarse file.
- **One Ncrit, no trips.** Scirocco uses one Ncrit for both surfaces and
  has no forced transition, so `XTR` stays at 1.0.
- **Polar types.** Only type 1 (fixed Mach and Reynolds number).
- **No polar dump files.** XFOIL's binary dump files are not written.
- **Different numbers.** Results come from Scirocco, so they are close to
  XFOIL's but not identical. For example, NACA 2412 at Re = 10⁶ and α = 4°
  gives CL = 0.7202 and CD = 0.00693, against XFOIL's 0.7146 and 0.00693.
