---
title: Use the XFOIL-style command line
description: fauXfoil runs Scirocco behind XFOIL's menus, prompts, output and polar files.
---

# Use the XFOIL-style command line

`fauXfoil` is an interactive program that looks and behaves like XFOIL 6.99:
the same menus, prompts, command abbreviations, iteration printout, plots
and polar file layout. Scirocco does the analysis. If you already drive XFOIL
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

## Plots

When `fauXfoil` is built with Qt 6, plots open in a window laid out like
XFOIL's: white lines on black, the same axes, legends and colours. The
window is native, so it needs no X11 on macOS. Resize it freely; plots are
redrawn as vectors at the screen's full resolution. Lettering uses
fauXfoil's own single-stroke font, which is drawn with the same thin pen as
the lines, in the style of XFOIL's plotter text.

![Cp plot of NACA 2412 at Re = 10⁶ and α = 4°](/img/fauxfoil-cp.png)

- **Cp(x).** Drawn after every `ALFA`, `CL` or `CLI` point, and by `CPX`.
  Viscous Cp is blue on the upper surface and red on the lower, with the
  inviscid Cp dashed (`CPIO` toggles it) and the displacement surfaces
  drawn on the airfoil below. `CPMI` sets the top of the Cp axis and
  `GRID` adds a grid.
- **Sequences.** During `ASEQ` and `CSEQ` the Cp curves of all points are
  overlaid with a table of results. `SEQP` plots the polar instead while
  `PACC` is on.
- **Pressure vectors.** `CPV` draws the airfoil at its angle of attack with
  a pressure vector at each node.
- **Boundary layer.** `VPLO` opens XFOIL's menu: `H`, `DT`, `DB`, `UE`,
  `CF`, `CD`, `N`, `CT`, `RT` and `RTL` plot one variable on both surfaces
  and the wake; `X` and `Y` set axis limits, `XPLO` switches between x and
  arc length, `RESE` resets the limits and `GRID` toggles the grid.
- **Polars.** `PPLO` plots the stored polars: CL against CD, CL and CM
  against α, and CL against transition location.
- **Geometry.** `GDES` plots the buffer airfoil with its area, thickness,
  camber, leading-edge radius and trailing-edge angle, and redraws it after
  every change. `PLOT` redraws it; `GRID`, `TICK` and `GPAR` toggle the
  grid, node ticks and parameters.
- **Hardcopy.** `HARD` writes the current plot to `plot.ps` as a PostScript
  page, black on white; later hardcopies add pages, as in XFOIL.

![Polar plot of NACA 2412 at Re = 10⁶](/img/fauxfoil-polar.png)

As in XFOIL, `PLOP` and `G` turn plotting off. `HARD` still works then,
because plots are built either way.

### Build with plots

The plot window needs Qt 6 (the Widgets module). CMake uses it when it finds
it and builds a text-only `fauXfoil` otherwise; the configure step says
which.

```bash
brew install qt                  # macOS
sudo apt install qt6-base-dev    # Debian, Ubuntu
```

Set `SCIROCCO_CLI_GRAPHICS=OFF` to build without Qt even when it's
installed. On Linux, plots need an X11 or Wayland display. Without one, or
with `FAUXFOIL_NOGRAPHICS=1` in the environment, `fauXfoil` runs text-only.
`FAUXFOIL_SNAPSHOT_DIR=dir` saves every plot shown as `dir/plotNNN.png`.

## Scripts

`fauXfoil` reads commands from standard input, so XFOIL batch scripts run
unchanged, including the usual `PLOP`, `G`, `<return>` opening that turns
off plotting:

```bash
fauXfoil < polar_script.txt
```

As with XFOIL, a file name on the command line is loaded at startup
(`fauXfoil myfoil.dat`).

## Commands

| Menu | Commands |
|---|---|
| Top level | `QUIT`, `OPER`, `GDES`, `NACA`, `LOAD`, `SAVE`, `PSAV`, `REVE`, `NORM`, `XYCM`, `PCOP`, `PANE`, `PPAR`, `PLOP`, `NAME`, `NINC` |
| `OPER` | `VISC`, `RE`, `MACH`, `TYPE` (type 1), `ITER`, `INIT`, `ALFA`, `CL`, `CLI`, `ASEQ`, `CSEQ`, `!`, `VPAR`, `SEQP`, `PACC`, `PGET`, `PWRT`, `PSUM`, `PLIS`, `PDEL`, `PPLO`, `GRID`, `CPIO`, `CPX`, `CPV`, `VPLO`, `HARD`, `CPMI`, `CPWR`, `DUMP`, `CPMN`, `NAME`, `NINC` |
| `VPAR` | `SHOW`, `N`, `NT`, `NB`, `WAKE`, `INIT`, `LAG`, `GB`, `CTR`, `REST` |
| `VPLO` | `H`, `DT`, `DB`, `UE`, `CF`, `CD`, `N`, `CT`, `RT`, `RTL`, `X`, `Y`, `XPLO`, `RESE`, `GRID`, `HARD` |
| `GDES` | `GSET`, `EXEC`, `ADEG`, `ARAD`, `TRAN`, `SCAL`, `UNIT`, `DERO`, `FLAP`, `PLOT`, `GRID`, `TICK`, `GPAR`, `HARD`, `NAME`, `NINC` |

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

- **Plots.** There is no cursor input, so zooming (`Z`, `U`, `BLOW`),
  annotation (`ANNO`) and plot sizing (`SIZE`) are unavailable, as are
  boundary-layer profiles (`BL`), reference data overlays (`CREF`, `OVER`)
  and the Tollmien-Schlichting frequency lines on the `N` plot. When no
  plot window shows the Cp plot, inviscid results are printed instead.
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
