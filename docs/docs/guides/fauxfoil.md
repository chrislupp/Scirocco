---
title: Use the XFOIL-style command line
description: fauXFOIL runs Scirocco behind XFOIL's menus, prompts, output and polar files.
---

# Use the XFOIL-style command line

`fauXFOIL` is an interactive program that looks and behaves like XFOIL 6.99:
the same menus, prompts, command abbreviations, iteration printout, plots
and polar file layout. Scirocco does the analysis. If you already drive XFOIL
by hand or from a script, you can use `fauXFOIL` the same way.

It's built by default (`SCIROCCO_BUILD_CLI=ON`) and installed next to the
library:

```bash
cmake -S scirocco -B scirocco/build -DCMAKE_BUILD_TYPE=Release
cmake --build scirocco/build -j --target fauXFOIL
scirocco/build/fauXFOIL
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

When `fauXFOIL` is built with Qt 6, plots open in a window laid out like
XFOIL's: white lines on black, the same axes, legends and colours. The
window is native, so it needs no X11 on macOS. Resize it freely; plots are
redrawn as vectors at the screen's full resolution. Lettering uses
fauXFOIL's own single-stroke font, which is drawn with the same thin pen as
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
- **Plot size.** `SIZE` in `OPER`, `VPLO` and `GDES` sets the plot-object
  size (default 10) for the plots that follow.

### Cursor input

Commands that ask XFOIL for cursor input take mouse clicks and key presses
in the plot window. The window takes keyboard focus while it waits, shows a
crosshair, and where XFOIL shows a **Done** button so does fauXFOIL
(clicking it is the same as typing `D`).

- **Zoom.** `Z` (any menu) zooms the window onto two clicked corners; `U`
  unzooms. Each corner is printed in plot and page ("abs", inches)
  coordinates, as XFOIL prints them.
- **Annotation.** `ANNO` opens XFOIL's annotation menu: text (`C`), slanted
  text (`S`), Greek "math" characters (`M`, where `a` is α, `b` is β, ...),
  point symbols (`P`, 14 types chosen with `T`), lines (`L`) and arrows
  (`A`), and the character width factor (`W`).
- **Velocity profiles.** `BL i` draws `i` boundary-layer velocity profiles
  (21 by default) equally spaced round the airfoil under the Cp plot and
  lists δ*, Re_θ and H_k at each; `BLC` draws them where you click, and
  `BLWT` sets their scale. The profile shapes are textbook families fitted
  to the solver's shape factor: the Pohlhausen quartic when laminar, a
  power law when turbulent.
- **Blowup.** `BLOW` in `VPLO` sets the axis limits to a clicked box; in
  `GDES` it plots only that region of the airfoil. `WIND` in `GDES` moves
  the region with the cursor: `I` zooms in, `O` out, and `P` or a click
  pans. `RESE` restores the full view.
- **Distance.** `DIST` in `GDES` prints the distance between two clicked
  points.

### Reference data and plot options

- `CREF` overlays measured Cp(x) from a file of "x Cp" lines on each Cp
  plot. As in XFOIL, the file name and a label are asked for each time; an
  empty answer keeps the previous file.
- `FREF` lists measured Mach, Re, α, CL, CD and CM beside the computed
  values.
- `RGET f` reads a reference polar in XFOIL's format (CD-CL, α-CL, α-CM and
  xtr-CL blocks, each ended by a `999 999` line) and plots it with symbols
  on the polar plot; `RDEL` deletes it. `PPAX` sets the polar plot's axis
  limits, and `APLO` plots the airfoil stored with each polar.
- In `VPLO`, `SYMB` marks the stations, `LABE` hides the header and titles,
  `CLIP` clips curves to the axes, `DUMP f` writes the plotted variable to a
  file and `OVER f` overlays "x value" data from a file.
- In `GDES`, `TICK` and `TSIZ` mark the points, `INPL` plots the airfoil at
  one inch per chord, `OVER f` overlays another airfoil with its
  parameters, and `CPLO` plots the surface curvature.

![Polar plot of NACA 2412 at Re = 10⁶](/img/fauxfoil-polar.png)

As in XFOIL, `PLOP` and `G` turn plotting off. `HARD` still works then,
because plots are built either way.

### Build with plots

The plot window needs Qt 6 (the Widgets module). CMake uses it when it finds
it and builds a text-only `fauXFOIL` otherwise; the configure step says
which.

```bash
brew install qt                  # macOS
sudo apt install qt6-base-dev    # Debian, Ubuntu
```

Set `SCIROCCO_CLI_GRAPHICS=OFF` to build without Qt even when it's
installed. On Linux, plots need an X11 or Wayland display. Without one, or
with `FAUXFOIL_NOGRAPHICS=1` in the environment, `fauXFOIL` runs text-only.
`FAUXFOIL_SNAPSHOT_DIR=dir` saves every plot shown as `dir/plotNNN.png`.

## Scripts

`fauXFOIL` reads commands from standard input, so XFOIL batch scripts run
unchanged, including the usual `PLOP`, `G`, `<return>` opening that turns
off plotting:

```bash
fauXFOIL < polar_script.txt
```

As with XFOIL, a file name on the command line is loaded at startup
(`fauXFOIL myfoil.dat`).

## Commands

| Menu | Commands |
|---|---|
| Top level | `QUIT`, `OPER`, `GDES`, `NACA`, `LOAD`, `SAVE`, `PSAV`, `REVE`, `NORM`, `XYCM`, `PCOP`, `PANE`, `PPAR`, `PLOP`, `NAME`, `NINC` |
| `OPER` | `VISC`, `RE`, `MACH`, `TYPE` (type 1), `ITER`, `INIT`, `ALFA`, `CL`, `CLI`, `ASEQ`, `CSEQ`, `!`, `VPAR`, `SEQP`, `PACC`, `PGET`, `PWRT`, `PSUM`, `PLIS`, `PDEL`, `PPLO`, `APLO`, `PPAX`, `RGET`, `RDEL`, `GRID`, `CPIO`, `CREF`, `FREF`, `CPX`, `CPV`, `VPLO`, `ANNO`, `HARD`, `SIZE`, `CPMI`, `BL`, `BLC`, `BLWT`, `CPWR`, `DUMP`, `CPMN`, `NAME`, `NINC` |
| `VPAR` | `SHOW`, `N`, `NT`, `NB`, `WAKE`, `INIT`, `LAG`, `GB`, `CTR`, `REST` |
| `VPLO` | `H`, `DT`, `DB`, `UE`, `CF`, `CD`, `N`, `CT`, `RT`, `RTL`, `DUMP`, `OVER`, `X`, `Y`, `XPLO`, `BLOW`, `RESE`, `SIZE`, `ANNO`, `HARD`, `GRID`, `SYMB`, `LABE`, `CLIP` |
| `GDES` | `GSET`, `EXEC`, `ADEG`, `ARAD`, `TRAN`, `SCAL`, `UNIT`, `DERO`, `FLAP`, `DIST`, `CPLO`, `PLOT`, `INPL`, `BLOW`, `RESE`, `WIND`, `TSIZ`, `TICK`, `GRID`, `GPAR`, `OVER`, `SIZE`, `ANNO`, `HARD`, `NAME`, `NINC` |
| Any menu | `Z`, `U` |

Abbreviations work as in XFOIL (`a` for `ALFA`, `as` for `ASEQ`, `v` for
`VISC`, and so on). The menus list every XFOIL command; one that `fauXFOIL`
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
  fauXFOIL and the Scirocco version instead of XFOIL. `PACC` appends to an existing
  polar file after checking its Mach, Reynolds number and Ncrit, as XFOIL
  does.
- `CPWR` writes x and Cp; `DUMP` writes s, x, y, Ue, δ*, θ, Cf, H, H*, P,
  m and K on the surface and wake.

## How it differs from XFOIL

- **Plots.** The `N` plot has no Tollmien-Schlichting frequency lines
  (`FRPL`): XFOIL draws them from Orr-Sommerfeld stability data that
  Scirocco does not have. When no plot window shows the Cp plot, inviscid
  results are printed instead.
- **Not yet available.** Polar sorting and editing (`PSOR`, `PREM`, `PNAM`,
  `ASET`), Cp-min and hinge-moment columns (`CINC`, `HINC`), flap hinge
  moments and point velocities (`FMOM`, `FNEW`, `VELS`), the `GDES`
  geometry editors other than those listed above, and inverse design
  (`MDES`, `QDES`).
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
