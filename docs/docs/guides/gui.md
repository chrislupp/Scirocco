---
title: Browse results in Scirocco Viewer
description: A desktop viewer for the HDF5 files Hdf5Writer writes, with Cp, boundary-layer and history plots.
---

# Browse results in Scirocco Viewer

Scirocco Viewer is a desktop program for looking through the HDF5 files that
[`Hdf5Writer`](hdf5.md) writes. Open a file and step through its solves:
the airfoil and wake, pressure arrows on the surface, the Cp curve,
boundary-layer distributions and the cl, cd and cm history. It can also
export the sequence as a video.

![Scirocco Viewer showing a NACA 2412 at -2 degrees with Cp arrows and the Cp chart](/img/gui-overview.png)

The viewer reads the file only; it does not run Scirocco. Use it on the
output of a sweep, an optimization history or any other run that writes
one group per solve.

## Build it

The viewer is a separate Qt project in the repository's `gui/` directory.
It needs Qt 6.6 or later (Quick, Quick Controls, Svg and Charts) and HDF5
with its C++ library:

```bash
cmake -S gui -B gui/build -DCMAKE_BUILD_TYPE=Release
cmake --build gui/build -j
```

Set `CMAKE_PREFIX_PATH` to your Qt installation if CMake cannot find it
(with Homebrew, `$(brew --prefix qt)`). On macOS the build is an app bundle,
`gui/build/SciroccoViewer.app`; elsewhere it is `gui/build/SciroccoViewer`.
Scirocco's own build does not include the viewer.

## Open a file

Pass the file on the command line, or start the viewer and press
<kbd>Ctrl</kbd>+<kbd>O</kbd> (<kbd>⌘</kbd>+<kbd>O</kbd> on macOS) or
double-click the empty window:

```bash
gui/build/SciroccoViewer.app/Contents/MacOS/SciroccoViewer alpha_sweep.h5
```

Each `iteration_NNN` group in the file is one step. The slider and arrows
at the bottom move between them, and the readout beside them shows the
current α, Cl and Cd.

| Key | Action |
| --- | --- |
| <kbd>←</kbd> / <kbd>→</kbd> | Previous / next step |
| <kbd>Home</kbd> / <kbd>End</kbd> | First / last step |
| <kbd>Space</kbd> | Play or pause the sequence |
| <kbd>Ctrl</kbd>+<kbd>E</kbd> (<kbd>⌘</kbd>+<kbd>E</kbd>) | Export video |

## What it shows

The buttons down the left edge choose what is on screen. Hover over one to
see its name.

- **Charts** opens a panel under the airfoil with three plots: **Cp**
  (viscous, with the inviscid curve dashed), **Distributions** (cf, ue, θ,
  δ\*, Hk, Reθ and the amplification or shear-lag variable n/Cτ along the
  surface) and **History** (cl, cd and cm against step).
- **Results** lists the operating point, the force coefficients and whether
  the solve converged.
- **HUD** overlays the step, α, Cl and Cd on the airfoil, in four text
  sizes; useful in exported videos.
- **Panel nodes** marks the panel end points on the airfoil.
- **Cp arrows** draws the pressure coefficient as arrows normal to the
  surface. Arrows for positive Cp point at the surface and arrows for
  suction point away from it. Length follows |Cp| and colour follows Cp.
  The panel sets the line width, the length scale and the colour map.
- **Export video** writes the sequence to an MP4 file, in the light or dark
  theme, at a chosen frame rate and resolution. It needs
  [FFmpeg](https://ffmpeg.org/) installed.

The last button switches between the light and dark themes.

## Command-line options

| Option | What it does |
| --- | --- |
| `--theme light` | Start in the light theme (default `dark`) |
| `--set name=value` | Set a window property before showing it, for example `activeSection=1` (charts), `showCpArrows=true` or `showHud=true`; repeat for more than one |
| `--screenshot out.png` | Save the window to an image once it has drawn, then quit |

`--set` and `--screenshot` together draw a picture without touching the
mouse. With `QT_QPA_PLATFORM=offscreen` they need no display, and
`QT_SCALE_FACTOR=2` doubles the resolution. `scripts/gui_screenshot.sh`
draws the picture on this page that way.
