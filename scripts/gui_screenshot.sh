#!/usr/bin/env bash
# Regenerate the Scirocco Viewer screenshot on the docs site.
#
#   scripts/gui_screenshot.sh
#
# The committed HEAD is built in a temporary worktree, so the checkout is
# untouched: Scirocco's HDF5 example writes a NACA 2412 alpha sweep, and the
# viewer opens it with the charts, Cp arrows and HUD showing and saves its
# window to docs/static/img/gui-overview.png. It draws offscreen at twice
# the window size, so it needs no display. Needs Qt 6.6+ and HDF5; set
# CMAKE_PREFIX_PATH if CMake cannot find them.
set -euo pipefail

root=$(git rev-parse --show-toplevel)
cd "$root"

work=$(mktemp -d)
cleanup() {
    git -C "$root" worktree remove --force "$work/src" >/dev/null 2>&1 || true
    rm -rf "$work"
}
trap cleanup EXIT
git worktree add --quiet --detach "$work/src" HEAD

if [ -z "${CMAKE_PREFIX_PATH:-}" ] && command -v brew >/dev/null; then
    export CMAKE_PREFIX_PATH
    CMAKE_PREFIX_PATH=$(brew --prefix qt 2>/dev/null || true)
fi

cmake -S "$work/src/scirocco" -B "$work/build-scirocco" \
    -DCMAKE_BUILD_TYPE=Release -DSCIROCCO_HDF5_OUTPUT=ON \
    -DSCIROCCO_BUILD_TESTS=OFF -DSCIROCCO_NATIVE_ARCH=OFF >/dev/null
cmake --build "$work/build-scirocco" -j --target hdf5_output_example >/dev/null
mkdir -p "$work/run"
(cd "$work/run" && "$work/build-scirocco/hdf5_output_example" >/dev/null)

cmake -S "$work/src/gui" -B "$work/build-gui" -DCMAKE_BUILD_TYPE=Release \
    >/dev/null
cmake --build "$work/build-gui" -j >/dev/null

viewer="$work/build-gui/SciroccoViewer"
[ -f "$viewer" ] && [ -x "$viewer" ] || viewer="$work/build-gui/SciroccoViewer.app/Contents/MacOS/SciroccoViewer"
QT_QPA_PLATFORM=offscreen QT_SCALE_FACTOR=2 "$viewer" \
    --set activeSection=1 --set showCpArrows=true --set showHud=true \
    --screenshot "$work/gui-overview.png" "$work/run/alpha_sweep.h5"
cp "$work/gui-overview.png" docs/static/img/gui-overview.png
echo "wrote docs/static/img/gui-overview.png"
