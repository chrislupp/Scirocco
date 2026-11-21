#!/usr/bin/env bash
# Regenerate the fauXFOIL screenshots on the docs site for a coming release.
#
#   scripts/fauxfoil_screenshots.sh minor      # major | minor | patch
#
# The bump is the label the release PR will carry; the version it gives is
# worked out as release_version.py works it out. The committed HEAD is built
# in a temporary worktree with that version, so the checkout is untouched,
# and the plots are written to docs/static/img/fauxfoil-{cp,polar}.png.
# Needs Qt 6 (Widgets); set CMAKE_PREFIX_PATH if CMake cannot find it.
set -euo pipefail

bump="${1:-}"
case "$bump" in
    major | minor | patch) ;;
    *) echo "usage: $0 major|minor|patch" >&2; exit 2 ;;
esac

root=$(git rev-parse --show-toplevel)
cd "$root"
version=$(python3 - "$bump" <<'EOF'
import re, subprocess, sys
sys.path.insert(0, "scripts")
from release_version import calculate
cmake = open("scirocco/CMakeLists.txt").read()
current = re.search(r"project\(scirocco\s+VERSION\s+(\d+\.\d+\.\d+)\b", cmake)[1]
tags = subprocess.check_output(["git", "tag", "--list"], text=True).splitlines()
print(calculate(current, ["release", sys.argv[1]], tags)["new_version"])
EOF
)
echo "fauXFOIL screenshots for $version"

work=$(mktemp -d)
cleanup() {
    git -C "$root" worktree remove --force "$work/src" >/dev/null 2>&1 || true
    rm -rf "$work"
}
trap cleanup EXIT
git worktree add --quiet --detach "$work/src" HEAD
perl -pi -e "s/(project\(scirocco\s+VERSION\s+)\d+\.\d+\.\d+/\${1}$version/" \
    "$work/src/scirocco/CMakeLists.txt"

if [ -z "${CMAKE_PREFIX_PATH:-}" ] && command -v brew >/dev/null; then
    export CMAKE_PREFIX_PATH
    CMAKE_PREFIX_PATH=$(brew --prefix qt 2>/dev/null || true)
fi
cmake -S "$work/src/fauxfoil" -B "$work/build" -DCMAKE_BUILD_TYPE=Release \
    -DFAUXFOIL_GRAPHICS=ON -DFAUXFOIL_BUILD_TESTS=OFF >/dev/null
cmake --build "$work/build" -j >/dev/null

# The Cp plot of the first point is the first plot drawn, the polar the last.
mkdir -p "$work/shots" "$work/run"
(
    cd "$work/run"
    printf '%s\n' 'NACA 2412' OPER 'VISC 1e6' 'ALFA 4' \
        PACC '' '' 'ASEQ -4 12 1' PPLO '' '' '' QUIT |
        QT_QPA_PLATFORM=offscreen FAUXFOIL_SNAPSHOT_DIR="$work/shots" \
            "$work/build/fauXFOIL" >"$work/session.log" 2>&1
)
grep -q "fauXFOIL Version $version\$" "$work/session.log" || {
    echo "fauXFOIL did not report version $version" >&2; exit 1
}
shots=("$work"/shots/plot*.png)
if [ "${#shots[@]}" -lt 2 ]; then
    echo "expected the Cp and polar plots, got ${#shots[@]} plot(s)" >&2
    tail -20 "$work/session.log" >&2
    exit 1
fi
cp "${shots[0]}" docs/static/img/fauxfoil-cp.png
cp "${shots[${#shots[@]} - 1]}" docs/static/img/fauxfoil-polar.png
echo "wrote docs/static/img/fauxfoil-cp.png and fauxfoil-polar.png"
