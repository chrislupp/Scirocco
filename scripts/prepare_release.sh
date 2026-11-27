#!/usr/bin/env bash
# Prepare a release and open its PR into main.
#
#   scripts/prepare_release.sh minor        # stable: major | minor | patch
#   scripts/prepare_release.sh minor rc     # prerelease: alpha | beta | rc
#
# Branches release/<tag> from origin/develop, merges origin/main in,
# regenerates the fauXFOIL screenshots (stable releases), makes every
# versioned change (scripts/release_plan.py prepare: CMake version,
# CHANGELOG, docs snapshot), builds the docs, checks the result, commits,
# pushes and opens the PR with its labels. The PR's CI then tests the exact
# commit that ships; merging it publishes the release (release.yml).
# Needs gh, npm and, for screenshots, Qt 6 (see fauxfoil_screenshots.sh).
#
# For testing: PREPARE_RELEASE_BASE starts from another commit instead of
# origin/develop, and PREPARE_RELEASE_NO_PUSH=1 stops after the commit.
set -euo pipefail
base="${PREPARE_RELEASE_BASE:-origin/develop}"

bump="${1:-}"
qualifier="${2:-}"
case "$bump" in major | minor | patch) ;; *)
    echo "usage: $0 major|minor|patch [alpha|beta|rc]" >&2; exit 2 ;;
esac
case "$qualifier" in "" | alpha | beta | rc) ;; *)
    echo "usage: $0 major|minor|patch [alpha|beta|rc]" >&2; exit 2 ;;
esac

root=$(git rev-parse --show-toplevel)
cd "$root"
if ! git diff --quiet || ! git diff --cached --quiet; then
    echo "The working tree has uncommitted changes; commit or stash them." >&2
    exit 1
fi

git fetch --quiet --tags origin develop main
labels=("$bump")
if [ -n "$qualifier" ]; then labels+=(prerelease "$qualifier"); else labels+=(release); fi
labels_json=$(printf '%s\n' "${labels[@]}" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read().split()))')

git switch --quiet --detach "$base"
tag=$(python3 - "$labels_json" <<'EOF'
import json, re, subprocess, sys
sys.path.insert(0, "scripts")
from release_version import calculate
cmake = open("scirocco/CMakeLists.txt").read()
current = re.search(r"project\(scirocco\s+VERSION\s+(\d+\.\d+\.\d+)\b", cmake)[1]
tags = subprocess.check_output(["git", "tag", "--list"], text=True).split()
print(calculate(current, json.loads(sys.argv[1]), tags)["tag_version"])
EOF
)
branch="release/$tag"
echo "Preparing $tag on $branch"
git switch --quiet -c "$branch"
git merge --quiet --no-edit origin/main || {
    echo "origin/main does not merge cleanly into develop; resolve that first." >&2
    exit 1
}

if [ -z "$qualifier" ]; then
    scripts/fauxfoil_screenshots.sh "$bump"
fi
npm --prefix docs ci --silent
python3 scripts/release_plan.py prepare "$bump" ${qualifier:+"$qualifier"}
npm --prefix docs run --silent typecheck
npm --prefix docs run --silent build
python3 -m unittest discover -s scripts -p 'test_*.py'

git add -A scirocco/CMakeLists.txt CHANGELOG.md docs
git commit --quiet -m "Release: $tag"
RELEASE_LABELS="$labels_json" python3 scripts/release_plan.py check --base origin/main

header=$(RELEASE_LABELS="$labels_json" python3 scripts/release_plan.py check --base origin/main |
    sed -n 's/^changelog_header=//p')
notes=$(python3 scripts/release_plan.py notes "$header")
if [ "${PREPARE_RELEASE_NO_PUSH:-0}" = 1 ]; then
    echo "Prepared $tag on $branch; not pushed (PREPARE_RELEASE_NO_PUSH=1)."
    exit 0
fi
git push --quiet -u origin "$branch"
label_args=()
for label in "${labels[@]}"; do label_args+=(--label "$label"); done
# shellcheck disable=SC2016 # the backticks are Markdown, not a command
gh pr create --base main --head "$branch" --title "Release $tag" "${label_args[@]}" \
    --body "$(printf 'Releases **%s**. Merging publishes it (see `.github/workflows/release.yml`); use a merge commit.\n\n%s\n' "$tag" "$notes")"
