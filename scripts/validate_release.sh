#!/usr/bin/env bash
# Validate the prepared checkout without publishing. Run from the repo root.
# Requires CMake/Eigen/GTest/HDF5, PhiloteCpp, npm, Docker, and the Python
# reference + smoke-client dependencies. PHILOTE_CPP_DIR is a build-tree package.
set -euo pipefail
: "${PHILOTE_CPP_DIR:?Set PHILOTE_CPP_DIR to the pinned PhiloteCpp build}"
review_build="${RELEASE_BUILD_DIR:-build-release}"
python3 -m unittest discover -s scripts -p 'test_*.py'
python3 -m unittest discover -s mfoil -p 'test_*.py'
cmake -S scirocco -B "$review_build" -DCMAKE_BUILD_TYPE=Release \
    -DSCIROCCO_NATIVE_ARCH=OFF -DSCIROCCO_BUILD_PHILOTE=ON \
    -DPhiloteCpp_DIR="$PHILOTE_CPP_DIR" -DSCIROCCO_BUILD_TESTS=ON \
    -DSCIROCCO_BUILD_SENSITIVITY_DEBUG_TESTS=ON -DSCIROCCO_HDF5_OUTPUT=ON
cmake --build "$review_build" -j2
ctest --test-dir "$review_build" --output-on-failure --no-tests=error
npm --prefix docs ci
npm --prefix docs run typecheck
npm --prefix docs run build
# Keep these images local until *both* pass. Publication retags these exact IDs.
for kind in analysis geometry; do
    if [[ "$kind" == analysis ]]; then
        dockerfile=docker/scirocco-server/Dockerfile
    else
        dockerfile=docker/kulfan-server/Dockerfile
    fi
    docker build -f "$dockerfile" -t "scirocco-release:$kind" .
    bash scripts/smoke_image.sh "scirocco-release:$kind" "$kind"
done
