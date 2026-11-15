#!/usr/bin/env bash
# Exercise a locally loaded image; use an ephemeral host port.
set -euo pipefail
image="$1"
kind="$2"
case "$kind" in
    analysis) port=50051 ;;
    geometry) port=50052 ;;
    *) echo "Unknown server kind: $kind" >&2; exit 1 ;;
esac
container=$(docker run --rm -d -p "127.0.0.1::$port" "$image")
cleanup() {
    docker logs "$container" || true
    docker stop "$container" >/dev/null 2>&1 || true
}
trap cleanup EXIT
address=$(docker port "$container" "$port/tcp")
python3 scripts/smoke_server.py --address "$address" --kind "$kind"
