---
title: Docker
description: Run the Scirocco servers in containers instead of building them.
---

# Docker

Both servers come as Docker images, so you can skip the Philote-Cpp build.
You still need the [Python environment](./index.md#set-up) on the client
side.

## Published images

| Image | Port | Server |
| --- | --- | --- |
| [`chrislupp/scirocco`](https://hub.docker.com/r/chrislupp/scirocco) | 50051 | Analysis |
| [`chrislupp/scirocco-kulfan`](https://hub.docker.com/r/chrislupp/scirocco-kulfan) | 50052 | CST geometry |

| Tag | Tracks | Updated |
| --- | --- | --- |
| `<version>`, for example `0.2.0` | That release | Never |
| `<major>.<minor>`, `<major>`, `latest` | The newest stable release in that series | On each stable release |
| `next` | The `develop` branch | On each push to `develop` that passes the image smoke test |

Use a version tag, not `latest`, so your results don't change under you.
The commands on this site use the tag that matches the docs you're reading:
`next` in the unreleased docs, and the release number in a release's docs.

<DockerRun />

Stop the servers with `docker stop scirocco scirocco-kulfan`.

The images are built for `linux/amd64`. On Apple silicon, Docker Desktop
runs them under emulation: they work, but solve more slowly than a native
build, and Docker prints a platform warning when it starts them.

Use `next` to try changes before they are released. It moves with
`develop` and can change between pulls, so pull it again before each run:

```bash
docker pull chrislupp/scirocco:next
docker run --rm -p 50051:50051 chrislupp/scirocco:next
```

The two images are published separately, so their `next` tags can briefly
come from different commits.

## Build the images

From the repository root:

```bash
docker build -f docker/scirocco-server/Dockerfile -t scirocco-server .
docker build -f docker/kulfan-server/Dockerfile  -t kulfan-server .

docker run --rm -p 50051:50051 scirocco-server
docker run --rm -p 50052:50052 kulfan-server
```

Both images run as a non-root user on `debian:stable-slim`.

## Connect

Containers look just like local servers. Connect to the mapped ports, then
run the [first analysis](./index.md#run-one-analysis):

```bash
python examples/tutorials/quickstart.py
```

The images are built without HDF5, so the `hdf5_file` option isn't
available in them. Build the servers yourself to record histories.
