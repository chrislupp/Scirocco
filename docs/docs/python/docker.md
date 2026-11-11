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

Images are published for tagged releases. Use a version tag, not `latest`:

```bash
docker run --rm -p 50051:50051 chrislupp/scirocco:<version>
docker run --rm -p 50052:50052 chrislupp/scirocco-kulfan:<version>
```

Until the first image is published, build them yourself.

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
