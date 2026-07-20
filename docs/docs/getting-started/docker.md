---
sidebar_position: 3
---

# Docker

Versioned Docker images are published to Docker Hub on tagged releases. Until
the first public image is available, use the source-build instructions below.
Do not rely on an unversioned `latest` tag for a reproducible tutorial.

## Published images

| Image | Port | Hosts |
| --- | --- | --- |
| [`chrislupp/scirocco`](https://hub.docker.com/r/chrislupp/scirocco) | 50051 | `SiroccoDiscipline` (aerodynamic analysis) |
| [`chrislupp/scirocco-kulfan`](https://hub.docker.com/r/chrislupp/scirocco-kulfan) | 50052 | `KulfanDiscipline` (CST geometry) |

For a published release, replace `<version>` with its release number and run
the analysis server:

```bash
docker pull chrislupp/scirocco:<version>
docker run --rm -p 50051:50051 chrislupp/scirocco:<version>
```

Pull and run the Kulfan geometry server:

```bash
docker pull chrislupp/scirocco-kulfan:<version>
docker run --rm -p 50052:50052 chrislupp/scirocco-kulfan:<version>
```

Both images run as a non-root user on `debian:stable-slim`.

## Build from source

From the repository root:

```bash
docker build -f docker/scirocco-server/Dockerfile -t scirocco-server .
docker build -f docker/kulfan-server/Dockerfile  -t kulfan-server .
```

Run the locally-built images:

```bash
docker run --rm -p 50051:50051 scirocco-server
docker run --rm -p 50052:50052 kulfan-server
```

With both containers running, install the tested Python environment and run
the same smoke-test used by the source build:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r examples/tutorials/requirements.txt
python examples/tutorials/quickstart.py
```

## Connecting from Python

Once a container is running, connect to it the same way as a locally-launched
server—just target the mapped host port. This fragment creates the analysis
component; connect its `x` input to a local Python geometry component or to the
Kulfan component shown in the [Python quickstart](./quickstart):

```python
import grpc
from philote_mdo.openmdao import RemoteExplicitComponent

channel = grpc.insecure_channel("localhost:50051")
analysis = RemoteExplicitComponent(
    channel=channel,
    naca="0012",
    npanel=100,
)
```

See [Server setup](../optimization/server) for the full discipline interface
and [CST design](../optimization/cst-design) for a worked optimization
example.
