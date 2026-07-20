# syntax=docker/dockerfile:1.7

# Kulfan geometry discipline gRPC server (port 50052).
#
#   docker build --target kulfan-server -t kulfan-server:dev \
#     -f docker/kulfan-server/Dockerfile .

ARG DEBIAN_BASE=debian:stable-slim
ARG PHILOTE_REPO=https://github.com/MDO-Standards/Philote-Cpp.git
ARG PHILOTE_REF=cd19666a23ddca009ee8b27ecdd40c0a0a2ea4b4

# ---------------------------------------------------------------------------
# builder: compile PhiloteCpp + Scirocco, copy binary to /out
# ---------------------------------------------------------------------------
FROM ${DEBIAN_BASE} AS builder

ARG PHILOTE_REPO
ARG PHILOTE_REF

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        pkg-config \
        ca-certificates \
        libeigen3-dev \
        libgrpc++-dev \
        libprotobuf-dev \
        protobuf-compiler \
        protobuf-compiler-grpc \
    && rm -rf /var/lib/apt/lists/*

# Build PhiloteCpp (not vendored in this repo) and consume its build tree.
RUN git clone --recurse-submodules "${PHILOTE_REPO}" /src/Philote-Cpp \
    && git -C /src/Philote-Cpp checkout "${PHILOTE_REF}" \
    && git -C /src/Philote-Cpp submodule update --init --recursive \
    && cmake -S /src/Philote-Cpp -B /src/Philote-Cpp/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTS=OFF \
        -DBUILD_EXAMPLES=OFF \
    && cmake --build /src/Philote-Cpp/build -j"$(nproc)"

# Build Scirocco. -DSCIROCCO_NATIVE_ARCH=OFF keeps the image portable across
# x86-64 hosts (matches CI).
COPY scirocco /src/scirocco
RUN cmake -S /src/scirocco -B /src/scirocco/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DSCIROCCO_BUILD_TESTS=OFF \
        -DSCIROCCO_BUILD_EXAMPLES=OFF \
        -DSCIROCCO_NATIVE_ARCH=OFF \
        -DPhiloteCpp_DIR=/src/Philote-Cpp/build \
    && cmake --build /src/scirocco/build -j"$(nproc)" \
        --target scirocco_kulfan_server \
    && mkdir -p /out/bin \
    && cp /src/scirocco/build/src/philote/scirocco_kulfan_server /out/bin/

# ---------------------------------------------------------------------------
# runtime
# ---------------------------------------------------------------------------
FROM ${DEBIAN_BASE} AS kulfan-server

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgrpc++1.51 \
        libprotobuf32 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --create-home --shell /usr/sbin/nologin app

COPY --from=builder /out/bin/scirocco_kulfan_server /usr/local/bin/scirocco_kulfan_server

USER app
WORKDIR /home/app

EXPOSE 50052
ENTRYPOINT ["/usr/local/bin/scirocco_kulfan_server"]
CMD ["0.0.0.0:50052"]
