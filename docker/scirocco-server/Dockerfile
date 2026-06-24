# syntax=docker/dockerfile:1.7

# Scirocco discipline gRPC server (port 50051).
#
#   docker build --target scirocco-server -t scirocco-server:dev \
#     -f docker/scirocco-server/Dockerfile .

ARG DEBIAN_BASE=debian:stable-slim
ARG PHILOTE_REPO=https://github.com/chrislupp/Philote-Cpp.git
ARG PHILOTE_REF=v0.5.0

# ---------------------------------------------------------------------------
# builder: compile PhiloteCpp + Scirocco, install binary to /out
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

# Build and install PhiloteCpp (not vendored in this repo).
RUN git clone --recurse-submodules "${PHILOTE_REPO}" /src/Philote-Cpp \
    && git -C /src/Philote-Cpp checkout "${PHILOTE_REF}" \
    && git -C /src/Philote-Cpp submodule update --init --recursive \
    && cmake -S /src/Philote-Cpp -B /src/Philote-Cpp/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_TESTS=OFF \
        -DBUILD_EXAMPLES=OFF \
    && cmake --build /src/Philote-Cpp/build -j"$(nproc)" \
    && cmake --install /src/Philote-Cpp/build

# Build Scirocco. -DSCIROCCO_NATIVE_ARCH=OFF keeps the image portable across
# x86-64 hosts (matches CI).
COPY scirocco /src/scirocco
RUN cmake -S /src/scirocco -B /src/scirocco/build \
        -DCMAKE_BUILD_TYPE=Release \
        -DSCIROCCO_BUILD_TESTS=OFF \
        -DSCIROCCO_BUILD_EXAMPLES=OFF \
        -DSCIROCCO_NATIVE_ARCH=OFF \
    && cmake --build /src/scirocco/build -j"$(nproc)" \
        --target scirocco_server \
    && mkdir -p /out/bin \
    && cp /src/scirocco/build/src/philote/scirocco_server /out/bin/

# ---------------------------------------------------------------------------
# runtime
# ---------------------------------------------------------------------------
FROM ${DEBIAN_BASE} AS scirocco-server

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgrpc++1.51 \
        libprotobuf32 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --system --create-home --shell /usr/sbin/nologin app

COPY --from=builder /out/bin/scirocco_server /usr/local/bin/scirocco_server

USER app
WORKDIR /home/app

EXPOSE 50051
ENTRYPOINT ["/usr/local/bin/scirocco_server"]
CMD ["0.0.0.0:50051"]
