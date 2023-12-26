FROM docker.io/library/debian:bookworm-slim AS base
FROM base AS builder

RUN \
    --mount=type=cache,sharing=private,target=/var/cache/apt \
    apt-get update && apt-get install -y g++-11 cmake make ninja-build git

RUN git config --global url.'https://'.insteadOf 'git://'

COPY . /usr/local/src/ants
WORKDIR /build

ARG CC=gcc-11 CXX=g++-11

RUN cmake \
    -GNinja \
    -DBUILD_TESTING=ON \
    -DRUN_LONG_TESTS=OFF \
    -DRUN_SHORT_TESTS=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=/opt/ants \
    /usr/local/src/ants
RUN cmake --build . --parallel
WORKDIR /build/ANTS-build
RUN cmake --install .

FROM base

COPY --from=builder /opt/ants /opt/ants
ENV PATH="/opt/ants/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants/lib:$LD_LIBRARY_PATH"

LABEL org.opencontainers.image.authors="ANTsX team" \
      org.opencontainers.image.url="https://stnava.github.io/ANTs/" \
      org.opencontainers.image.source="https://github.com/ANTsX/ANTs" \
      org.opencontainers.image.licenses="Apache License 2.0" \
      org.opencontainers.image.title="Advanced Normalization Tools" \
      org.opencontainers.image.description="ANTs is part of the ANTsX ecosystem (https://github.com/ANTsX). \
ANTs Citation: https://pubmed.ncbi.nlm.nih.gov/24879923"

