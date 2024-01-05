FROM ubuntu:jammy as builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
                    apt-transport-https \
                    bc \
                    build-essential \
                    ca-certificates \
                    gnupg \
                    ninja-build \
                    git \
                    software-properties-common \
                    wget

# install miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-1-Linux-$(uname -m).sh \
    && /bin/bash Miniconda3-py310_23.11.0-1-Linux-$(uname -m).sh -b -p /opt/conda \
    && rm Miniconda3-py310_23.11.0-1-Linux-$(uname -m).sh
ENV PATH=/opt/conda/bin:$PATH

# install cmake binary using conda for multi-arch support
RUN conda update -c defaults conda
RUN conda install -c conda-forge cmake

ADD . /tmp/ants/source
RUN mkdir -p /tmp/ants/build \
    && cd /tmp/ants/build \
    && mkdir -p /opt/ants \
    && git config --global url."https://".insteadOf git:// \
    && cmake \
      -GNinja \
      -DBUILD_TESTING=ON \
      -DRUN_LONG_TESTS=OFF \
      -DRUN_SHORT_TESTS=ON \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/ants \
      /tmp/ants/source \
    && cmake --build . --parallel \
    && cd ANTS-build \
    && cmake --install .

# Need to set library path to run tests
ENV LD_LIBRARY_PATH="/opt/ants/lib:$LD_LIBRARY_PATH"

RUN cd /tmp/ants/build/ANTS-build \
    && cmake --build . --target test

FROM ubuntu:jammy
COPY --from=builder /opt/ants /opt/ants

LABEL maintainer="ANTsX team" \
      description="ANTs is part of the ANTsX ecosystem (https://github.com/ANTsX). \
ANTs Citation: https://pubmed.ncbi.nlm.nih.gov/24879923"

ENV PATH="/opt/ants/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants/lib:$LD_LIBRARY_PATH"
RUN apt-get update \
    && apt install -y --no-install-recommends bc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /data

CMD ["/bin/bash"]
