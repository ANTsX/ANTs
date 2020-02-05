# Use Ubuntu 16.04 LTS
FROM nvidia/cuda:9.1-runtime-ubuntu16.04

# Prepare environment
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    curl \
                    bzip2 \
                    ca-certificates \
                    build-essential \
                    autoconf \
                    libtool \
                    pkg-config \
                    wget \
                    libboost-all-dev \
                    unzip \
                    libgl1-mesa-dev \
                    libglu1-mesa-dev \
                    freeglut3-dev \
                    mesa-utils \
                    g++ \
                    jq \
                    tar \
                    zip \
                    gcc \
                    make \
                    python \
                    zlib1g-dev \
                    imagemagick \
                    software-properties-common \
                    git && \
    curl -sL https://deb.nodesource.com/setup_10.x | bash - && \
    apt-get install -y --no-install-recommends \
      nodejs && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
ADD https://cmake.org/files/v3.11/cmake-3.11.4-Linux-x86_64.sh /cmake-3.11.4-Linux-x86_64.sh

ADD . /tmp/ants/source
RUN mkdir /opt/cmake \
  && sh /cmake-3.11.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license \
  && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake \
    && cd /tmp/ants/source \
    && mkdir -p /tmp/ants/build \
    && cd /tmp/ants/build \
    && mkdir -p /opt/ants-latest \
    && git config --global url."https://".insteadOf git:// \
    && cmake -DBUILD_TESTING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/opt/ants-latest /tmp/ants/source \
    && make -j2 \
    && cd ANTS-build \
    && make install \
    && rm -rf /tmp/ants \
    && rm -rf /opt/cmake /usr/local/bin/cmake

ENV ANTSPATH="/opt/ants-latest/bin" \
    PATH="/opt/ants-latest/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants-latest/lib:$LD_LIBRARY_PATH"
