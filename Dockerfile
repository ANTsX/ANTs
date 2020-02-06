FROM ubuntu:bionic-20200112 as builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    software-properties-common \
                    build-essential \
                    apt-transport-https \
                    ca-certificates \
                    gnupg \
                    software-properties-common \
                    wget \
                    ninja-build \
                    git \
                    zlib1g-dev

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
    | apt-key add - \
  && apt-add-repository -y 'deb https://apt.kitware.com/ubuntu/ bionic main' \
  && apt-get update \
  && apt-get -y install cmake

ADD . /tmp/ants/source
RUN mkdir -p /tmp/ants/build \
    && cd /tmp/ants/build \
    && mkdir -p /opt/ants-latest \
    && git config --global url."https://".insteadOf git:// \
    && cmake \
      -GNinja \
      -DBUILD_TESTING=OFF \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX=/opt/ants-latest \
      /tmp/ants/source \
    && cmake --build . --parallel \
    && cd ANTS-build \
    && cmake --install .

FROM ubuntu:bionic-20200112
COPY --from=builder /opt/ants-latest /opt/ants-latest

ENV ANTSPATH="/opt/ants-latest/bin" \
    PATH="/opt/ants-latest/bin:$PATH" \
    LD_LIBRARY_PATH="/opt/ants-latest/lib:$LD_LIBRARY_PATH"
