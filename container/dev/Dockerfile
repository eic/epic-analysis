FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# dependencies
RUN \
apt-get update && \
apt-get upgrade -yq && \
apt-get install -yq \
    binutils \
    cmake \
    curl \
    dpkg-dev \
    emacs \
    g++ \
    gcc \
    gfortran \
    git \
    graphviz-dev \
    libavahi-compat-libdnssd-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libftgl-dev \
    libglew1.5-dev \
    libglew-dev \
    libglu1-mesa-dev \
    libgsl0-dev \
    libjpeg-dev \
    libkrb5-dev \
    libldap2-dev \
    libmysqlclient-dev \
    libpcre3-dev \
    libpng-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxext-dev \
    libxft-dev \
    libxi-dev \
    libxml2-dev \
    libxmu-dev \
    libxpm-dev \
    libxt-dev \
    make \
    python \
    python-dev \
    rsync \
    tcl \
    wget \
    vim \
    xlibmesa-glu-dev

# general environment
ENV LS_COLORS="di=01;31"

# root
WORKDIR /opt
RUN \
curl -sL https://root.cern/download/root_v6.24.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz | tar -C /opt -zxf -
ENV ROOTSYS=/opt/root
ENV PATH=$ROOTSYS/bin:$PATH
ENV LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

# pythia
WORKDIR /opt
RUN \
curl -sL https://pythia.org/download/pythia83/pythia8306.tgz | tar -C /opt -zxf - && \
cd pythia8306 && \
./configure --prefix=/opt/pythia && \
make install && \
cd - && \
rm -rf /opt/pythia8306
ENV PYTHIA8=/opt/pythia
ENV LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

# delphes
WORKDIR /opt
RUN \
git clone https://github.com/delphes/delphes && \
cd delphes && \
make
ENV DELPHES_HOME=/opt/delphes
ENV LD_LIBRARY_PATH=$DELPHES_HOME:$LD_LIBRARY_PATH
ENV PYTHONPATH=$DELPHES_HOME/python:$PYTHONPATH
ENV LIBRARY_PATH=$DELPHES_HOME:$LIBRARY_PATH

# minIO client
WORKDIR /opt
RUN \
wget -Pbin https://dl.min.io/client/mc/release/linux-amd64/mc && \
chmod +x bin/mc
ENV PATH=$PATH:/opt/bin

# copy sample root file for tests (source path is w.r.t. Dockerfile)
COPY datarec/fromDuane/dis-10x100-xm25.trimmed.root /data/example_10x100.root

# non-privileged user
RUN \
mkdir -p /tmp/.X11-unix && \
useradd -u 1000 --create-home athena && \
usermod --shell /bin/bash athena
USER athena
WORKDIR /home/athena
