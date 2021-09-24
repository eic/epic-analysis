FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get install -yq \
    binutils \
    cmake \
    curl \
    dpkg-dev \
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
    xlibmesa-glu-dev

RUN \
mkdir -p /opt && \
curl -sL https://root.cern/download/root_v6.24.02.Linux-ubuntu20-x86_64-gcc9.3.tar.gz | tar -C /opt -zxf -

 
ENV ROOTSYS=/opt/root
ENV PATH=$ROOTSYS/bin:$PATH
ENV LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
ENV PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

ENV PYTHIA8=/opt/pythia
ENV LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

RUN \
curl -sL https://pythia.org/download/pythia83/pythia8306.tgz | tar -C /opt -zxf - && \
  cd /opt/pythia8306 && \
  ./configure --prefix=/opt/pythia && \
  make install && \
  cd - && \
  rm -rf /opt/pythia8245 && \
git clone https://github.com/delphes/delphes /opt/delphes && \
  cd /opt/delphes && \
  make && \
  cd -

ENV DELPHES_HOME=/opt/delphes
ENV LD_LIBRARY_PATH=$DELPHES_HOME:$LD_LIBRARY_PATH
ENV PYTHONPATH=$DELPHES_HOME/python:$PYTHONPATH
ENV LIBRARY_PATH=$DELPHES_HOME:$LIBRARY_PATH

