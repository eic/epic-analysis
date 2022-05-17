### based on https://hub.docker.com/r/eicweb/jug_xl
FROM eicweb/jug_xl:4.0-deathvalley-1.5T-stable

### ld paths
ENV LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH

### install delphes
WORKDIR /opt
RUN \
git clone https://github.com/delphes/delphes && \
cd delphes && \
make
ENV DELPHES_HOME=/opt/delphes
ENV LD_LIBRARY_PATH=$DELPHES_HOME:$LD_LIBRARY_PATH
ENV PYTHONPATH=$DELPHES_HOME/python:$PYTHONPATH
ENV LIBRARY_PATH=$DELPHES_HOME:$LIBRARY_PATH

### copy sample root file for tests (source path is w.r.t. Dockerfile)
# COPY datarec/fromDuane/dis-10x100-xm25.trimmed.root /data/example_10x100.root

### non-privileged user
RUN \
mkdir -p /tmp/.X11-unix && \
useradd -u 1000 --create-home sidis && \
usermod --shell /bin/bash sidis
USER sidis
WORKDIR /home/sidis
