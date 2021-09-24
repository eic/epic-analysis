#!/bin/bash
# start docker shell, for testing image builds

homedir=/home/athena
swdir=$homedir/largex-eic

docker run \
  -it \
  --rm \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -e DISPLAY \
  -v $(pwd):$swdir \
  -w $swdir \
  largexeic_ci:latest bash \
  -c "cat $swdir/docker/.shellrc >> $homedir/.bashrc && exec bash"
