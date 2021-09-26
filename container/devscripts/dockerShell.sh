#!/bin/bash
# start docker shell, for testing image builds

#image=largexeic_ci:latest
image=largexeic_dev:latest

homedir=/home/athena
swdir=$homedir/largex-eic

docker run \
  -it \
  --rm \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -e DISPLAY \
  -v $(pwd):$swdir \
  -w $swdir \
  $image bash \
  -c "cat $swdir/container/devscripts/.shellrc >> $homedir/.bashrc && exec bash"
