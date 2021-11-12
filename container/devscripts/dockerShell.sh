#!/bin/bash
# start docker shell, for testing image builds

image=largexeic_dev:latest
#image=cjdilks/largex-eic:dev

homedir=/home/athena
swdir=$homedir/largex-eic

docker run \
  -it \
  --rm \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -e DISPLAY \
  -v $(pwd):$swdir \
  -w $swdir \
  $image bash
  #-u $(id -u):$(id -g) \
