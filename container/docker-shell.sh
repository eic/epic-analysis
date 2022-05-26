#!/bin/bash
# start docker shell, for testing image builds

# image=sidiseic_dev:latest
image=cjdilks/sidis-eic:latest

homedir=/home/sidis
swdir=$homedir/sidis-eic

docker run \
  -it \
  --rm \
  -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
  -e DISPLAY \
  -v $(pwd):$swdir \
  -w $swdir \
  $image bash
  #-u $(id -u):$(id -g) \
