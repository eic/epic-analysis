#!/bin/bash
# start an interactive shell in the singularity container

# determine top-level directory
if [ -z "${BASH_SOURCE[0]}" ]; then topDir=$(dirname $(realpath $0))/..
else topDir=$(dirname $(realpath ${BASH_SOURCE[0]}))/..; fi

# set image file
imgDir=$topDir/container/img
imgFile=$imgDir/sidis-eic.sif

# start shell
args=$@
singularity exec $imgFile bash -c "cd $topDir && source environ.sh && bash $args"
