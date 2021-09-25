#!/bin/bash
# start an interactive shell in the singularity container

# set image directory to definitely be ./img
if [ -z "${BASH_SOURCE[0]}" ]; then imgDirDflt=$(dirname $(realpath $0))/img
else imgDirDflt=$(dirname $(realpath ${BASH_SOURCE[0]}))/img; fi

# set image file
imgFile=$imgDirDflt/largex-eic.sif

# start shell
singularity shell $imgFile
