#!/bin/bash
# download and build the singularity image


# set default image directory
if [ -z "${BASH_SOURCE[0]}" ]; then imgDirDflt=$(dirname $(realpath $0))/img
else imgDirDflt=$(dirname $(realpath ${BASH_SOURCE[0]}))/img; fi


# usage
if [ $# -lt 1 ]; then 
  echo """
  download and build the singularity image

  USAGE: $0 [image directory] [(optional) temp directory]
    
    - [image directory]: directory to store the Singularity image
      - run \"$0 .\" to use the default directory, which is:
        $imgDirDflt
      - otherwise specify your preferred directory
      - note that the image size is of order 2 GB
    
    - [(optional) temp directory]: temporary directory
      - used while building the SIF file, may need a few GB for a few minutes
      - your default directory is probably $([[ -z "$TMPDIR" ]] && echo "/tmp" || echo "$TMPDIR")
        specify a different directory, if there is not enough space there
  """
  exit 1
fi


# cleanup
if [ -d "$imgDirDflt" ] || [ -f "$imgDirDflt" ]; then
  echo ""
  echo "This is not your first installation attempt, clean up first by"
  echo "answering y/n to the following prompts:"
  echo ""
  rm -ri $imgDirDflt
fi


# arguments
if [ "$1" == "." ]; then
  imgDir=$imgDirDflt
  mkdir -p $imgDir
else
  imgDir=$1
  if [ ! -d "$imgDir" ]; then
    echo "ERROR: $imgDir does not exist"
    exit 1
  fi
  ln -sf $imgDir $imgDirDflt
fi
if [ $# -ge 2 ]; then export TMPDIR=$2; fi


# print settings
imgFile=$imgDir/largex-eic.sif
echo """
image directory = $imgDir
image file = $imgFile
temp directory = $([[ -z "$TMPDIR" ]] && echo "/tmp" || echo "$TMPDIR")
"""


# pull image from dockerhub
singularity pull $imgFile docker://cjdilks/largex-eic:dev && \
echo "SUCCESS!"

