#!/bin/bash

# do `source env.sh /path/to/delphes/repo`; otherwise pass no
# argument to use the default location
if [ $# -eq 1 ]; then
  delphesDir=$1
else
  # default delphes location
  delphesDir=${HOME}/builds/delphes
fi

echo "delphes repository expected at $delphesDir"

pushd $delphesDir > /dev/null
if [ -f "DelphesEnv.sh" ]; then
  source DelphesEnv.sh
  echo "DELPHES_HOME set to $DELPHES_HOME"
  export PATH=$PATH:$DELPHES_HOME
  echo "success!"
else
  echo "ERROR: DelphesEnv.sh not found, is the repo really at ${delphesDir}?"
fi
popd > /dev/null
