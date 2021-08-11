#!/bin/bash

# run `source env.sh /path/to/delphes/repo`, otherwise pass no
# argument to use the default location; it may be useful to
# set your own default location
if [ $# -eq 1 ]; then
  delphesDir=$1
else
  # default delphes location
  delphesDir=${HOME}/builds/delphes
fi
echo "delphes repository expected at $delphesDir"

pushd $delphesDir > /dev/null
if [ -f "DelphesEnv.sh" ]; then

  # source delphes environment
  source DelphesEnv.sh
  echo "DELPHES_HOME set to $DELPHES_HOME"
  export PATH=$PATH:$DELPHES_HOME
  echo "success!"
  popd > /dev/null

  # symlink delphes external stuff (so that we can run analysis code here)
  ln -sf $DELPHES_HOME/external ./

else
  echo "ERROR: DelphesEnv.sh not found, is the repo really at ${delphesDir}?"
  popd > /dev/null
fi

