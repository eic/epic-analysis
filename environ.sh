#!/bin/bash

### Delphes
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
export DELPHES_HOME=$(pwd)/deps/delphes
if [ -f "$DELPHES_HOME/DelphesEnv.sh" ]; then
  cd $DELPHES_HOME
  source DelphesEnv.sh # source Delphes environment
  export PATH=$PATH:$DELPHES_HOME
  cd -
  echo "Delphes found at $DELPHES_HOME"
else
  echo "WARNING: Delphes is not found at $DELPHES_HOME"
fi

### MSTWPDF
export MSTWPDF_HOME=$(pwd)/mstwpdf
echo "MSTWPDF found at $MSTWPDF_HOME"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MSTWPDF_HOME"
