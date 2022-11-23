#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks


### SIDIS-EIC
if [ -z "${BASH_SOURCE[0]}" ]; then
  export SIDIS_EIC_HOME=$(dirname $(realpath $0))
else
  export SIDIS_EIC_HOME=$(dirname $(realpath ${BASH_SOURCE[0]}))
fi
echo "SIDIS_EIC_HOME = $SIDIS_EIC_HOME"

### DELPHES
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH
export DELPHES_HOME=$SIDIS_EIC_HOME/deps/delphes
if [ -f "$DELPHES_HOME/DelphesEnv.sh" ]; then
  # set $LIBRARY_PATH to "", if unbound (for CI `eic/run-cvmfs-osg-eic-shell` payloads)
  LIBRARY_PATH_BIND=${LIBRARY_PATH:-}
  [ -z "$LIBRARY_PATH_BIND" ] && export LIBRARY_PATH=""
  # source Delphes environment
  source $DELPHES_HOME/DelphesEnv.sh
  export PATH=$PATH:$DELPHES_HOME
  echo "Delphes found at $DELPHES_HOME"
else
  echo "WARNING: Delphes is not found at $DELPHES_HOME"
fi

### MSTWPDF
export MSTWPDF_HOME=$SIDIS_EIC_HOME/deps/mstwpdf
echo "MSTWPDF found at $MSTWPDF_HOME"
export LD_LIBRARY_PATH=$MSTWPDF_HOME:$LD_LIBRARY_PATH

### ADAGE
export ADAGE_HOME=$SIDIS_EIC_HOME/deps/adage
echo "ADAGE found at $ADAGE_HOME"
export LD_LIBRARY_PATH=$ADAGE_HOME/lib:$LD_LIBRARY_PATH

# prioritize SIDIS_EIC libraries
export LD_LIBRARY_PATH=$SIDIS_EIC_HOME/lib:$LD_LIBRARY_PATH
