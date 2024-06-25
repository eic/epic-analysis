#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks


### EPIC-ANALYSIS
if [ -z "${BASH_SOURCE[0]}" ]; then
  export EPIC_ANALYSIS_HOME=$(dirname $(realpath $0))
else
  export EPIC_ANALYSIS_HOME=$(dirname $(realpath ${BASH_SOURCE[0]}))
fi
echo "EPIC_ANALYSIS_HOME = $EPIC_ANALYSIS_HOME"

### DELPHES
export DELPHES_HOME=$EPIC_ANALYSIS_HOME/deps/delphes
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
export MSTWPDF_HOME=$EPIC_ANALYSIS_HOME/deps/mstwpdf
echo "MSTWPDF found at $MSTWPDF_HOME"
export LD_LIBRARY_PATH=$MSTWPDF_HOME:$LD_LIBRARY_PATH

### ADAGE
export ADAGE_HOME=$EPIC_ANALYSIS_HOME/deps/adage
echo "ADAGE found at $ADAGE_HOME"
export LD_LIBRARY_PATH=$ADAGE_HOME/lib:$LD_LIBRARY_PATH

# prioritize EPIC_ANALYSIS libraries
export LD_LIBRARY_PATH=$EPIC_ANALYSIS_HOME/lib:$LD_LIBRARY_PATH
