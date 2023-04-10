#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# wrapper script to execute delphes


# default card
cardfile=deps/delphes_EIC/ATHENA.tcl

# delphes executable 
exeDelphes=$DELPHES_HOME/DelphesHepMC3

#####################################################


# check environment
if [ -z "$DELPHES_HOME" ]; then
  echo "ERROR: you must source environ.sh first"
  exit 1
fi

# arguments
if [ $# -lt 1 ]; then
  echo "USAGE $0 [hepmc or hepmc.gz file] [(optional)output file name] [(optional)tcl card]"
  exit 1
fi
infile=$1
if [ $# -ge 2 ]; then outfile=$2
else
  if [[ $infile =~ ^datagen ]]; then
    outfile=$(echo $infile|sed 's/^datagen/datarec/'|sed 's/\.hepmc.*$/.root/g')
    mkdir -p $(dirname $outfile)
  else
    outfile=$(echo $infile|sed 's/^.*\//datarec\//g'|sed 's/\.hepmc.*$/.root/g')
  fi
fi
if [ $# -ge 3 ]; then cardfile=$3; fi
echo "infile = $infile"
echo "outfile = $outfile"
echo "cardfile = $cardfile"

# cleanup
if [ -f "$outfile" ]; then
  echo "output root file exists...remove before continuing..."
  rm -i $outfile
fi

# run delphes
if [[ $infile == *"hepmc.gz"* ]]; then
  echo "hepmc file is gunzipped"
  cat $infile | gunzip | $exeDelphes $cardfile $outfile
else
  echo "hepmc file is not gunzipped"
  $exeDelphes $cardfile $outfile $infile
fi
