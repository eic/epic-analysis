#!/bin/bash
# wrapper script to execute delphes


# default card
cardfile=cards/delphes_EICmatrixv2_3T.tcl


# check environment
if [ -z "$DELPHES_HOME" ]; then
  echo "ERROR: you must source env.sh first"
  exit 1
fi

# arguments
if [ $# -lt 1 ]; then
  echo "USAGE $0 [hepmc or hepmc.gz file] [(optional)output file name] [(optional)tcl card]"
  exit 1
fi
infile=$1
if [ $# -ge 2 ]; then outfile=$2
else outfile=$(echo $infile|sed 's/datagen/datarec/g'|sed 's/.hepmc.*$/.root/g'); fi
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
  cat $infile | gunzip | DelphesHepMC3 $cardfile $outfile
else
  echo "hepmc file is not gunzipped"
  DelphesHepMC3 $cardfile $outfile $infile
fi
