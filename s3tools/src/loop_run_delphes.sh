#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# runs Delphes in a loop over all files in the directories specified in $*
# - optionally runs multi-threaded: one thread per directory
# - assumes input is `datagen/_____` and will output to `datarec/_____`

inputs=$*

function status { echo ""; echo "[+] $1"; }

status "running Delphes (one thread per directory)"
function runDelphes { 
  echo "delphes log: " > $1/delphes.log
  for infile in $1/*.hepmc{,.gz}; do 
    if [ ! -f "$infile" ]; then continue; fi
    echo "RUN DELPHES ON $infile" >> $1/delphes.log
    deps/run_delphes.sh $infile 2>&1 >> $1/delphes.log
  done
  status "DONE running Delphes on directory $1"
}

for input in $inputs; do
  status "RUNNING DELPHES on files in $input"
  runDelphes $input & # comment out `&` if you want to run single threaded
done

status "WAIT FOR DELPHES"
echo "  - quit by running: while [ 1 ]; do pkill DelphesHepMC3; done"
echo "  - monitor progress in log files in another window with:"
echo "    find $genDir -name \"delphes.log\" | xargs tail -F"
wait
status "DONE RUNNING DELPHES"
