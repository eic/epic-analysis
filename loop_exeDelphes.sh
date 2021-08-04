#!/bin/bash
# for cross check study
# - will run `nthread` parallel threads

nthread=4

mkdir -p tmp/log
njob=0
for f in datagen/arc/crossCheck*.hepmc; do
  if [ $njob -lt $nthread ]; then
    outfile=$(echo $f|sed 's/datagen/datarec/g;s/\.hepmc/.root/')
    logfile=tmp/log/$(echo $f|sed 's/^.*\///g').log
    exeDelphes.sh $f $outfile 2>&1 > $logfile &
    echo "processing $outfile (log $logfile)"
    let njob++
  else
    wait
    njob=0
  fi
done
