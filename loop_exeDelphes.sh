#!/bin/bash
# for cross check study

for f in datagen/arc/crossCheck*.hepmc; do
  outfile=$(echo $f|sed 's/datagen/datarec/g;s/\.hepmc/.root/')
  exeDelphes.sh $f $outfile
done
