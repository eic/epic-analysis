#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# get generated cross section from all hepmc files found in specified directory

# settings #########################
limiter=50000
evgenDir=S3/eictest/ATHENA/EVGEN/DIS/NC
####################################

# ensure evgenDir ends with `/` (so `mc find` results are correct)
evgenDir=$evgenDir/
evgenDir=$(echo $evgenDir | sed 's/\/\//\//g') # remove '//'
evgenDir_esc=$(echo $evgenDir | sed 's/\//\\\//g') # replace '/' -> '\/', for sed

# get list of hepmc files
#mc find $evgenDir --name "*.hepmc" > evgenList.tmp
mc find $evgenDir --name "*.hepmc" | grep -vE 'vtxfix|novtx' > evgenList.tmp ###### matches Brian's logfile list
#mc find $evgenDir --name "*.hepmc" | grep vtxfix > evgenList.tmp    ## TODO: different cross section for `vtxfix` or `novtx`?
#mc find $evgenDir --name "*.hepmc" | grep novtx > evgenList.tmp     ## only `vtxfix` is in canyonlands 2.0

# cull
sed -i '/2_022.hepmc/d' evgenList.tmp

# pretty print the list
echo "pretty print of hepmc file names:"
echo "----------------------------------"
cat evgenList.tmp | sed 's/^.*\///g' | sed -s 's/_/ /g'| column -t
echo "----------------------------------"

# loop through hepmc files
while read hepmc; do
  echo "READ $hepmc"
  xsecDir=$(dirname $(echo $hepmc | sed "s/$evgenDir_esc/datarec\/xsec\//g"))
  outFile=$xsecDir/$(basename $hepmc).xsec
  mkdir -p $xsecDir
  echo "  outFile = $outFile"
  s3tools/src/get-cross-section.sh $hepmc $limiter | tee $outFile
done < evgenList.tmp
rm evgenList.tmp
#tree datarec/xsec
