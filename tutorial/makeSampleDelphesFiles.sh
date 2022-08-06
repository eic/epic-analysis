#!/bin/bash
# produce sample delphes files, for tutorials

# SETTINGS --------------------------------------------------------------------
minQ2=1                            # minimum Q2
ebeamEn=10                         # electron energy
pbeamEn=100                        # proton energy
maxNumFiles=10                     # maximum number of hepmc files to process
genDir=datagen/forTutorial         # directory of pythia output hepmc files
recDir=datarec/forTutorial         # directory of delphes output trees
configFile=tutorial/delphes.config # config file location
# -----------------------------------------------------------------------------

if [ -z "$S3_SECRET_KEY" -o -z "$S3_ACCESS_KEY" ]; then
  echo "ERROR: need to set env vars S3_SECRET_KEY and S3_ACCESS_KEY"
  exit 1
fi
s3tools/add-host.sh

mkdir -p $genDir $recDir

function sep { echo "----- $1 -----"; }

sep "build list"
s3tools/generate-hepmc-list.sh ${ebeamEn}x${pbeamEn} $minQ2 $maxNumFiles | tee hepmc.list.tmp

sep "download"
while read hepmc; do mc cp $hepmc $genDir/; done < hepmc.list.tmp
sep "list $genDir/"
ls -lh $genDir

sep "run delphes"
for hepmc in $genDir/*.hepmc.gz; do ./run_delphes.sh $hepmc; done
sep "list $recDir/"
ls -lh $recDir

sep "make config file"
s3tools/make-fastsim-local-config.sh ${ebeamEn}x${pbeamEn} $minQ2 $recDir $configFile
echo ">> config file:"
cat $configFile

rm hepmc.list.tmp

echo ""
echo ""
sep "DONE"
echo ""
echo "-> HEPMC FILES in $genDir/:"
ls -lh $genDir
echo ""
echo "-> DELPHES OUTPUT TREES in $recDir/:"
ls -lh $recDir
echo ""
echo "-> CONFIG FILE: $configFile"
echo ""
