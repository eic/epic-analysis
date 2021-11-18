#!/bin/bash
# get generated cross section from all hepmc files found in specified directory

# settings #########################
limiter=1000
evgenDir=S3/eictest/ATHENA/EVGEN/DIS/NC
####################################

# ensure evgenDir ends with `/` (so `mc find` results are correct)
evgenDir=$evgenDir/
evgenDir=$(echo $evgenDir | sed 's/\/\//\//g') # remove '//'
evgenDir_esc=$(echo $evgenDir | sed 's/\//\\\//g') # replace '/' -> '\/'

# get list of hepmc files
#mc find $evgenDir --name "*.hepmc" > evgenList.tmp
#mc find $evgenDir --name "*.hepmc" | grep -vE 'vtxfix|novtx' > evgenList.tmp ###### matches Brian's list
#mc find $evgenDir --name "*.hepmc" | grep vtxfix > evgenList.tmp
#mc find $evgenDir --name "*.hepmc" | grep novtx > evgenList.tmp

# cull
sed -i '/2_022.hepmc/d' evgenList.tmp

# pretty print the list
echo "pretty print of hepmc file names:"
echo "----------------------------------"
cat evgenList.tmp | sed 's/^.*\///g' | sed -s 's/_/ /g'| column -t
echo "----------------------------------"

# loop through hepmc files
while read hepmc; do
  xsecDir=$(dirname $(echo $hepmc | sed "s/$evgenDir_esc/datarec\/xsec\//g"))
  echo $xsecDir 
  # mkdir -p xsecDir ...
done < evgenList.tmp
#rm evgenList.tmp


