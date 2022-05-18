#!/bin/bash
# generate list of hepmc files

# arguments and usage
if [ $# -lt 2 ]; then echo """
USAGE: $0 [energy] [Q2min] [limit]
- [energy] is beam energies, e.g., 10x100
- [Q2min] is the minimum Q2
- [limit] is used to limit the number of files (default=0, take all)
- output is to stdout, so you can either:
  - redirect it to a file
  - pipe through grep for further filtering
  - you may need to edit the grep filter in the script
"""; exit 2; fi
energy=$1
minQ2=$2
limit=0
if [ $# -ge 3 ]; then limit=$3; fi

evgenDir="S3/eictest/ATHENA/EVGEN/DIS/NC/$energy/minQ2=$minQ2"
if [ $limit -gt 0 ]; then
  mc ls $evgenDir | grep -E 'hepmc.gz$' | grep -v GiB | grep vtxfix | head -n$limit | sed "s;^.* ;$evgenDir/;g"
else
  mc ls $evgenDir | grep -E 'hepmc.gz$' | grep -v GiB | grep vtxfix | sed "s;^.* ;$evgenDir/;g"
fi
