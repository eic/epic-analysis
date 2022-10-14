#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for fastsim files
# - see also s3tools/make-fastsim-S3-config.sh, to automate S3 downloading and delphes execution
###################

# usage:
if [ $# -lt 5 ]; then
  echo """
  USAGE: $0 [energy] [Q2min] [Q2max] [source_directory] [output_file_name]

   - [energy]: 5x41 5x100 10x100 10x275 18x275
   - [Q2min]: 1 10 100 1000
   - [Q2max]: 0 100 ...
   - [source_directory]: location of fastsim files
   - [output_file_name]: output config file name

  """
  exit 2
fi
energy=$1
Q2min=$2
Q2max=$3
sourceDir=$4
configFile=$5

# cd to the main directory 
pushd $(dirname $(realpath $0))/..

# build a config file
function status { echo ""; echo "[+] $1"; }
status "build config file..."
> $configFile.list
crossSection=$(s3tools/read-xsec-table.sh "pythia8:$energy/minQ2=$Q2min")
s3tools/generate-local-list.sh "$sourceDir" $crossSection $Q2min $Q2max | grep -v UNKNOWN | tee -a $configFile.list
s3tools/generate-config-file.rb $configFile $energy $configFile.list

# output some info
popd
status "done building config file at:"
echo "     $configFile"
echo ""
if [ -n "$(grep UNKNOWN $configFile.list)" ]; then
  >&2 echo "ERROR: missing some cross sections"
  exit 1
fi
