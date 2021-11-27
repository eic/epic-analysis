#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for fastsim files
###################

# usage:
if [ $# -lt 4 ]; then
  echo """
  USAGE: $0 [energy] [Q2min] [source_directory] [output_file_name]

   - [energy]: 5x41 5x100 10x100 10x275 18x275
   - [Q2min]: 1 10 100 1000
   - [source_directory]: location of fastsim files
   - [output_file_name]: output config file name

  """
  exit 2
fi
energy=$1
Q2min=$2
sourceDir=$3
configFile=$4

# cd to the main directory 
pushd $(dirname $(realpath $0))/..

# build a config file
function status { echo ""; echo "[+] $1"; }
status "build config file..."
> $configFile
crossSection=$(s3tools/read-xsec-table.sh $energy $Q2min)
s3tools/generate-local-list.sh "$sourceDir" 0 $crossSection $Q2min | tee -a $configFile

# PATCH: convert config file to one-line-per-Q2min format
status "reformatting config file to one-line-per-Q2min format..."
mv -v $configFile{,.bak}
s3tools/reformat-config.sh $configFile{.bak,}

# output some info
popd
status "done building config file at:"
echo "     $configFile"
status "run root macros with parameters:"
echo "     '(\"$configFile\",$(echo $energy|sed 's/x/,/'))'"
echo ""
if [ -n "$(grep UNKNOWN $configFile)" ]; then
  >&2 echo "ERROR: missing some cross sections"
  exit 1
fi
