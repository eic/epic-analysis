#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for a specific release,
# supporting streaming or downloading from S3
# - the config file consists of file names (or URLs), with Q2 minima and cross sections
###################

# RELEASE TAG AND RECO DIR: ###########################
release="deathvalley-v1.0"
releaseDir="S3/eictest/ATHENA/RECO/$release/DIS/NC"
#######################################################

# usage:
if [ $# -lt 3 ]; then
  echo """
  USAGE: $0 [energy] [local_dir] [mode] [limit(optional)] [config_file(optional)]

   - [energy]: 5x41      | - see below for available datasets
               5x100     | - data from different Q2minima are combined,
               10x100    |   weighted by cross sections
               10x275
               18x275
               
   - [local_dir]: output directory name: datarec/[local_dir]

   - [mode]:   s - make config file for streaming from S3
               d - download from S3, then make the local config file
               c - just make the local config file, for local files
   
   - [limit]   integer>0 : only stream/download this many files per Q2 min
               0         : stream/download all files
               default=5

   - [config_file]  name of the config file; if not specified, the
                    config file will be in datarec/[local_dir]

   See script for local and remote file path settings; they are
   configured for a specific set of data, but you may want to change
   them.

  CURRENT RELEASE: $release

  AVAILABLE DATA ON S3 (press ^C to abort S3 query):
  """
  mc tree $releaseDir
  exit 2
fi
energy=$1
locDir=$2
mode=$3
limit=5
configFile=""
if [ $# -ge 4 ]; then limit=$4; fi
if [ $# -ge 5 ]; then configFile=$5; fi

# cd to the main directory 
pushd $(dirname $(realpath $0))/..

# settings #############################################################
sourceDir="$releaseDir/$energy"
targetDir="datarec/$locDir/$release/$energy"
Q2minima=( 1000 100 10 1 )
Q2max=0 # no maximum
########################################################################

# download files from S3
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" ]; then
  status "downloading files from S3..."
  for Q2min in ${Q2minima[@]}; do
    if [ $limit -gt 0 ]; then
      s3tools/generate-s3-list.sh "$sourceDir/minQ2=$Q2min" | head -n$limit | s3tools/download.sh "$targetDir/minQ2=$Q2min"
    else
      s3tools/generate-s3-list.sh "$sourceDir/minQ2=$Q2min" |                 s3tools/download.sh "$targetDir/minQ2=$Q2min"
    fi
  done
fi

# build a config file
status "build config file..."
mkdir -p $targetDir
if [ -z "$configFile" ]; then configFile=$targetDir/files.config; fi
> $configFile.list
for Q2min in ${Q2minima[@]}; do
  crossSection=$(s3tools/read-xsec-table.sh "pythia8:$energy/minQ2=$Q2min")
  if [ "$mode" == "d" -o "$mode" == "c" ]; then
    s3tools/generate-local-list.sh "$targetDir/minQ2=$Q2min" $crossSection $Q2min $Q2max | grep -v UNKNOWN | tee -a $configFile.list
  elif [ "$mode" == "s" ]; then
    if [ $limit -gt 0 ]; then
      s3tools/generate-s3-list.sh "$sourceDir/minQ2=$Q2min" $crossSection $Q2min $Q2max | grep -v UNKNOWN | head -n$limit | tee -a $configFile.list
    else
      s3tools/generate-s3-list.sh "$sourceDir/minQ2=$Q2min" $crossSection $Q2min $Q2max | grep -v UNKNOWN | tee -a $configFile.list
    fi
  else
    echo "ERROR: unknown mode"
    exit 1
  fi
done
s3tools/generate-config-file.rb $configFile $energy $configFile.list

# output some info
#status "files in target directory:"
#tree $targetDir
popd
status "done building config file at:"
echo "     $configFile"
echo ""
if [ -n "$(grep UNKNOWN $configFile.list)" ]; then
  >&2 echo "ERROR: missing some cross sections"
  exit 1
fi
