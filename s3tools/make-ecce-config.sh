#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for a specific release,
# supporting streaming or downloading from S3
# - the config file consists of file names (or URLs), with Q2 minima and cross sections
###################

# RELEASE TAG AND RECO DIR: ###########################
release="prop.6/prop.6.0"
releaseDir="S3/eictest/ECCE/MC/$release/SIDIS/pythia6"
#######################################################

# usage:
if [ $# -lt 2 ]; then
  echo """
  USAGE: $0 [energy] [mode(d/s/c)] [limit(optional)] [outputFile(optional)]

   - [energy]: 5x41
               5x100
               10x100
               10x275
               18x275

       - NOTE: for release 'prop.6/prop.6.0', only 5x41 and 18x275 are available
               
   - [mode]:   s - make config file for streaming from S3
               d - download from S3, then make the local config file
               c - just make the local config file, for local files
   
   - [limit]   integer>0 : only stream/download this many files per Q2 min
               0         : stream/download all files
               default=5
   
   - [outputFile]: output file name (optional)
                   - default name is based on release version 
                   - relative paths will be relative to main dir
   
   Examples: $0 5x41 d       # download
             $0 18x275 s     # stream

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
mode=$2
limit=5
outFile=""
if [ $# -ge 3 ]; then limit=$3; fi
if [ $# -ge 4 ]; then outFile=$4; fi

# cd to the main directory 
pushd $(dirname $(realpath $0))/..

# settings #############################################################
sourceDir="$releaseDir/ep-$energy"
targetDir="datarec/ecce/$release/$energy"
Q2min=1 # FIXME: assumed, so far this script only looks at the general Q2 
        # production, and it doesn't matter if this is the *correct* Q2min;
        # this Q2min only matters when you want to combine datasets with
        # different Q2 minima (see `make-athena-config.sh`)
########################################################################

# download files from S3
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" ]; then
  status "downloading files from S3..."
  if [ $limit -gt 0 ]; then
    s3tools/generate-s3-list.sh "$sourceDir" | head -n$limit | s3tools/download.sh "$targetDir"
  else
    s3tools/generate-s3-list.sh "$sourceDir" |                 s3tools/download.sh "$targetDir"
  fi
fi

# build a config file
status "build config file..."
mkdir -p $targetDir
if [ -z "$outFile" ]; then configFile=$targetDir/files.config
else configFile=$outFile; fi
> $configFile
crossSection=$(s3tools/read-xsec-table.sh $energy $Q2min)
if [ "$mode" == "d" -o "$mode" == "c" ]; then
  s3tools/generate-local-list.sh "$targetDir" 0 $crossSection $Q2min | tee -a $configFile
elif [ "$mode" == "s" ]; then
  if [ $limit -gt 0 ]; then
    s3tools/generate-s3-list.sh "$sourceDir" 0 $crossSection $Q2min | head -n$limit | tee -a $configFile
  else
    s3tools/generate-s3-list.sh "$sourceDir" 0 $crossSection $Q2min | tee -a $configFile
  fi
else
  echo "ERROR: unknown mode"
  exit 1
fi

# PATCH: convert config file to one-line-per-Q2min format
status "reformatting config file to one-line-per-Q2min format..."
mv -v $configFile{,.bak}
s3tools/reformat-config.sh $configFile{.bak,}

# output some info
#status "files in target directory:"
#tree $targetDir
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
