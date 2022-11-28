#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for a specific release,
# supporting streaming or downloading from S3
# - the config file consists of file names (or URLs), with Q2 minima and cross sections
###################

# RELEASE TAG AND RECO DIR: ###########################
detector_config="epic_arches"
# detector_config="epic_brycecanyon"
release="22.11.2"
releaseDir="S3/eictest/EPIC/RECO/$release/$detector_config/DIS/NC"
filter='eicrecon'
# filter='juggler'
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
  S3 Directory:    $releaseDir

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
Q2minima=( 1000 100  10  1  )
Q2maxima=( 0    1000 100 10 )
########################################################################

# download files from S3
function q2subdir { echo $1/minQ2=$2; }
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" ]; then
  status "downloading files from S3..."
  for Q2min in ${Q2minima[@]}; do
    echo " sourceDir = `q2subdir $sourceDir $Q2min`"
    echo " targetDir = `q2subdir $targetDir $Q2min`"
    s3tools/generate-s3-list.sh  `q2subdir $sourceDir $Q2min`   | \
      { if [ $limit -gt 0 ]; then head -n$limit; else cat; fi } | \
      grep -E $filter                                           | \
      s3tools/download.sh `q2subdir $targetDir $Q2min`
  done
fi

# build a config file
status "build config file..."
mkdir -p $targetDir
if [ -z "$configFile" ]; then configFile=$targetDir/files.config; fi
> $configFile.list
for (( i=0; i<${#Q2minima[@]}; i++)); do
  Q2min=${Q2minima[$i]}
  Q2max=${Q2maxima[$i]}
  crossSection=$(s3tools/read-xsec-table.sh "pythia8:$energy/minQ2=$Q2min")
  case $mode in
    d) listScript=s3tools/generate-local-list.sh; listDir=`q2subdir $targetDir $Q2min`; ;;
    c) listScript=s3tools/generate-local-list.sh; listDir=`q2subdir $targetDir $Q2min`; ;;
    s) listScript=s3tools/generate-s3-list.sh;    listDir=`q2subdir $sourceDir $Q2min`; ;;
    *) echo "ERROR: unknown mode" >&2; exit 1;    ;;
  esac
  $listScript $listDir $crossSection $Q2min $Q2max            | \
    grep -v UNKNOWN                                           | \
    { if [ $limit -gt 0 ]; then head -n$limit; else cat; fi } | \
    tee -a $configFile.list
done
s3tools/generate-config-file.rb $configFile $energy $configFile.list

# finalize
popd
status "done building config file at:"
echo "     $configFile"
echo ""
if [ -n "$(grep UNKNOWN $configFile.list)" ]; then
  >&2 echo "ERROR: missing some cross sections"
  exit 1
fi
