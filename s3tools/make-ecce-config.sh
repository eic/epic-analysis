#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for a specific release,
# supporting streaming or downloading from S3
# - the config file consists of file names (or URLs), with Q2 minima and cross sections
###################

#### RELEASE TAG #########################################
tag=latest
##########################################################

### use tag to specify the release information and settings
case $tag in
  latest)
    release="22.1"
    releaseDir="S3/eictest/EPIC/Campaigns/$release/SIDIS/pythia6"
    eventEvalDir="" # latest files are in the top-level directory, not `eval_0000{0,2}`
    ;;
  legacy) ### older example ECCE release
    release="prop.5/prop.5.1/SIDIS";
    releaseDir="S3/eictest/ECCE/MC/$release/pythia6"
    eventEvalDir="eval_00000"
    ;;
  *)
    echo "ERROR: unknown production tag"
    exit 1
esac
# common settings for all releases
eventEvalFileRegex='.*g4event_eval.root'

# usage:
if [ $# -lt 3 ]; then
  echo "Querying S3 for available data directories..."
  echo """
  USAGE: $0 [energy] [local_dir] [mode] [limit(optional)] [config_file(optional)]

   - [energy]: beam energies; data from differing Q2 ranges are combined
               automatically, weighted by cross sections

              AVAILABLE ENERGIES ON S3 in $releaseDir
              ========================
$(mc ls $releaseDir | sed 's;.* ep-;                ;' | sed 's;/$;;' | sed 's;-.*;;g' | uniq )
              ========================
               
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

  CURRENT RELEASE: 
    release:      $release
    releaseDir:   $releaseDir
    eventEvalDir: $eventEvalDir

  NOTE: Lambda directories are ignored, but if you want one, append
        '-Lambda' to [energy]; for example:
          $0 18x275-Lambda ecce.lambdas s 3
  """
  exit 2
fi
energy=$1
locDir=$2
mode=$3
limit=5
configFile=""
if [ $# -ge 4 ]; then limit=$4; fi
if [ $# -ge 5 ]; then configFile=$5; fi

### get list of subdirectories associated to this beam energy; each subdirectory has a different Q2 range
echo "Querying S3 for available data directories..."
if [[ "$energy" =~ "-" ]]; then
  subdirList=$(mc ls $releaseDir | grep $energy | awk '{print $NF}' | sed 's;/$;;')
else
  subdirList=$(mc ls $releaseDir | grep $energy | awk '{print $NF}' | sed 's;/$;;' | grep -v Lambda)
fi
printf "\nSubdirectories:\n"
for subdir in $subdirList; do echo "  $subdir"; done

# cd to the main directory 
pushd $(dirname $(realpath $0))/..

# function to get the sourceDir, given subdir $1
function getSourceDir {
  echo "$releaseDir/$1/$eventEvalDir" | sed 's;/$;;' | sed 's;//;/;g';
}

# function to get the Q2 minimum, given subdir $1
function getQ2min {
  if   [[ "$1" =~ "q2-low"  ]]; then echo 1;
  elif [[ "$1" =~ "q2-high" ]]; then echo 100;
  else echo 1  # general Q2 
  fi
}
function getQ2max {
  if [[ "$1" =~ "q2-low"  ]]; then echo 100;
  else echo 0 # no maximum
  fi
}

# set destination directory
targetDir="datarec/$locDir/$release/$energy"

# print settings
printf "\nsource directories:\n"
for subdir in $subdirList; do echo "  $(getSourceDir $subdir)   Q2min = $(getQ2min $subdir)   Q2max = $(getQ2max $subdir)"; done
echo "targetDir = $targetDir"

# download files from S3
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" ]; then
  status "downloading files from S3..."
  for subdir in $subdirList; do
    sourceDir=$(getSourceDir $subdir)
    if [ $limit -gt 0 ]; then
      s3tools/generate-s3-list.sh "$sourceDir" | grep -E $eventEvalFileRegex | head -n$limit | s3tools/download.sh "$targetDir/$subdir"
    else
      s3tools/generate-s3-list.sh "$sourceDir" | grep -E $eventEvalFileRegex |                 s3tools/download.sh "$targetDir/$subdir"
    fi
  done
fi

# build a config file
status "build config file..."
mkdir -p $targetDir
if [ -z "$configFile" ]; then configFile=$targetDir/files.config; fi
> $configFile.list
for subdir in $subdirList; do
  crossSection=$(s3tools/read-xsec-table.sh "pythia6:$subdir")
  sourceDir=$(getSourceDir $subdir)
  Q2min=$(getQ2min $subdir)
  Q2max=$(getQ2max $subdir)
  if [ "$mode" == "d" -o "$mode" == "c" ]; then
    s3tools/generate-local-list.sh "$targetDir/$subdir" $crossSection $Q2min $Q2max | tee -a $configFile.list
  elif [ "$mode" == "s" ]; then
    if [ $limit -gt 0 ]; then
      s3tools/generate-s3-list.sh "$sourceDir" $crossSection $Q2min $Q2max | grep -E $eventEvalFileRegex | head -n$limit | tee -a $configFile.list
    else
      s3tools/generate-s3-list.sh "$sourceDir" $crossSection $Q2min $Q2max | grep -E $eventEvalFileRegex | tee -a $configFile.list
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
