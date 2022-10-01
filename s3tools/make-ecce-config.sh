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
if [ $# -lt 2 ]; then
  echo "Querying S3 for available data directories..."
  echo """
  USAGE: $0 [energy] [mode(d/s/c)] [limit(optional)] [outputFile(optional)]

   - [energy]: beam energies; data from differing Q2 ranges are combined
               automatically, weighted by cross sections

              AVAILABLE ENERGIES ON S3 in $releaseDir
              ========================
$(mc ls $releaseDir | sed 's;.* ep-;                ;' | sed 's;/$;;' | sed 's;-.*;;g' | uniq )
              ========================
               
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

  CURRENT RELEASE: 
    release:      $release
    releaseDir:   $releaseDir
    eventEvalDir: $eventEvalDir
  """
  exit 2
fi
energy=$1
mode=$2
limit=5
outFile=""
if [ $# -ge 3 ]; then limit=$3; fi
if [ $# -ge 4 ]; then outFile=$4; fi

### get list of subdirectories associated to this beam energy; each subdirectory has a different Q2 range
echo "Querying S3 for available data directories..."
subdirListUnsorted=$(mc ls $releaseDir | grep $energy | awk '{print $NF}' | sed 's;/$;;')
### sort in order of decreasing Q2 minimum
function searchSubdirList { echo $subdirListUnsorted | sed 's; ;\n;g' | grep $*; }
subdirList="`searchSubdirList high` `searchSubdirList low` `searchSubdirList -v q2`" # take all 3 directories
# subdirList="`searchSubdirList high` `searchSubdirList -v q2`" # skip the 'q2-low` directory
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

# set destination directory
targetDir="datarec/ecce/$release/$energy"

# print settings
printf "\nsource directories:\n"
for subdir in $subdirList; do echo "  $(getSourceDir $subdir)   Q2min = $(getQ2min $subdir)"; done
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
if [ -z "$outFile" ]; then configFile=$targetDir/files.config
else configFile=$outFile; fi
> $configFile
for subdir in $subdirList; do
  crossSection=$(s3tools/read-xsec-table.sh "pythia6:$subdir")
  sourceDir=$(getSourceDir $subdir)
  Q2min=$(getQ2min $subdir)
  if [ "$mode" == "d" -o "$mode" == "c" ]; then
    s3tools/generate-local-list.sh "$targetDir/$subdir" 0 $crossSection $Q2min | tee -a $configFile
  elif [ "$mode" == "s" ]; then
    if [ $limit -gt 0 ]; then
      s3tools/generate-s3-list.sh "$sourceDir" 0 $crossSection $Q2min | grep -E $eventEvalFileRegex | head -n$limit | tee -a $configFile
    else
      s3tools/generate-s3-list.sh "$sourceDir" 0 $crossSection $Q2min | grep -E $eventEvalFileRegex | tee -a $configFile
    fi
  else
    echo "ERROR: unknown mode"
    exit 1
  fi
done

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
