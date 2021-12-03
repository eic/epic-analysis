#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for fastsim files
# - automates S3 downloading and delphes execution
# - see also `s3tools/make-fastsim-config.sh` to generate a config file for
#   a specific directory of Delphes output root files
###################

# usage:
if [ $# -lt 3 ]; then
  echo """
  USAGE: $0 [energy] [local_dir] [mode] [limit(optional)]

  provides automation for downloading hepmc files from S3, running
  Delphes on them (one thread per Q2min), and finally the generation
  of a config file

   - [energy]: 5x41 5x100 10x100 10x275 18x275

   - [local_dir]: local directory name
     - hepmc files will be downloaded to datagen/[local_dir]
     - Delphes output root files will appear in datarec/[local_dir]
     - subdirectories with energy and Q2min are created within
     - if you need to put files on a different disk, create or symlink
       datagen/[local_dir] and datarec/[local_dir] beforehand

   - [mode]:   a - run all the modes below
               d - only download hepmc files
               f - only run the delphes fast simulation
               c - only generate the config file

   - [limit]:  integer>0 : only download this many files per Q2 min
               0         : download all files (default)
               - default=5
               - limit only applies to downloading

  """
  exit 2
fi
energy=$1
locDir=$2
mode=$3
limit=5
if [ $# -ge 4 ]; then limit=$4; fi

# settings
Q2minima=( 1000 100 10 1 ) # should be decreasing order
genDir=datagen/$locDir
recDir=datarec/$locDir
configFile=$recDir/delphes.config

# download hepmc files from S3
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" -o "$mode" == "a" ]; then
  status "downloading files from S3..."
  for Q2min in ${Q2minima[@]}; do
    if [ $limit -gt 0 ]; then
      s3tools/generate-hepmc-list.sh $energy $Q2min | head -n$limit | s3tools/download.sh "$genDir/$energy/minQ2=$Q2min"
    else
      s3tools/generate-hepmc-list.sh $energy $Q2min |                 s3tools/download.sh "$genDir/$energy/minQ2=$Q2min"
    fi
  done
fi

# run delphes
if [ "$mode" == "f" -o "$mode" == "a" ]; then
  status "clean Delphes output directories"
  for Q2min in ${Q2minima[@]}; do
    rm -rv $recDir/$energy/minQ2=$Q2min
  done
  status "running Delphes (one thread per Q2min)"
  function runDelphes { 
    echo "delphes log: " > $1/delphes.log
    for infile in $1/*.hepmc.gz; do 
      echo "RUN DELPHES ON $infile" >> $1/delphes.log
      ./exeDelphes.sh $infile 2>&1 >> $1/delphes.log  # AQUI: seems to just be outputting to top-level datarec
    done
    status "DONE running Delphes on directory $1"
  }
  for Q2min in ${Q2minima[@]}; do
    status "RUNNING DELPHES on energy=$1 Q2min=$Q2min..."
    runDelphes "$genDir/$energy/minQ2=$Q2min" & # comment out `&` if you want to run single threaded
  done
  status "WAIT FOR DELPHES"
  echo "  - quit by running: while [ 1 ]; do pkill DelphesHepMC3; done"
  echo "  - monitor progress in log files in another window with:"
  echo "    find $genDir -name \"delphes.log\" | xargs tail -F"
  wait
  status "DONE RUNNING DELPHES"
fi

# generate config file
if [ "$mode" == "c" -o "$mode" == "a" ]; then
  for Q2min in ${Q2minima[@]}; do
    status "make config file for $recDir/$energy/minQ2=$Q2min"
    s3tools/make-fastsim-config.sh $energy $Q2min $recDir/$energy/minQ2=$Q2min{,/delphes.config}
  done
  status "concatenate config files to target config file: $configFile"
  > $configFile
  for Q2min in ${Q2minima[@]}; do
    cat $recDir/$energy/minQ2=$Q2min/delphes.config >> $configFile
  done
fi

# output some info
status "done building config file at:"
echo "     $configFile"
status "run root macros with parameters:"
echo "     '(\"$configFile\",$(echo $energy|sed 's/x/,/'))'"
echo ""
if [ -n "$(grep UNKNOWN $configFile)" ]; then
  >&2 echo "ERROR: missing some cross sections"
  exit 1
fi
