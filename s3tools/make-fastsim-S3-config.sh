#!/bin/bash

###################
# TOP-LEVEL SCRIPT to automate the creation of a config file for fastsim files
# - automates S3 downloading and delphes execution
# - see also `s3tools/make-fastsim-local-config.sh` to generate a config file for
#   a specific directory of Delphes output root files
###################

# usage:
if [ $# -lt 3 ]; then
  echo """
  USAGE: $0 [energy] [local_dir] [mode] [limit(optional)] [config_file(optional)]

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

   - [config_file]  name of the config file; if not specified, the
                    config file will be in datarec/[local_dir]

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

# settings
Q2minima=( 1000 100 10 1 )
Q2max=0 # no max
genDir=datagen/$locDir/$energy
recDir=datarec/$locDir/$energy
if [ -z "$configFile" ]; then configFile=$recDir/delphes.config; fi

# download hepmc files from S3
function status { echo ""; echo "[+] $1"; }
if [ "$mode" == "d" -o "$mode" == "a" ]; then
  status "downloading files from S3..."
  for Q2min in ${Q2minima[@]}; do
    if [ $limit -gt 0 ]; then
      s3tools/generate-hepmc-list.sh $energy $Q2min | head -n$limit | s3tools/download.sh "$genDir/minQ2=$Q2min"
    else
      s3tools/generate-hepmc-list.sh $energy $Q2min |                 s3tools/download.sh "$genDir/minQ2=$Q2min"
    fi
  done
fi

# run delphes
if [ "$mode" == "f" -o "$mode" == "a" ]; then
  status "clean Delphes output directories"
  for Q2min in ${Q2minima[@]}; do
    rm -rvf $recDir/minQ2=$Q2min
  done
  status "running Delphes (one thread per Q2min)"
  function runDelphes { 
    echo "delphes log: " > $1/delphes.log
    for infile in $1/*.hepmc.gz; do 
      echo "RUN DELPHES ON $infile" >> $1/delphes.log
      deps/run_delphes.sh $infile 2>&1 >> $1/delphes.log
    done
    status "DONE running Delphes on directory $1"
  }
  for Q2min in ${Q2minima[@]}; do
    status "RUNNING DELPHES on energy=$1 Q2min=$Q2min..."
    runDelphes "$genDir/minQ2=$Q2min" & # comment out `&` if you want to run single threaded
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
  listFiles=""
  for Q2min in ${Q2minima[@]}; do
    configDir=$recDir/minQ2=$Q2min
    configFilePart=$configDir/delphes.config
    listFiles="$listFiles $configFilePart.list"
    status "make config file for $configDir"
    s3tools/make-fastsim-local-config.sh $energy $Q2min $Q2max $configDir $configFilePart
  done
  s3tools/generate-config-file.rb $configFile $energy $listFiles
fi

# output some info
if [ "$mode" == "c" -o "$mode" == "a" ]; then
  status "done building config file at:"
  echo "     $configFile"
fi
