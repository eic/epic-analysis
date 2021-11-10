#!/bin/bash
# generate config file for a local directory of simulation files

# arguments and usage
if [ $# -eq 0 ]; then echo """
USAGE: $0 [directory] [columns]
- [directory] is a local directory
- [columns] are for additional columns you want to append
- output is to stdout, so you can either:
  - redirect it to a file
  - pipe through grep for further filtering
  - you may need to edit the grep filter in the script
"""; exit 2; fi
dataDir=$1
shift
columns=$@

# get file list from local directory ########################
ls $dataDir |\
while read fileName; do
  echo "$dataDir/$fileName $columns"
done
