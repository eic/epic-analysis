#!/bin/bash
# generate config file for a directory on S3

# arguments and usage
if [ $# -eq 0 ]; then echo """
USAGE: $0 [directory] [columns]
- [directory] is a directory on S3
- [columns] are for additional columns you want to append;
  currently we use columns [q2min] [crossSection]
- output is to stdout, so you can either:
  - redirect it to a file
  - pipe through grep for further filtering
  - you may need to edit the grep filter in the script
"""; exit 2; fi
remoteDir=$1
shift
columns=$@

# format URLs
hostURL="https://dtn01.sdcc.bnl.gov:9000"
hostURL_esc=$(echo $hostURL | sed 's/\//\\&/g')
remoteDirURL=$(echo $remoteDir | sed "s/^S3/s3$hostURL_esc/" | sed 's/\//\\&/g')

# check environment
if [ -z "$S3_SECRET_KEY" -o -z "$S3_ACCESS_KEY" ]; then
  echo "NOTE: consider setting env vars S3_SECRET_KEY and S3_ACCESS_KEY"
else
  mc config host add S3 $hostURL $S3_ACCESS_KEY $S3_SECRET_KEY
fi

# get file list from S3 ###################################
mc ls $remoteDir | sed 's/^.*] //' | awk '{print $2}' |\
  grep -vE 'raw.root$|ecal.root$|hcal.root$' |\
  sed "s/$/ $columns/" |\
  sed "s/^/ $remoteDirURL\//"
