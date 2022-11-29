#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks

# download list of files from S3

# arguments and usage
if [ $# -lt 1 ]; then
  echo "USAGE: $0 path/to/target/directory file1 file2 ..."
  exit 2
fi
targetPath=$1
mkdir -p $targetPath
shift

# start downloader script
downloadScript=$targetPath/get-files.sh
echo "#!/bin/bash" > $downloadScript

# read list of source files and generate downloader script
while read sourceURL; do
  sourcePath=$(echo $sourceURL | awk '{print $1}' | sed 's/https:\/\///g' | sed 's/^[^\/]*\//S3\//g')
  echo "mc cp $sourcePath $targetPath; echo \"\"" >> $downloadScript
done
echo "generated downloader script $downloadScript"

# download
echo "now starting download..."
chmod u+x $downloadScript
$downloadScript
echo "done."
