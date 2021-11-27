#!/bin/bash

#############
# reformat a config file into the current required format, which is one line per Q2 minimum
#
# Expected input format:
#   - one line per file
#   - columns: [filename] [entries] [cross-section] [Q2min]
# 
# Expected output format:
#   - one line per Q2min, sorted by Q2min
#   - columns: [Q2min] [cross-section] [filename_1] [entries_1] [filename_2] [entries_2] ... [filename_N] [entries_N]
#   - for the case where there are multiple files for a certain Q2min, only the cross section of the first one will be
#     used in the reformatted versio
#
# This is a temporary patch script, and may eventually be removed
# 
#############

if [ $# -ne 2 ]; then
  echo "USAGE: $0 [input-file] [output-file]"
  exit 2
fi
inFile=$1
outFile=$2
> $outFile
echo "reformatting $inFile -> $outFile"

# rearrange columns: [filename] [entries] [cross-section] [Q2min] -> [Q2min] [cross-section] [filename] [entries]
> $inFile.tmp
while read filename entries xsec Q2min; do
  echo $Q2min $xsec $filename $entries >> $inFile.tmp
done < $inFile

# sort by Q2min in decreasing order
sort -nr -o $inFile.tmp{,}

# merge lines for each value of Q2min
Q2minTmp="-1"
while read Q2min xsec filename entries; do
  if [ "$Q2min" -eq "$Q2minTmp" ]; then
    #echo "  append $filename $entries"
    lineOut="$lineOut $filename $entries"
  else
    if [ "$Q2minTmp" -ne "-1" ]; then
      echo "  -> output line for Q2min=$Q2minTmp"
      echo $lineOut >> $outFile
    fi
    Q2minTmp="$Q2min"
    #echo "  append $filename $entries"
    lineOut="$Q2min $xsec $filename $entries"
  fi
done < $inFile.tmp
echo "  -> output line for Q2min=$Q2minTmp"
echo $lineOut >> $outFile

# cleanup
rm $inFile.tmp
