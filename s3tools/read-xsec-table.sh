#!/bin/bash
# read the cross section table, and find the associated cross section

# arguments
table=datarec/xsec/xsec.dat
if [ $# -lt 2 ]; then echo """
  USAGE: $0 [energy] [Q2min] [xsec.dat file (optional)]

  EXAMPLE: $0 18x275 100
  default dat file: $table
  """
  exit 2
fi
if [ $# -ge 3 ]; then table=$3; fi
query="$1/minQ2=$2" 
grep -w $query $table | awk '{print $2}'
