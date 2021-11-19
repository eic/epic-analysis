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

# grep for query ("label")
query="$1/minQ2=$2" 
xsec=$(grep -w $query $table | awk '{print $2}')

# not found error
if [ -z "$xsec" ]; then
  >&2 echo "ERROR: cannot find cross section for $query in $table"
  xsec="UNKNOWN"
fi

# return
echo $xsec
