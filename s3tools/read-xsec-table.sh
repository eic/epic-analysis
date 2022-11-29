#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks

# read the cross section table, and find the associated cross section

# arguments
table=datarec/xsec/xsec.dat
if [ $# -lt 1 ]; then echo """
  USAGE: $0 [search_string] [xsec.dat file (optional)]
  EXAMPLE: $0 18x275/minQ2=100
  default dat file: $table
  """
  exit 2
fi
if [ $# -ge 2 ]; then table=$2; fi

# grep for query ("label")
query="$1" 
xsec=$(grep -E "^${query} " $table | awk '{print $2}')

# not found error
if [ -z "$xsec" ]; then
  >&2 echo "ERROR: cannot find cross section for $query in $table"
  xsec="UNKNOWN"
fi

# return
echo $xsec
