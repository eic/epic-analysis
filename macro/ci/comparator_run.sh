#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# wrapper for CI usage of comparator.C

set -e

echo "CALLED $0"
if [ $# -ne 14 ]; then
  echo "ERROR: improper args (see script)"
  exit 1
fi

aname=$1
pname=$2
recon=$3
xvar=$4
yvar=$5
title1=$6;    name1=$7
title2=$8;    name2=$9
title3=${10}; name3=${11}
title4=${12}; name4=${13}
outname=${14}

dets="$name1 $name2 $name3 $name4"

# rename aname -> pname
for det in $dets; do
  cp out/$det.{$aname,$pname}.$recon.root
done

# run comparator.C
args=""
args+="\"$title1\",\"out/$name1.$pname.$recon.root\","
args+="\"$title2\",\"out/$name2.$pname.$recon.root\","
args+="\"$title3\",\"out/$name3.$pname.$recon.root\","
args+="\"$title4\",\"out/$name4.$pname.$recon.root\","
args+="\"out/comparison.$outname.$pname.$recon\",\"$xvar\",\"$yvar\""
export CLING_DEBUG=1
root -b -q "macro/ci/comparator.C($args)"

echo "DONE $0"
