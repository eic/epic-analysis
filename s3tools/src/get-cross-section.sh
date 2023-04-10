#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# get generated cross section from hepmc file

# arguments
limiter=10000
if [ $# -eq 0 ]; then echo """
  USAGE: $0 [hepmc file] [event num limiter (default=$limiter)]

  - generated cross section value is asymptotic; use the last event's value
  - if reading a file on S3, use the [event num limiter] argument to specify which event
    number to stop at; this is needed since 'tail' does not exist in MinIO
    client, which means hepmc files can only be streamed from first event to last
    event
  - if reading a local file, 'tail' is used instead, and [event num limiter] is ignored
  - the last ~50 events' cross sections are printed, to give an idea of local convergence
  - the columns after \"GenCrossSection\" are: 
    - double cross_section; // cross section in pb.
    - double cross_section_error; // error associated with this cross section.
    - long accepted_events; ///< The number of events generated so far.
    - long attempted_events; ///< The number of events attempted so far.

"""
  exit 2
fi
if [ $# -ge 2 ]; then limiter=$2; fi
hepmcFile=$1

# grep for cross section | tail
if [[ $hepmcFile =~ ^S3\/.* ]]; then # S3 file
  mc cat $hepmcFile | grep GenCrossSection | head -n$limiter | tail -n50
else # local file
  tail -n1000 $hepmcFile | grep GenCrossSection | tail -n50
fi
