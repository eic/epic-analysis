#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# run a given ROOT macro with valgrind

if [ $# -eq 0 ]; then
  echo "USAGE: $0 some_analysis_macro.C'(...)'"
  exit 2
fi

valgrind \
  --tool=memcheck \
  --track-origins=yes \
  --leak-check=full \
  --log-file="valgrind.log" \
  --suppressions=$ROOTSYS/etc/root/valgrind-root.supp \
  root.exe $*
