#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# add host to MinIO client config
# - requires env vars S3_SECRET_KEY and S3_ACCESS_KEY
if [ -z "$S3_SECRET_KEY" -o -z "$S3_ACCESS_KEY" ]; then
  echo "ERROR: need to set env vars S3_SECRET_KEY and S3_ACCESS_KEY"
  exit 1
fi
hostURL="https://eics3.sdcc.bnl.gov:9000"
mc config host add S3 $hostURL $S3_ACCESS_KEY $S3_SECRET_KEY
