#!/bin/bash

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# update or download and install delphes

### download/update source
source environ.sh

echo "[+] downloading onnx runtime source"
wget https://github.com/microsoft/onnxruntime/releases/download/v1.14.1/onnxruntime-linux-x64-1.14.1.tgz
mv onnxruntime-linux-x64-1.14.1.tgz deps
tar -zxvf deps/onnxruntime-linux-x64-1.14.1.tgz

