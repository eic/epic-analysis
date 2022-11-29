#!/bin/bash
# update or download and install delphes

### download/update source
source environ.sh
echo "[+] downloading pybind11 source"
git submodule add -b stable ../../pybind/pybind11 deps/pybind11
git submodule update --init
