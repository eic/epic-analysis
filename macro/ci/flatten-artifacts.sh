#!/bin/bash
# flatten artifact subdirectories (`mv *.images/* ./`)
if [ $# -ne 1 ]; then echo "USAGE: $0 [artifacts-dir]"; exit 2; fi
pushd $1
ls -d */ | while read d; do
  echo "FLATTEN $d"
  pushd $d
  if [ -n "$(ls -d */ | grep '\.images')" ]; then
    mv -v *.images/* ./
    rm -r *.images
  fi
  popd
done
popd
