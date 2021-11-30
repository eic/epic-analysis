#!/bin/bash
# organize artifacts by combining directories
# - execute in CI pipelines from top-level directory, after all artifacts
#   have been collected into one directory, specified by argument

# arguments
if [ $# -ne 1 ]; then echo "USAGE: $0 [artifacts-dir]"; exit 2; fi
pushd $1

##########################
# merge directories of fastsim and fullsim files, for comparison
# - adds `.fastsim` or `.fullsim` suffixes to filenames
##########################
echo "------------- merge fastsim and fullsim -------------------"
# loop through *fastsim* directories
ls | grep fastsim |\
while read dirFast; do

  # set fullsim and output directory names
  echo "--------------------"
  dirFull=$(echo $dirFast | sed 's/fastsim/fullsim/g')
  dirOut=$(echo $dirFast | sed 's/fastsim//g' | sed 's/--/-/g')
  mkdir -p $dirOut

  # move fastsim artifacts to output directory, and add suffix to file name
  echo "MOVE AND RENAME artifacts in $dirFast -> $dirOut"
  pushd $dirFast
  mv -v *.images/* ./; rm -r *.images
  for file in *; do
    mv -v $file ../$dirOut/$(echo $file | sed 's/^.*\./&fastsim./g')
  done
  popd

  # repeat for fullsim artifacts
  echo "MOVE AND RENAME artifacts in $dirFull -> $dirOut"
  pushd $dirFull
  mv -v *.images/* ./; rm -r *.images
  for file in *; do
    mv -v $file ../$dirOut/$(echo $file | sed 's/^.*\./&fullsim./g')
  done
  popd

  # remove empty directories
  echo "CLEANUP"
  rm -rv $dirFast $dirFull

done


########################
# merge recon-* directories, adding recon method to suffix
########################
echo "------------- merge recon methods -------------------"
# loop through recon* directories
ls | grep -E '^recon-' |\
while read dirRecon; do
  method=$(echo $dirRecon | sed 's/-plots$//g' | sed 's/^.*-//g')
  dirOut=$(echo $dirRecon | sed "s/-$method//g" | sed 's/^recon-/&methods-/g')
  echo "MOVE ARTIFACTS IN $dirRecon/ FOR METHOD \"$method\" TO $dirOut/"
  mkdir -p $dirOut
  pushd $dirRecon
  mv -v *.images/* ./; rm -r *.images
  for file in *; do
    mv -v $file ../$dirOut/$(echo $file | sed "s/^.*\./&$method./g")
  done
  popd
  rm -r $dirRecon
done


############
popd
