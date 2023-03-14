#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

# arguments
outputBaseName = 'analysis'
if ARGV.length < 1
  puts """
  USAGE:
    #{$0} [output directory] [output basename (optional)]

    - [output directory]: the output directory containing ROOT part files,
                          which are assumed to be at: 

                            [output directory]/parts/*.root

    - [output basename]: optionally specify an output basename, so the output
                         file name will be:

                            [output directory]/[output basename].root
                            
                            default [output basename]: #{outputBaseName}

    NOTE: for more control, run: hpc/src/merge_analysis_files.exe
  """
  exit 2
end
outputDir      = ARGV[0]
outputBaseName = ARGV[1] if ARGV.length>1

# run merge_analysis_files.exe
partFileList   = Dir.glob "#{outputDir}/parts/*.root"  # get list of input part files
outputFileName = "#{outputDir}/#{outputBaseName}.root" # set output file name
system "hpc/src/merge_analysis_files.exe #{outputFileName} #{partFileList.join ' '}" # merge
