#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks, Connor Pecar

require 'fileutils'

# arguments
if ARGV.length < 3
  puts """
  USAGE:
    #{$0} [analysis_macro] [directory of config.part files] [output subdirectory] [additional macro args]...

    - [analysis_macro]: the analysis ROOT macro to run on each file
                        NOTE: the first two arguments of the macro must be:
                              - input config file
                              - output file prefix

    - [directory of config.part files]: directory containing config.part files to run

    - [output subdirectory]: output subdirectory for resulting ROOT files
                             NOTE: this will be a subdirectory of out/
                             NOTE: the ROOT files in this directory will be REMOVED

    - [additional macro args]: all remaining arguments will be sent to the macro
                               NOTE: prepend strings with an underscore (_)
                                     to ensure they will have quotation marks in
                                     the macro call
  """
  exit 2
end
RootMacro, InDir, OutSubDir = ARGV[0..2]
AdditionalArgs              = ARGV[3..-1]

# set directories
OutDir         = "out/#{OutSubDir}"
ShellScriptDir = "#{OutDir}/scripts"
LogDir         = "hpc/log"
LogSubDir      = "#{LogDir}/#{Time.now.to_i.to_s}" # use unix time for log subdirectory

# ask a yes/no question to the user
def ask(question)
  print "#{question}\n> "
  answer = $stdin.gets.chomp.split('').first
  return false if answer.nil?
  answer.downcase == "y"
end

# make/clean directories
cleanDirs = [
  "#{OutDir}/parts",
  ShellScriptDir,
]
puts "\nCleanup ...\nRemoving the following directories:"
cleanDirs.each{ |dir| puts "- #{dir}" }
correct = ask("\nIs this OK? [y/N]")
puts "you answered " + (correct ? "yes" : "no; stopping!")
exit unless correct
cleanDirs.each do |dir|
  FileUtils.rm_rf   dir, verbose: true
  FileUtils.mkdir_p dir, verbose: true
end
FileUtils.mkdir_p LogSubDir, verbose: true

# glob config.part files
partFileList = Dir.glob "#{InDir}/*.config.part"
sep = '-'*50
puts 'config.part files:', sep, partFileList, sep

# add explicit quotes around string arguments
def enquote(str)
  "\"#{str}\""
end
RootMacroArgs = AdditionalArgs.map do |arg|
  if arg.match? /^_/
    enquote arg.sub(/^_/,'')
  else
    arg
  end
end

# locate eic-shell # FIXME: this may not be correct for everyone!
eicShellPrefix = ENV['EIC_SHELL_PREFIX']
if eicShellPrefix.nil?
  $stderr.puts "ERROR: unknown $EIC_SHELL_PREFIX, cannot find eic-shell"
  exit 1
end
eicShell = eicShellPrefix.split('/')[0..-2].join('/') + '/eic-shell'

# generate slurm config
slurmConfigN = "#{OutDir}/scripts/run.slurm"
commandListFile = "#{OutDir}/scripts/commandlist.slurm"
File.open(slurmConfigN, 'w') do |slurmConfig|
  File.open(commandListFile, 'w') do |commandList|
    # header
    # Split the InDir string and extract the desired text
    desired_text = InDir.split("___")[1].gsub('/', '')
    slurmConfig.puts """#!/bin/bash
#SBATCH --job-name=epic-analysis-#{desired_text}
#SBATCH --account=eic
#SBATCH --partition=production
#SBATCH --mem-per-cpu=500
#SBATCH --time=24:00:00"""
    
    # job array
    nfiles = partFileList.count()
    slurmConfig.puts """#SBATCH --array=1-#{nfiles}"""                  
    # output
    slurmConfig.puts """#SBATCH --output=#{LogSubDir}/%x-%j-%N.out
#SBATCH --error=#{LogSubDir}/%x-%j-%N.err
    """
    
    slurmConfig.puts """
srun $(sed -n ${SLURM_ARRAY_TASK_ID}p #{commandListFile})
    """
                       
    # loop over config.part files
    partFileList.each do |partFile|

      # get root file name from this config.part file
      rootFile = File.open(partFile)
                   .readlines
                   .map(&:chomp)
                   .grep_v(/^:/)
                   .grep(/\.root$/)
                   .first
      if rootFile.nil?
        $sterr.puts "ERROR: cannot find ROOT file in #{partFile}"
        next
      end

      # generate output file prefix
      outFilePrefix = [
        "#{OutSubDir}/parts/",
        File.basename(partFile,'.config.part'),
        '__',
        File.basename(rootFile,'.root'),
      ].join

      # set ROOT macro arguments
      macroArgs = [
        enquote(partFile),
        enquote(outFilePrefix),
        *RootMacroArgs,
      ].join(',')

      # set log file name
      logName = LogSubDir + '/' + File.basename(outFilePrefix)

      # generate wrapper shell script; this is needed so we can run from outside eic-shell
      shellScriptName = ShellScriptDir + '/make__' + File.basename(outFilePrefix) + '.sh'
      File.open(shellScriptName,'w') do |script|
        script.puts '#!/bin/bash'
        script.puts """
# call as: eic-shell -- #{shellScriptName}
source environ.sh
root -b -q #{RootMacro}'(#{macroArgs})'"""
      end
      FileUtils.chmod 'u+x', shellScriptName
      commandList.puts """#{eicShellPrefix}/../eic-shell -- #{shellScriptName}
"""
    end # loop over config.part files
  end # close slurm command list file
end # close slurmConfig file

# print some things to check
puts """
#{sep}

Sample shell script:
#{sep}
#{`cat #{Dir.glob("#{ShellScriptDir}/*.sh").first}`}
#{sep}
Assuming host eic-shell is: #{eicShell}
... if this is incorrect, fix #{$0} (help needed) ...
"""

# print out what to do next
puts """
#{sep}
Produced slurm config file: #{slurmConfigN}

Make sure everything is correct, then submit jobs by running
the following command from OUTSIDE of eic-shell:

  sbatch #{slurmConfigN}

Log files will be written to: #{LogSubDir}/
Monitor errors with:

  more #{LogSubDir}/*.err

When slurm jobs are done, merge the ROOT files by running
(from inside eic-shell):

  hpc/merge.rb #{OutDir}

#{sep}
"""
