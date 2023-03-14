#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

require 'fileutils'

# arguments
if ARGV.length != 2
  $stderr.puts """
  USAGE: #{$0} [config_file] [directory_name]
    - [config_file] will be split into several files, one
      config file per ROOT file; the new files will be
      written to [directory_name]
  """
  exit 2
end
ConfigFileN, OutDirN = ARGV

# clean output directory
puts "\nCleaning #{OutDirN} ..."
FileUtils.mkdir_p OutDirN
FileUtils.rm_f Dir.glob("#{OutDirN}/*config.part"), verbose: true


# create config.template file, which we will use to generate config.part files
# - the template file is a copy of the essential parts of the config file
# - special strings that match the regex /__.*__/ will be modified when generating
#   config.part files from the template
templateFileN   = ConfigFileN + '.template'
templateFile    = File.open templateFileN, 'w'
configFileHash  = Hash.new
numEventsHash   = Hash.new
readingFileList = false
rootFileNum     = 0
q2key           = Proc.new{ configFileHash[':q2min'] + '_' + configFileHash[':q2max'] }
puts "\nParsing #{ConfigFileN} and counting the number of events ..."
File.open(ConfigFileN).readlines.each do |line_in|

  # remove comments and newlines
  line = line_in.gsub(/#.*/,'').chomp

  # parse line: if key-value pair
  if line.match? /^:/

    # if already reading a list of ROOT files, reset some settings for this Q2 range
    if readingFileList or configFileHash.size==0
      configFileHash[':q2min']        = '1.0'
      configFileHash[':q2max']        = '0.0'
      configFileHash[':crossSection'] = '0.0'
      templateFile.puts ":endGroup\n\n" if readingFileList
      readingFileList = false
    end

    # add this key-value pair to the hash
    key = line.split.first
    val = line.split[1..-1].join ' '
    configFileHash[key] = val

    # add to config.template file
    templateFile.puts line

  # parse line: if a ROOT file, produce a new config.part file
  elsif line.match? /\.root/
    rootFile = line.split.first

    # get the number of events in this ROOT file
    if numEventsHash[q2key.call].nil?
      numEventsHash[q2key.call] = {
        :numEvents => 0,
      }
    end
    puts "Counting events in #{rootFile}"
    numEvents = `hpc/src/count_events.exe #{rootFile}`.chomp.to_i
    numEventsHash[q2key.call][:numEvents] += numEvents
    puts "  => #{numEvents} events"

    # add template line for numEvents, if reading a new list of ROOT files
    unless readingFileList
      templateFile.puts "__numEvents__ #{q2key.call}"
    end

    # add ROOT file to template and prepare to parse the next line
    templateFile.puts "__ROOT_#{rootFileNum}__ #{rootFile}"
    readingFileList = true
    rootFileNum += 1

  end # end parsing
end # end loop over config file lines
templateFile.puts ":endGroup"
templateFile.close


# generate config.part files from the template
puts "\nParsing #{templateFileN} to generate config.part files ..."
rootFileNum.times do |i|

  # start new config.part file
  partFileN = OutDirN + '/' + File.basename(ConfigFileN).sub(/\.config$/,'') + ".%07d.config.part" % [i]
  rootFile = 'ERROR: UNKNOWN'
  File.open(partFileN,'w') do |partFile|

    # parse template, and write to config.part file
    File.open(templateFileN).readlines.each do |line|
      if line.match? /^__numEvents__/
        numEvents = numEventsHash[line.split.last][:numEvents]
        line = ":numEvents #{numEvents}"
      elsif line.match? /^__ROOT_#{i}__/
        rootFile = line.split.last
        line = rootFile
      elsif line.match? /^__ROOT/
        next
      end
      partFile.puts line
    end

  end
  puts "#{partFileN} => #{rootFile}"

end
puts "\nDone. Config files written to #{OutDirN}/"
