#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

require 'fileutils'

# arguments
if ARGV.length != 3
  $stderr.puts """
  USAGE: #{$0} [config_file] [directory_name] [root_files_per_config]
    - [config_file] will be split into several files, one
      config file per specified number of ROOT files; the new files will be
      written to [directory_name]
    - [root_files_per_config] specifies how many ROOT files to include in each config file
  """
  exit 2
end
ConfigFileN, OutDirN, RootFilesPerConfig = ARGV

# Validate RootFilesPerConfig
RootFilesPerConfig = RootFilesPerConfig.to_i
if RootFilesPerConfig <= 0
  $stderr.puts "Error: root_files_per_config must be a positive integer"
  exit 3
end

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

part_file_index = 0
root_files_accumulated = 0
rootFile = nil  # Declare rootFile here
created_files = []  # Array to store the names of created files

rootFileNum.times do |i|
  if root_files_accumulated == 0
    # Start a new config.part file
    @partFileN = OutDirN + '/' + File.basename(ConfigFileN).sub(/\.config$/, '') + ".%07d.config.part" % [part_file_index]
    part_file_index += 1
    @partFile = File.open(@partFileN, 'w')
  end

  # Increment the count of accumulated ROOT files
  root_files_accumulated += 1

  # Process the template file
  File.open(templateFileN).each do |line|
    if line.match?(/^__numEvents__/)
      numEvents = numEventsHash[line.split.last][:numEvents]
      line = ":numEvents #{numEvents}"
    elsif line.match?(/^__ROOT_#{i}__/)
      rootFile = line.split.last  # Update rootFile here
      line = rootFile
    elsif line.match?(/^__ROOT/)
      next
    end
    @partFile.puts line
  end

  # Check if we have accumulated enough ROOT files or if this is the last ROOT file
  if root_files_accumulated == RootFilesPerConfig || i == rootFileNum - 1
    @partFile.close
    puts "#{@partFileN} => #{rootFile}"  # rootFile is now defined
    created_files << @partFileN  # Store the file name
    root_files_accumulated = 0
  end
end

puts "\nDone. Config files written to #{OutDirN}/"

# The ROOT files created above may have repeated Q2 blocks 
# Refactor the saved config files with the correct Q2 ordering

created_files.each do |file_name|
  # Open the existing file for reading
  existing_file = File.open(file_name, 'r')
  grouped_lines = Hash.new { |hash, key| hash[key] = { cross_section: nil, lines: [] } }

  # Variables to store common values
  eleBeamEn = ionBeamEn = crossingAngle = totalCrossSection = nil

  # Read the file and group lines by q2min, also capturing common values
  existing_file.each do |line|
    case line
    when /^:eleBeamEn\s+(\d+)/
      eleBeamEn = $1
    when /^:ionBeamEn\s+(\d+)/
      ionBeamEn = $1
    when /^:crossingAngle\s+([-\d]+)/
      crossingAngle = $1
    when /^:totalCrossSection\s+(\d+\.\d+)/
      totalCrossSection = $1
    when /^:crossSection\s+(\d+\.\d+)/
      current_cross_section = $1
      grouped_lines[@current_q2min][:cross_section] = current_cross_section if current_cross_section
    when /^:q2min\s+(\d+\.\d+)/
      @current_q2min = line.strip
    when /^s3/
      grouped_lines[@current_q2min][:lines] << line.strip if @current_q2min
    end
  end

  existing_file.close

  # Write the new file
  File.open(file_name, 'w') do |file|
    # Write the dynamic lines at the top
    file.puts(":eleBeamEn #{eleBeamEn}")
    file.puts(":ionBeamEn #{ionBeamEn}")
    file.puts(":crossingAngle #{crossingAngle}")
    file.puts(":totalCrossSection #{totalCrossSection}")
    file.puts("")

    # Iterate over grouped lines and write them
    # To read s3 files in slurm, you need to read them from xrootd
    # See https://eic.github.io/epic-prod/documentation/faq.html
    grouped_lines.each do |q2min, data|
      file.puts("# Q2 range")
      file.puts(q2min)
      file.puts(":crossSection #{data[:cross_section]}")
      data[:lines].each do |line|
        modified_line = line.gsub('s3https://eics3.sdcc.bnl.gov:9000/eictest/', 'root://dtn-eic.jlab.org//work/eic2/')
        file.puts(modified_line)
      end
      file.puts("")
    end
  end

  puts "Refactored config file: #{file_name}"
end
