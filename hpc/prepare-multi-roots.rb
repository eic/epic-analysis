#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

require 'csv'
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

rootFile = nil  # Declare rootFile here
created_files = []  # Array to store the names of created files

template_line_data = []
current_q2_min = ""
current_q2_xsec = ""
eleBeamEn = ""
ionBeamEn = ""
crossingAngle = ""

# Parse through the templateFile
# While getting the numEvents for each ROOT file, check with a locally stored
# database for improved speed. If the numEvents is not found, add it to the database

# File.open(templateFileN, 'r') do |templateFile|
#   templateFile.each do |line|
#     line.strip!
#     next if line.empty?
    
    
#     # Capture the values of :eleBeamEn, :ionBeamEn, and :crossingAngle
#     if line.start_with?(":eleBeamEn")
#       eleBeamEn = line
#     elsif line.start_with?(":ionBeamEn")
#       ionBeamEn = line
#     elsif line.start_with?(":crossingAngle")
#       crossingAngle = line
#     elsif line.start_with?(":q2min")
#       current_q2_min = line
#     elsif line.start_with?(":crossSection")
#       current_q2_xsec = line
#     elsif line.include?("__ROOT_")
#       rootFileN = line.split(" ")[1]
#       puts "Counting events in #{rootFileN}"
#       numEvents = `hpc/src/count_events.exe #{rootFileN}`.chomp.to_i
#       puts "  => #{numEvents} events"
#       template_line_data << [current_q2_min, current_q2_xsec, numEvents, rootFileN]
#     end
#   end
# end

def create_directory_structure(base_path, subdirs)
  dir_path = File.join(base_path, *subdirs)
  FileUtils.mkdir_p(dir_path) unless Dir.exist?(dir_path)
  dir_path
end

def file_exists_in_csv?(filename, csv_file)
  return false unless File.exist?(csv_file)
  CSV.foreach(csv_file) do |row|
    return true if row[0] == filename
  end
  false
end

def get_num_events_from_csv(filename, csv_file)
  CSV.foreach(csv_file) do |row|
    return row[1].to_i if row[0] == filename
  end
end

base_path = 'hpc/nevents_databases'

File.open(templateFileN, 'r') do |templateFile|
  templateFile.each do |line|
    line.strip!
    next if line.empty?
    # Capture the values of :eleBeamEn, :ionBeamEn, and :crossingAngle
    if line.start_with?(":eleBeamEn")
      eleBeamEn = line
    elsif line.start_with?(":ionBeamEn")
      ionBeamEn = line
    elsif line.start_with?(":crossingAngle")
      crossingAngle = line
    elsif line.start_with?(":q2min")
      current_q2_min = line
    elsif line.start_with?(":crossSection")
      current_q2_xsec = line
    elsif line.include?("__ROOT_")
      rootFileN = line.split(" ")[1]
      basename = File.basename(rootFileN)

      # Extract directory names from rootFileN
      path_parts = rootFileN.split('/').slice(5, 7) # Adjust indices as needed
      dir_path = create_directory_structure(base_path, path_parts)
      csv_file = File.join(dir_path, 'data.csv')

      # Check if rootFileN exists in the CSV
      if file_exists_in_csv?(rootFileN, csv_file)
        numEvents = get_num_events_from_csv(rootFileN, csv_file)
        puts "Retrieved #{numEvents} events for #{basename} from database"
      else
        puts "Counting events in #{basename}"
        modified_rootFileN = rootFileN.gsub('s3https://eics3.sdcc.bnl.gov:9000/eictest/', 'root://dtn-eic.jlab.org//work/eic2/')
        numEvents = `hpc/src/count_events.exe #{modified_rootFileN}`.chomp.to_i
        puts "  => #{numEvents} events"

        # Append to CSV in the corresponding directory
        CSV.open(csv_file, 'a') do |csv|
          csv << [rootFileN, numEvents]
        end
      end

      template_line_data << [current_q2_min, current_q2_xsec, numEvents, rootFileN]
    end
  end
end

# Calculate the Weights for the Monte Carlo
# These will be manually inserted into the config files
q2bank = {}
max_q2_xsec_value = 0

template_line_data.each do |data|
  current_q2_min, current_q2_xsec, numEvents, _ = data

  # Numerical values for the q2min and q2xsec
  q2_min_value = current_q2_min.split(" ")[1].to_f
  q2_xsec_value = current_q2_xsec.split(" ")[1].to_f
  max_q2_xsec_value = q2_xsec_value if q2_xsec_value>max_q2_xsec_value
    
  unless q2bank.has_key?(q2_min_value.to_s)
    q2bank[q2_min_value.to_s] = {
      "xsec" => q2_xsec_value,
      "numEvents" => 0,
      "Weight" => 0
    }
  end
  q2bank[q2_min_value.to_s]["numEvents"] += numEvents
end

entriesTotal = 0
q2bank.each do |key, value|
  entriesTotal += value["numEvents"]
end

lumiTotal = entriesTotal / max_q2_xsec_value

q2bank.each do |keyi, valuei|
  q2min_i = keyi.to_f
  q2max_i = 1000000 # Default, currently only analyzing q2 value with only specified minimum
  lumiThis=0
  # check if Q2 range `j` contains the Q2 range `i`; if so, include its luminosity
  q2bank.each do |keyj, valuej|
     q2min_j = keyj.to_f
     q2max_j = q2max_i
     if q2min_i >= q2min_j and q2min_i<=q2max_j and q2max_i >= q2min_j and q2max_i <= q2max_j
         lumiThis += valuej["numEvents"] / valuej["xsec"]
     end
  end
  valuei["Weight"] = lumiTotal / lumiThis
  valuei["Weight"] = sprintf("%.10f", valuei["Weight"]) # reformat to remove things like 5E-5

end


i = 0
root_files_accumulated = 0
part_file_index = 0
created_files = []
# Generate config files for each small batch

loop do
  data = template_line_data[i]
  current_q2_min, current_q2_xsec, numEvents, _ = data

  # Numerical values for the q2min and q2xsec
  q2_min_value = current_q2_min.split(" ")[1].to_f
  q2_xsec_value = current_q2_xsec.split(" ")[1].to_f
  q2_weight = q2bank[q2_min_value.to_s]["Weight"]

  if root_files_accumulated.zero?
    @partFileN = OutDirN + '/' + File.basename(ConfigFileN).sub(/\.config$/, '') + ".%07d.config.part" % [part_file_index]
    @partFile = File.open(@partFileN, 'w')
    @partFile.puts eleBeamEn
    @partFile.puts ionBeamEn
    @partFile.puts crossingAngle
    @partFile.puts ":totalCrossSection ???"
    part_file_index += 1
  end
  i += 1
  root_files_accumulated += 1

  @partFile.puts data[0]
  @partFile.puts data[1]
  @partFile.puts ":Weight #{q2_weight}"
  @partFile.puts data[3]

  @partFile.puts "\n"

  if root_files_accumulated == RootFilesPerConfig || i == rootFileNum
    # Because of Analysis::GetEventQ2Idx() in Analysis.cxx, we need to
    # include the Q2weights of all files in the simulation, even if the batch only
    # contains 1 range
    q2bank.each do |key, value|
      @partFile.puts ":q2min #{key}"
      @partFile.puts ":crossSection #{value['xsec']}"
      @partFile.puts ":Weight #{value['Weight']}"
      @partFile.puts "\n"
    end
    # Before closing the file, replace ??? with max_q2_xsec_value
    @partFile.close
    @partFile = File.open(@partFileN, 'r+')
    file_contents = @partFile.read  # Read the contents
    @partFile.close  # Close to reset the file pointer
    file_contents.gsub!(":totalCrossSection ???", ":totalCrossSection #{max_q2_xsec_value}")  # Replace placeholder
    File.write(@partFileN, file_contents)  # Write the updated contents back to the file

    puts "#{@partFileN} => #{data[2]}"  # rootFile is now defined
    created_files << @partFileN
    root_files_accumulated = 0
    break if i == rootFileNum
  end
end


puts "\nDone. Config files written to #{OutDirN}/"

# The ROOT files created above may have repeated Q2 blocks 
# Refactor the saved config files with the correct Q2 ordering

created_files.each do |file_name|
  # Open the existing file for reading
  existing_file = File.open(file_name, 'r')
  grouped_lines = Hash.new { |hash, key| hash[key] = { cross_section: nil, lines: [] , weight: nil} }

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
    when /^:Weight\s+(\d+\.\d+)/
      current_weight = $1
      grouped_lines[@current_q2min][:weight] = current_weight if current_weight
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
    # Sort grouped lines by extracting and converting q2min to a floating-point number
    sorted_grouped_lines = grouped_lines.sort_by do |q2min, _|
      # Extract the numerical part of the q2min string and convert to float
      q2min.match(/:q2min\s+(\d+\.\d+)/)[1].to_f
    end
    
    # Iterate over sorted grouped lines and write them
    sorted_grouped_lines.each do |q2min, data|
      file.puts("# Q2 range")
      file.puts(q2min)
      file.puts(":crossSection #{data[:cross_section]}")
      file.puts(":Weight #{data[:weight]}")
      data[:lines].each do |line|
        modified_line = line.gsub('s3https://eics3.sdcc.bnl.gov:9000/eictest/', 'root://dtn-eic.jlab.org//work/eic2/')
        file.puts(modified_line)
      end
      file.puts(":endGroup") # Forces AddFileGroup to be run in Analysis.cxx, even if FileGroup is empty
      file.puts("")
    end
  end

  puts "Refactored config file: #{file_name}"
end
