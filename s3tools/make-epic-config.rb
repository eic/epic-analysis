#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks

require 'optparse'
require 'ostruct'
require 'fileutils'

# default CLI options
options = OpenStruct.new
options.version    = 'epic.22.11.3'
options.energy     = '18x275'
options.locDir     = ''
options.mode       = 's'
options.limit      = 2
options.configFile = ''
options.detector   = 'arches'
options.radCor     = false

# global settings
CrossSectionTable = 'datarec/xsec/xsec.dat'
HostURL           = 'https://eics3.sdcc.bnl.gov:9000'


# production specifications, latest first
# PRODUCTION_VERSION => {
#   :comment         => Description about this version
#   :crossSectionID  => Proc(minQ2,maxQ2,radDir) -> label of row of `CrossSectionTable`
#   :releaseSubDir   => Proc() -> Directory for this PRODUCTION_VERSION
#   :energySubDir    => Proc() -> Subdirectory associated to user-specified beam energy
#   :dataSubDir      => Proc(*version dependent*) -> Subdirectory of `:energySubDir`
# }
prodSettings = {
  'epic.22.11.3' => {
    :comment         => 'Pythia 6, with & without radiative corrections',
    :crossSectionID  => Proc.new { |minQ2,maxQ2,radDir| "pythia6:ep_#{radDir}.#{options.energy}_q2_#{minQ2}_#{maxQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{options.version.sub('epic.','')}/epic_#{options.detector}/SIDIS/pythia6" },
    :energySubDir    => Proc.new { "ep_#{options.energy}" },
    :dataSubDir      => Proc.new { |radDir|
      if [options.energy,radDir]==['18x275','noradcor']  # correct for S3 disorganization
        "hepmc_ip6"
      else
        "hepmc_ip6/#{radDir}"
      end
    },
  },
  'epic.22.11.2' => {
    :comment         => 'Pythia 8',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{options.version.sub('epic.','')}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
}

# additional lists
typicalEnergyList = [
  "5x41",
  "5x100",
  "10x100",
  "10x275",
  "18x275",
]
detectorConfigList = [
  "arches",
  "brycecanyon",
]

# parse options
OptionParser.new do |o|
  o.banner = "USAGE: #{$0} [OPTIONS]..."
  o.separator ''
  o.separator 'OPTIONS:'
  o.on("-v", "--version [PRODUCTION_VERSION]",
       "Production campaign version, one of:",
       *prodSettings.keys.map{|e|"  #{e}"},
       "Default: #{options.version}"
      ) do |a|
        unless prodSettings.keys.include? a
          $stderr.puts "ERROR: unknown DETECTOR_VERSION '#{a}'"
          exit 1
        end
        options.version = a
      end
  o.separator ''
  o.on("-e", "--energy [ENERGY]",
       "Energy setting, one of:",
       *typicalEnergyList.map{|e|"  #{e}"},
       "Default: #{options.energy}",
      ) do |a|
        $stderr.puts "WARNING: '#{a}' is not a typical beam energy" unless typicalEnergyList.include? a
        options.energy = a
      end
  o.separator ''
  o.on("-o", "--output [OUTPUT_DIRECTORY]",
       "Local output 'datarec/' subdirectory name",
       "Default: [PRODUCTION_VERSION] (-> datarec/[PRODUCTION_VERSION])"
      ) do |a|
        options.locDir = a
      end
  o.separator ''
  o.on("-m", "--mode [MODE]",
       "specify what to do:",
       "  s - make config file for streaming from S3",
       "  d - download from S3, then make the local config file",
       "  c - just make the local config file, for downloaded local files",
       "Default: #{options.mode}"
      ) do |a|
        unless ['s','d','c'].include? a
          $stderr.puts "ERROR: unknown MODE '#{a}'"
          exit 1
        end
        options.mode = a
      end
  o.separator ''
  o.on("-l", "--limit [LIMIT]", Integer,
       "Specify a limit on the number of files",
       "  if an integer>0: take this many per Q2 bin",
       "  set to very high number to take all files (a lot of data!)",
       "Default: #{options.limit}"
      ) do |a|
        options.limit = a
      end
  o.separator ''
  o.on("-c", "--config [CONFIG_FILE]",
       "Name of config file; if not specified,",
       "it will be in datarec/[OUTPUT_DIRECTORY]/"
      ) do |a|
        options.configFile = a
      end
  o.separator ''
  o.on("-d", "--detector [DETECTOR_VERSION]",
       "Detector configuration version, one of:",
       *detectorConfigList.map{|e|"  #{e}"},
       "Default: #{options.detector}"
      ) do |a|
        unless detectorConfigList.include? a
          $stderr.puts "ERROR: unknown DETECTOR_VERSION '#{a}'"
          exit 1
        end
        options.detector = a
      end
  o.separator ''
  o.on("-r", "--[no-]radcor",
       "Decide whether to read radiative-corrected versions:",
       "  --radcor:     radiative corrections ON",
       "  --no-radcor:  radiative corrections OFF",
       "Default: --no-radcor"
      ) do |a|
        options.radCor = a
      end
  o.separator ''
  o.on_tail("-h", "--help",
            "Show this message"
           ) do
             puts o
             exit 2
           end
end.parse!( ARGV.length>0 ? ARGV : ['--help'] )
puts "OPTIONS: #{options}"

# get release and energy subdirectories, for the user-specified release version
prod = prodSettings[options.version]
prod[:releaseDir] = prod[:releaseSubDir].call
prod[:energyDir]  = prod[:releaseDir] + '/' + prod[:energySubDir].call
puts "Release Dir: #{prod[:releaseDir]}"
puts "Energy Dir:  #{prod[:energyDir]}"
# system "mc tree #{prod[:releaseDir]}"

# set target `locDir` directory
prod[:targetDir] = "datarec/#{options.locDir.empty? ? options.version : options.locDir}"


## helper functions
# get a list of files on S3 at `dir`
def mc_ls(dir)
  `mc ls #{dir}`
    .split(/\n/)
    .map{ |line| line.split.last }
end
# get the cross section from `CrossSectionTable`
getCrossSection = Proc.new do |searchPattern|
  rows = File.readlines(CrossSectionTable).grep(/^#{searchPattern} /)
  if rows.size == 1
    rows.first.split[1].to_f
  else
    $stderr.puts "ERROR: missing cross section for search pattern '#{searchPattern}'"    if rows.size == 0
    $stderr.puts "ERROR: duplicated cross section for search pattern '#{searchPattern}'" if rows.size > 1
    0.0
  end
end



# RELEASE VERSION DEPENDENT CODDE ##################
# Fill the following additional elements of `prod`:
# -> prod[:q2ranges]  list of pairs [minQ2,maxQ2] for each Q2 range
# -> prod[:dataDirs]  list of data directories for each Q2 range
# -> prod[:fileLists] list of files for each Q2 range
# -> prod[:radDir] the name of the radiative corrections directory (if applicable)
####################################################

if ['epic.22.11.3'].include? options.version
  # set source data directory
  prod[:radDir] = options.radCor ? 'radcor' : 'noradcor'
  dataDir = prod[:energyDir] + '/' + prod[:dataSubDir].call(prod[:radDir])
  puts "Data Dir: #{dataDir}"
  # set target directory
  prod[:targetDir] += "/"+prod[:radDir]
  puts "Target Dir: #{prod[:targetDir]}"
  # check if there are any files
  if mc_ls("#{dataDir} | head").empty?
    $stderr.puts "ERROR: no files in this Data Dir"
    puts "Available Data Directory Tree (be patient...)"
    system "mc tree #{prod[:releaseDir]}"
    exit 1
  end
  # get the Q2 ranges
  puts "Getting the full list of files... be patient..."
  fullList = mc_ls(dataDir).grep(/\.root$/)
  prod[:q2ranges] = fullList
    .map{ |file| file.gsub(/.*_q2_/,'').sub(/_run.*/,'') }
    .uniq
    .map{ |range| range.split('_').map &:to_i }
  puts "Q2 ranges: #{prod[:q2ranges]}"
  # set dataDirs (the same for each Q2 range)
  prod[:dataDirs] = prod[:q2ranges].map{ |q2range| dataDir }
  # get a list of files for each Q2 range
  puts "File names for each Q2 range:"
  prod[:fileLists] = prod[:q2ranges].map do |minQ2, maxQ2|
    fileList = fullList
      .grep(/q2_#{minQ2}_#{maxQ2}/)
      .first(options.limit)
    puts "--- #{minQ2} < Q2 < #{maxQ2}"
    fileList.each{ |file| puts "  #{file}" }
    fileList
  end


elsif ['epic.22.11.2'].include? options.version
  # print target directory
  puts "Target Dir: #{prod[:targetDir]}"
  # get list of Q2 subdirectories
  q2dirList = mc_ls prod[:energyDir]
  if q2dirList.empty?
    $stderr.puts "ERROR: energy not found"
    puts "Available energies"
    system "mc ls #{prod[:releaseDir]}"
    exit 1
  end
  puts "Q2 subdirectories: #{q2dirList}"
  # get the Q2 ranges
  prod[:q2ranges] = q2dirList.map do |dir|
    [ dir.split('=').last.sub(/\/$/,'').to_i, 0 ]
  end
  puts "Q2 ranges: #{prod[:q2ranges]}"
  # get a list of files for each Q2 range
  puts "File names for each Q2 range:"
  prod[:dataDirs] = []
  prod[:fileLists] = prod[:q2ranges].map do |minQ2, maxQ2|
    puts "--- #{minQ2} < Q2 < #{maxQ2}"
    dataDir = prod[:energyDir] + '/' + prod[:dataSubDir].call(minQ2)
    prod[:dataDirs] << dataDir
    puts "Data Dir: #{dataDir}"
    fileList = mc_ls(dataDir)
      .first(options.limit)
    puts "Files:"
    fileList.each{ |file| puts "  #{file}" }
    fileList
  end
  prod[:radDir] = '' # not used

end # END RELEASE VERSION DEPENDENT CODDE ##################


# append the energy to the target directory
prod[:targetDir] += '/'+options.energy
FileUtils.mkdir_p prod[:targetDir]

# download or stream the files, and build config file lists
localFileTableName = "#{prod[:targetDir]}/files.config.list"
localFileTable = File.open localFileTableName, 'w'
prod[:q2ranges].zip(prod[:dataDirs],prod[:fileLists]).each do |q2range,dataDir,fileList|

  # get the list of file names and set the target directory
  minQ2, maxQ2 = q2range
  targetDir = "#{prod[:targetDir]}/q2_#{minQ2}_#{maxQ2}"

  # download the files
  if options.mode=='d'
    FileUtils.mkdir_p targetDir, verbose: true
    fileList.each do |file|
      system "mc cp '#{dataDir}/#{file}' #{targetDir}/"
      puts ""
    end
  end

  # get the cross section
  crossSection = getCrossSection.call( prod[:crossSectionID].call minQ2, maxQ2, prod[:radDir] )

  # build `localFileTable`
  fileList.each do |fileBase|
    file = options.mode=='s' ?
      "#{dataDir.sub(/^S3/,'s3'+HostURL)}/#{fileBase}" : # if streaming, make S3 URL
      "#{targetDir}/#{fileBase}"                         # otherwise, use local target path
    localFileTable.puts "#{file} #{crossSection} #{minQ2} #{maxQ2}"
  end

end
localFileTable.close


# convert `localFileTable` into a full config file
puts '.'*50
puts "File Config Table: #{localFileTableName}"
system "cat #{localFileTableName}"
puts '\''*50
configFile = options.configFile.empty? ?
  localFileTableName.sub(/\.list/,'') :
  options.configFile
system "s3tools/generate-config-file.rb #{configFile} #{options.energy} #{localFileTableName}"
