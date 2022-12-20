#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2022 Christopher Dilks

require 'optparse'
require 'ostruct'
# require 'awesome_print'

# default CLI options
options = OpenStruct.new
options.version    = '22.11.3'
options.energy     = '18x275'
options.locDir     = "epic.#{options.version}"
options.mode       = 'd'
options.limit      = 5
options.configFile = ''
options.detector   = 'arches'
options.radCor     = false

# production specifications, latest first
prodSettings = {
  '22.11.3' => {
    :comment    => 'Pythia 6, with & without radiative corrections',
    :releaseDir => Proc.new { "S3/eictest/EPIC/RECO/#{options.version}/epic_#{options.detector}/SIDIS/pythia6" },
    :energyDir  => Proc.new { "ep_#{options.energy}" },
    :dataDir    => Proc.new { |radDir|
      if [options.energy,radDir]==['18x275','noradcor']  # correct for S3 disorganization
        "hepmc_ip6"
      else
        "hepmc_ip6/#{radDir}"
      end
    },
  },
  '22.11.2' => {
    :comment    => 'Pythia 8',
    :releaseDir => Proc.new { "S3/eictest/EPIC/RECO/#{options.version}/epic_#{options.detector}/DIS/NC" },
    :energyDir  => Proc.new { "#{options.energy}" },
    :dataDir    => Proc.new { |minQ2| "minQ2=#{minQ2}" },
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
       "Default: #{options.locDir} (-> datarec/#{options.locDir})"
      ) do |a|
        options.locDir = a
      end
  o.separator ''
  o.on("-m", "--mode [MODE]",
       "specify what to do:",
       "  s - make config file for streaming from S3",
       "  d - download from S3, then make the local config file",
       "  c - just make the local config file, for local files",
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


## helper functions
# get a list of files on S3 at `dir`
def mc_ls(dir)
  `mc ls #{dir}`
    .split(/\n/)
    .map{ |line| line.split.last }
end


# get energy subdirectory
prod = prodSettings[options.version]
releaseDir = prod[:releaseDir].call
energyDir = releaseDir + '/' + prod[:energyDir].call
puts "Release Dir: #{releaseDir}"
puts "Energy Dir:  #{energyDir}"
# system "mc tree #{releaseDir}"


# VERSION DEPENDENT SPAGHETTI ##############
# If downloading or streaming, get the Q2 ranges and file lists from S3
# -> prod[:q2ranges]  list of pairs [minQ2,maxQ2]
# -> prod[:fileLists] list of files for each Q2 range in prod[:q2ranges]
############################################

if mode=='s' or mode=='d'

  if ['22.11.3'].include? options.version
    dataDir = energyDir + '/' + prod[:dataDir].call( options.radCor ? 'radcor' : 'noradcor' )
    puts "Data Dir: #{dataDir}"
    # check if there are any files
    if mc_ls("#{dataDir} | head").empty?
      $stderr.puts "ERROR: no files in this Data Dir"
      puts "Available Data Directory Tree (be patient...)"
      system "mc tree #{releaseDir}"
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

  elsif ['22.11.2'].include? options.version
    # get list of Q2 subdirectories
    q2dirList = mc_ls energyDir
    if q2dirList.empty?
      $stderr.puts "ERROR: energy not found"
      puts "Available energies"
      system "mc ls #{releaseDir}"
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
    prod[:fileLists] = prod[:q2ranges].map do |minQ2, maxQ2|
      puts "--- #{minQ2} < Q2 < #{maxQ2}"
      dataDir = energyDir + '/' + prod[:dataDir].call(minQ2)
      puts "Data Dir: #{dataDir}"
      fileList = mc_ls(dataDir)
        .first(options.limit)
      puts "Files:"
      fileList.each{ |file| puts "  #{file}" }
      fileList
    end

  end # if version
end # if mode == 's' or 'd'


# build full S3 file names
prod[:fileLists].map!{ |file| "#{dataDir}/#{file}" }

# download the files
# prod[:q2range].zip(prod[:fileLists]).each do |q2range,fileList|
#   mkdir "#{option.locDir}/#{option.energy}/q2_#{minQ2}_#{maxQ2}"
#   fileList.each do |file ...
# end

