#!/usr/bin/env ruby

# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright (C) 2023 Christopher Dilks

require 'optparse'
require 'ostruct'
require 'fileutils'

# default versions
VersionLatest   = 'epic.23.11.0'
VersionPrevious = 'epic.23.10.0'

# default CLI options
options = OpenStruct.new
options.version    = VersionLatest
options.energy     = '18x275'
options.locDir     = ''
options.mode       = 's'
options.limit      = 2
options.configFile = ''
options.detector   = 'craterlake'
options.radCor     = false
options.minQ2      = -1
options.maxQ2      = -1
options.numHepmc   = 0

# global settings
CrossSectionTable = 'datarec/xsec/xsec.dat'
HostURL           = 'https://eics3.sdcc.bnl.gov:9000'

# helpers
def versionNum(v) # options.version -> version number
  v
    .split('_').first
    .split('.')[1..-1]
    .join('.')
end
def ecceQ2range(minQ2,maxQ2) # return file path suffix, for ECCE Q2 ranges
  { [1,0]=>'', [1,100]=>'-q2-low', [100,0]=>'-q2-high' }[[minQ2,maxQ2]]
end

# production specifications, latest first
# PRODUCTION_VERSION => {
#   :comment         => Description about this version
#   :crossSectionID  => Proc(minQ2,maxQ2,radDir) -> label of row of `CrossSectionTable`
#   :releaseSubDir   => Proc() -> Directory for this PRODUCTION_VERSION
#   :energySubDir    => Proc() -> Subdirectory associated to user-specified beam energy
#   :dataSubDir      => Proc(*version dependent*) -> Subdirectory of `:energySubDir`
#   :fileExtension   => File extension (optional, defaults to 'root')
# }
prodSettings = {
  'epic.23.11.0' => {
    :comment         => 'Pythia 8: high-stats November 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.10.0' => {
    :comment         => 'Pythia 8: high-stats November 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.09.1' => {
    :comment         => 'Pythia 8: high-stats ??? 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.07.1' => {
    :comment         => 'Pythia 8: high-stats July 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.06.1' => {
    :comment         => 'Pythia 8: high-stats June 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.05.2' => {
    :comment         => 'Pythia 8: high-stats May 2023 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.05.1' => {
    :comment         => 'Pythia 8',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.23.03.0_pythia8' => {
    :comment         => 'Pythia 8, small sample, 10x100, minQ2=1000 only',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  # 'epic.23.03.0_pythia6' => { # FIXME: need cross section for Q2<1 bin
  #   :comment         => 'Pythia 6, small sample, 5x41, noradcor only, Q2<1 only',
  #   :crossSectionID  => Proc.new { |minQ2,maxQ2,radDir| "pythia6:ep_#{radDir}.#{options.energy}_q2_#{minQ2}_#{maxQ2}" },
  #   :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/SIDIS/pythia6" },
  #   :energySubDir    => Proc.new { "ep_#{options.energy}" },
  #   :dataSubDir      => Proc.new { |radDir| "hepmc_ip6/#{radDir}" },
  # },
  'epic.23.01.0' => {
    :comment         => 'Pythia 8',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'epic.22.11.3' => {
    :comment         => 'Pythia 6: high-stats November 2022 production, with & without radiative corrections',
    :crossSectionID  => Proc.new { |minQ2,maxQ2,radDir| "pythia6:ep_#{radDir}.#{options.energy}_q2_#{minQ2}_#{maxQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/SIDIS/pythia6" },
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
    :comment         => 'Pythia 8: high-stats November 2022 production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/RECO/#{versionNum(options.version)}/epic_#{options.detector}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'ecce.22.1' => {
    :comment         => 'Last ECCE Production, August 2022',
    :crossSectionID  => Proc.new { |minQ2,maxQ2| "pythia6:ep-#{options.energy}#{ecceQ2range(minQ2,maxQ2)}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/Campaigns/#{versionNum(options.version)}/SIDIS/pythia6" },
    :energySubDir    => Proc.new { "ep-#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2,maxQ2| ecceQ2range(minQ2,maxQ2) }, # (combined with :energySubDir)
  },
  'athena.deathvalley-v1.0' => {
    :comment         => 'ATHENA Proposal production',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/ATHENA/RECO/#{versionNum(options.version)}/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
  },
  'hepmc.pythia6' => {
    :comment         => 'HEPMC files from Pythia 6 for ePIC, with & without radiative corrections',
    :crossSectionID  => Proc.new { |minQ2,maxQ2,radDir| "pythia6:ep_#{radDir}.#{options.energy}_q2_#{minQ2}_#{maxQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/EPIC/EVGEN/SIDIS/pythia6" },
    :energySubDir    => Proc.new { "ep_#{options.energy}" },
    :dataSubDir      => Proc.new { |radDir| "hepmc_ip6/#{radDir}" },
    :fileExtension   => 'hepmc',
  },
  'hepmc.pythia8' => {
    :comment         => 'HEPMC files from Pythia 8, used for ATHENA proposal',
    :crossSectionID  => Proc.new { |minQ2| "pythia8:#{options.energy}/minQ2=#{minQ2}" },
    :releaseSubDir   => Proc.new { "S3/eictest/ATHENA/EVGEN/DIS/NC" },
    :energySubDir    => Proc.new { "#{options.energy}" },
    :dataSubDir      => Proc.new { |minQ2| "minQ2=#{minQ2}" },
    :fileExtension   => 'hepmc.gz',
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
  "craterlake",
]

# parse options
OptionParser.new do |o|
  o.banner = "USAGE: #{$0} [OPTIONS]..."
  o.separator ''
  o.separator 'OPTIONS:'
  o.on("-v", "--version [PRODUCTION_VERSION]",
       "Production campaign version",
       "Default: #{options.version}",
       "Choose one of the following:"
      ) do |a|
        ver = a
        ver = VersionLatest   if a == 'epic.latest'
        ver = VersionPrevious if a == 'epic.previous'
        unless prodSettings.keys.include? ver
          $stderr.puts "ERROR: unknown DETECTOR_VERSION '#{ver}'"
          exit 1
        end
        options.version = ver
      end
  o.separator ''
  o.separator 'epic.latest'.rjust(30)   + '  ::  ' + "Latest production: #{VersionLatest}"
  o.separator 'epic.previous'.rjust(30) + '  ::  ' + "Previous production: #{VersionPrevious}"
  prodSettings.map{ |k,v| o.separator k.rjust(30) + '  ::  ' + v[:comment] }
  o.separator ''
  o.separator ' '*20+'NOTE: if the version name starts with \'hepmc\', HEPMC files will be'
  o.separator ' '*20+'      downloaded to datagen/ and passed through Delphes fast simulation'
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
  o.separator 'Options useful for CI:'
  o.on("--minQ2 [MIN_Q2]", Float, "limit to Q2 bins with minQ2=[MIN_Q2]") { |a| options.minQ2 = a }
  o.on("--maxQ2 [MAX_Q2]", Float, "limit to Q2 bins with maxQ2=[MAX_Q2]") { |a| options.maxQ2 = a }
  o.on("--num-hepmc-events [NUM]", Integer, "limit the number of HEPMC events per file", "(for production version 'hepmc.pythia6' only)") { |a| options.numHepmc = a }
  o.separator ''
  o.separator ''
  o.on_tail("-h", "--help",
            "Show this message"
           ) do
             puts o
             exit 2
           end
end.parse!( ARGV.length>0 ? ARGV : ['--help'] )
puts "OPTIONS: #{options}"

# warn about large HEPMC files
if options.version=='hepmc.pythia6' and options.numHepmc<=0
  puts """
    WARNING: unfortunately, these HEPMC files are very large!!!
      Recommendation: limit the number of events per file with the
                      `--num-hepmc-events` option
  """
end
if options.numHepmc>0 and options.version!='hepmc.pythia6'
  $stderr.puts "WARNING: --num-hepmc-events option does not apply to production version '#{options.version}'"
end

# get release and energy subdirectories, for the user-specified release version
prod = prodSettings[options.version]
prod[:releaseDir] = prod[:releaseSubDir].call
prod[:energyDir]  = prod[:releaseDir] + '/' + prod[:energySubDir].call
puts "Release Dir: #{prod[:releaseDir]}"
puts "Energy Dir:  #{prod[:energyDir]}"
# system "mc tree #{prod[:releaseDir]}"

# set target `locDir` directory
prod[:targetDir] = "datarec/#{options.locDir.empty? ? options.version : options.locDir}"

# set the file extension, if specified
ext = prod[:fileExtension].nil? ? 'root' : prod[:fileExtension]

# settings for handling full simulation output `RECO` vs. fast simulation input `EVGEN` from event generation
readingEvGen = options.version.match? /^hepmc/
delphesCmd   = 's3tools/src/loop_run_delphes.sh'

## helper functions
# get a list of files on S3 at `dir`; use `preFilter` to filter out things that are not in the file name (e.g., file size)
def mc_ls(dir, preFilter='')
  ls = `mc ls #{dir}`.split(/\n/)
  ls = ls.grep_v(preFilter) unless preFilter==''
  ls.map{ |line| line.split.last }
end
# download a file from S3 (and do not clobber)
def mc_cp(srcfile,tgtdir)
  tgtfile = "#{tgtdir}/#{File.basename srcfile}"
  if File.exist? tgtfile
    puts "File already exists: #{tgtfile}"
  else
    system "mc cp '#{srcfile}' #{tgtdir}/"
    puts ""
  end
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
# remove Q2 bins as specified by options.minQ2 and options.maxQ2
cullQ2bins = Proc.new do |q2ranges|
  q2ranges.select! do |minQ2,maxQ2|
    (minQ2==options.minQ2 or options.minQ2<0) and (maxQ2==options.maxQ2 or options.maxQ2<0)
  end
end



# RELEASE VERSION DEPENDENT CODE ###################
# Organize a list of Q2 ranges and associated S3 files for each
# - The file tree layout and naming conventions on S3 varies as a function of productions
# - The following is a chain of `if` statements, grouping together productions with similarly structured file trees
# - Each if-block must fill the following additional elements of `prod`:
#   - prod[:q2ranges]  list of pairs [minQ2,maxQ2] for each Q2 range
#   - prod[:dataDirs]  list of data directories for each Q2 range
#   - prod[:fileLists] list of files for each Q2 range
#   - prod[:radDir]    the name of the radiative corrections directory (if applicable)
####################################################

# pattern: "ep_#{energy}/hepmc_ip6/" with Q2 range given in file name as "q2_#{minQ2}_#{maxQ2}"
if [
    'epic.23.03.0_pythia6',
    'epic.22.11.3',
    'hepmc.pythia6',
].include? options.version
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
  fullList = mc_ls(dataDir).grep(/\.#{ext}$/)
  prod[:q2ranges] = fullList
    .map{ |file| file.gsub(/.*_q2_/,'').sub(/_run.*/,'') }
    .uniq
    .map{ |range| range.split('_').map &:to_i }
  cullQ2bins.call prod[:q2ranges]
  puts "Q2 ranges: #{prod[:q2ranges]}"
  # set dataDirs (the same for each Q2 range)
  prod[:dataDirs] = prod[:q2ranges].map{ |q2range| dataDir }
  # get a list of files for each Q2 range
  puts "File names for each Q2 range:"
  prod[:fileLists] = prod[:q2ranges].map do |minQ2, maxQ2|
    fileList = fullList
      .grep(/q2_#{minQ2}_#{maxQ2}/)
      .first(options.limit)
    puts "--- #{minQ2} < Q2 < #{maxQ2>0?maxQ2:'inf'}"
    fileList.each{ |file| puts "  #{file}" }
    fileList
  end

# pattern: "#{energy}/minQ2=#{minQ2}/"
elsif [
  'epic.23.11.0',
  'epic.23.10.0',
  'epic.23.09.1',
  'epic.23.07.1',
  'epic.23.06.1',
  'epic.23.05.2',
  'epic.23.05.1',
  'epic.23.03.0_pythia8',
  'epic.23.01.0',
  'epic.22.11.2',
  'athena.deathvalley-v1.0',
  'hepmc.pythia8',
].include? options.version
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
  cullQ2bins.call prod[:q2ranges]
  puts "Q2 ranges: #{prod[:q2ranges]}"
  # get a list of files for each Q2 range
  puts "File names for each Q2 range:"
  prod[:dataDirs] = []
  prod[:fileLists] = prod[:q2ranges].map do |minQ2, maxQ2|
    puts "--- #{minQ2} < Q2 < #{maxQ2>0?maxQ2:'inf'}"
    dataDir = prod[:energyDir] + '/' + prod[:dataSubDir].call(minQ2)
    prod[:dataDirs] << dataDir
    puts "Data Dir: #{dataDir}"
    fileList = readingEvGen ?
      mc_ls(dataDir,/GiB/) .grep(/\.#{ext}$/) .grep(/vtxfix/) .first(options.limit) :
      mc_ls(dataDir)       .grep(/\.#{ext}$/) .first(options.limit)
    puts "Files:"
    fileList.each{ |file| puts "  #{file}" }
    fileList
  end
  prod[:radDir] = '' # not used

elsif [
  'ecce.22.1'
].include? options.version
  prod[:radDir] = '' # not used
  prod[:q2ranges] = mc_ls(prod[:releaseDir]).grep(/#{options.energy}/).grep_v(/Lambda/).map do |dir|
    if dir.match? /-q2-high/
      [100,0]
    elsif dir.match? /-q2-low/
      [1,100]
    else
      [1,0]
    end
  end
  cullQ2bins.call prod[:q2ranges]
  prod[:dataDirs] = prod[:q2ranges].map do |minQ2,maxQ2|
    prod[:energyDir] + prod[:dataSubDir].call(minQ2,maxQ2)
  end
  prod[:fileLists] = prod[:dataDirs].map do |dataDir|
    mc_ls(dataDir)
      .grep(/g4event_eval.root$/)
      .first(options.limit)
  end

else
  $stderr.puts "\nERROR: production version '#{options.version}' is missing in RELEASE VERSION DEPENDENT part of #{$0}; add it!"
  exit 1
end # END RELEASE VERSION DEPENDENT CODE ##################


# append the energy to the target directory
prod[:targetDir] += '/'+options.energy
FileUtils.mkdir_p prod[:targetDir]

# download or stream the files, and build config file lists
localFileTableName = options.configFile.empty? ?
  "#{prod[:targetDir]}/files.config.list" :
  options.configFile + '.list'
localFileTable = File.open localFileTableName, 'w'
prod[:q2ranges].zip(prod[:dataDirs],prod[:fileLists]).each do |q2range,dataDir,fileList|

  # get the list of file names and set the target directory
  minQ2, maxQ2 = q2range
  targetDir = "#{prod[:targetDir]}/q2_#{minQ2}_#{maxQ2}"

  # download the RECO files
  if options.mode=='d' and !readingEvGen
    puts "DOWNLOADING RECO FILES FROM S3..."
    FileUtils.mkdir_p targetDir, verbose: true
    fileList.each do |file|
      mc_cp "#{dataDir}/#{file}", targetDir
    end
  # or download the EVGEN files
  elsif readingEvGen
    puts "DOWNLOADING EVGEN FILES FROM S3..."
    genDir = targetDir.sub /^datarec/, 'datagen'
    delphesCmd += " #{genDir}"
    FileUtils.mkdir_p genDir, verbose: true
    fileList.each do |file|
      if options.version=="hepmc.pythia6" and options.numHepmc>0
        # truncate HEPMC file after `options.numHepmc` events (FIXME: could be more efficient)
        puts "Downloading #{file}, truncating after #{options.numHepmc} events..."
        lineNum = `mc cat #{dataDir}/#{file} | grep -En -m#{options.numHepmc+1} '^E' | tail -n1 | sed 's/:.*//g'`.to_i
        puts "  Event #{options.numHepmc} ends on line number #{lineNum-1}; now downloading..."
        system "mc head -n #{lineNum-1} #{dataDir}/#{file} > #{genDir}/#{file}"
        puts "  ...done"
      else
        # otherwise get the full file
        mc_cp "#{dataDir}/#{file}", genDir
      end
    end
  end

  # get the cross section
  crossSection = getCrossSection.call( prod[:crossSectionID].call minQ2, maxQ2, prod[:radDir] )

  # build `localFileTable`
  fileList.each do |fileBase|
    file = "#{targetDir}/#{fileBase}"
    if readingEvGen # if reading EVGEN, be sure to use `.root` extension
      file.sub! /\.#{ext}$/, '.root'
    elsif options.mode=='s' # if streaming, make URL
      file = "#{dataDir.sub(/^S3/,'s3'+HostURL)}/#{fileBase}"
    end
    localFileTable.puts "#{file} #{crossSection} #{minQ2} #{maxQ2}"
  end

end
localFileTable.close

# run Delphes, if reading EVGEN
if readingEvGen
  puts """
RUNNING DELPHES with the following command:
```
#{delphesCmd}
```"""
  system delphesCmd
end

# convert `localFileTable` into a full config file
puts '.'*50
puts "File Config Table: #{localFileTableName}"
system "cat #{localFileTableName}"
puts '\''*50
configFile = localFileTableName.sub(/\.list/,'')
system "s3tools/src/generate-config-file.rb #{configFile} #{options.energy} #{localFileTableName}"
