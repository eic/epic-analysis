#!/usr/bin/env ruby
require 'fileutils'
require 'yaml'
require 'optparse'

PWD=`pwd`.strip

# Parse command line options
options = {}
OptionParser.new do |opts|
  opts.on('--runcard FILE', 'YAML runcard file for configuration') do |f|
    options[:runcard] = f
  end
end.parse!

#################################################################################
#                          USER CONFIGURATION SECTION                            #
# Edit the parameters below to define your simulation campaign.                  #
# If a runcard is provided via --runcard, these defaults will be overridden.     #
#################################################################################
 
# Load configuration from runcard if provided
if options[:runcard]
  unless File.exist?(options[:runcard])
    puts "Error: Runcard file #{options[:runcard]} does not exist"
    exit 1
  end
  config = YAML.load_file(options[:runcard])
else
  config = {}
end

# ---------------------------- Project Metadata -------------------------------- #
PROJECT_NAME = config['project_name'] || "test"   # Prefix for output directory names
CAMPAIGNS    = config['campaigns'] || ["epic.25.08.0"]      # Simulation campaign tags (e.g., "epic.25.08.0")
DETECTORS    = config['detectors'] || ["epic_craterlake"]   # Detector configurations (must align with CAMPAIGNS)

# ------------------------------- Beam Settings -------------------------------- #
# Example: ["5x41", "10x100", "18x275"] for multiple energies per campaign
ENERGIES = config['energies'] || ["10x166"]              # Nested arrays: one list of energies per campaign

# ---------------------------- Simulation Mode --------------------------------- #
if config['is_beagle'].nil?
  IS_BEAGLE = false                  # default to false if not specified
else
  IS_BEAGLE = config['is_beagle']    # use specified value from runcard
end

# -------------------------- File and Job Parameters --------------------------- #
NFILES               = config['nfiles'] || 500_000       # Max number of files per Q² bin (large = all)
NROOT_FILES_PER_JOB  = config['nroot_files_per_job'] || 40            # ROOT files merged per job

# ------------------------------- Analysis Setup ------------------------------- #
PATH_TO_ANALYSIS_MACRO = config['path_to_analysis_macro'] || "macro/analysis_dihadron.C"   # ROOT analysis macro to run
PATH_TO_EIC_SHELL      = config['path_to_eic_shell'] || "#{ENV['EIC_SHELL_PREFIX']}/../"  # Path to eic-shell directory

#################################################################################
# End of user configuration section                                              #
#################################################################################

# Error Handling
# Check if project name is empty
if PROJECT_NAME.empty?
  puts "Error: PROJECT_NAME is empty"
  exit 1
end

# Check if NFILES is less than or equal to 0
if NFILES <= 0
  puts "Error: NFILES must be greater than 0"
  exit 1
end

# Check if NROOT_FILES_PER_JOB is less than or equal to 0
if NROOT_FILES_PER_JOB <= 0
  puts "Error: NROOT_FILES_PER_JOB must be greater than 0"
  exit 1
end

# Check if PATH_TO_ANALYSIS_MACRO file exists
unless File.exist?(PATH_TO_ANALYSIS_MACRO)
  puts "Error: File specified in PATH_TO_ANALYSIS_MACRO does not exist"
  exit 1
end

#################################################################################
#################################################################################
#################################################################################
#################################################################################
USERNAME = ENV['USER'] || ENV['USERNAME']

# Check for --overwrite argument
overwrite = ARGV.include?('--overwrite')

outdir=""
final_slurm_scripts=[]
# Loop over each campaign
CAMPAIGNS.each_with_index do |campaign, index|
  puts "Processing Campaign: #{campaign}"
  puts "--------------------------------------------------------------------"

  detector = DETECTORS[index]
  # Loop over energies for the current campaign
  ENERGIES.each do |energy|

    # create output directory
    outdir="#{PROJECT_NAME}___#{campaign}_#{energy}"
    dir_path = "#{PWD}/out/#{outdir}"

    # Check if the directory exists
    
    if Dir.exist?(dir_path) && !overwrite
      puts "Directory #{dir_path} already exists. Delete it? [y/N]"
      user_input = gets.chomp.downcase
      if user_input == 'y'
        FileUtils.rm_rf(dir_path)
        puts "Directory #{dir_path} has been deleted."
      else
        puts "Operation aborted by the user."
        exit
      end
    elsif Dir.exist?(dir_path) && overwrite
      FileUtils.rm_rf(dir_path)
      puts "Overwriting: Directory #{dir_path} has been deleted."
    end

    # Create the directory
    FileUtils.mkdir_p(dir_path)
    puts "Creating project --> #{dir_path}"

    # Copy the analysis macro to the output directory for record keeping
    FileUtils.cp(PATH_TO_ANALYSIS_MACRO, "#{dir_path}/#{File.basename(PATH_TO_ANALYSIS_MACRO)}")
    puts "Copied analysis macro to #{dir_path}/#{File.basename(PATH_TO_ANALYSIS_MACRO)}"


    # Grab the files from s3
    # If BeAGLE is used, automatically set the --target flag
    target_flag = IS_BEAGLE ? "--target He3" : ""
    puts "#{PWD}/s3tools/s3tool.rb -e #{energy} -o #{outdir} -l #{NFILES} -v #{campaign} #{target_flag}"
    puts `#{PWD}/s3tools/s3tool.rb -e #{energy} -o #{outdir} -l #{NFILES} -v #{campaign} #{target_flag}`

    if IS_BEAGLE
        puts "BeAGLE mode detected → using --target He3"
    else
        puts "Standard Pythia mode (no --target)"
    end       
      
    # create shell script to count nevents
     File.open("#{PWD}/out/#{outdir}/count-nevents.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Counting events for campaign #{campaign} with beam energy #{energy}...(may take a while)...\""
      file.puts "echo \"Results are stored for later access in hpc/nevents_databases for faster computation...\""
      file.puts "python3 #{PWD}/hpc/src/count_events.py #{campaign} #{detector} #{energy}"
    end

    puts "Shell script created at #{PWD}/out/#{outdir}/count-nevents.sh"

    # create shell script for generating the config files and slurm files
    File.open("#{PWD}/out/#{outdir}/make-configs.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "source #{PWD}/environ.sh"
      file.puts "#{PWD}/hpc/prepare-multi-roots.rb datarec/#{outdir}/#{energy}/files.config #{PWD}/out/#{outdir} #{NROOT_FILES_PER_JOB}"
      file.puts "echo \"y\" | #{PWD}/hpc/run-local-slurm.rb #{PATH_TO_ANALYSIS_MACRO} #{PWD}/out/#{outdir} #{outdir}"
    end

    puts "Shell script created at #{PWD}/out/#{outdir}/make-configs.sh"

    # create shell script for running the 'run.slurm' objective created by hpc/run-local-slurm.rb
    # This is a separate script that runs outside the eic-shell environment
    File.open("#{PWD}/out/#{outdir}/run-parallel.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Running run.slurm...\""
      file.puts "sbatch --wait #{PWD}/out/#{outdir}/scripts/run.slurm"
    end

    puts "Shell script created at #{PWD}/out/#{outdir}/run-parallel.sh"

    # create shell script for running the 'merge.rb' function to merge the TTrees
    File.open("#{PWD}/out/#{outdir}/merge.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Merging TTrees...\""
      file.puts "source #{PWD}/environ.sh"
      file.puts "#{PWD}/hpc/merge.rb #{PWD}/out/#{outdir}/"
    end

    puts "Shell script created at #{PWD}/out/#{outdir}/merge.sh"

    # create slurm script that will run the above shell scripts
    File.open("#{PWD}/out/#{outdir}/run-pipeline.slurm", 'w') do |file|
      file.puts """#!/bin/bash
#SBATCH --job-name=#{outdir}
#SBATCH --account=eic
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --time=24:00:00
#SBATCH --output=#{PWD}/out/#{outdir}/pipeline.out
#SBATCH --error=#{PWD}/out/#{outdir}/pipeline.err

bash #{PWD}/out/#{outdir}/count-nevents.sh
#{PATH_TO_EIC_SHELL}/eic-shell -- "
  cd #{PWD}
  #{PWD}/out/#{outdir}/make-configs.sh
  bash #{PWD}/out/#{outdir}/run-parallel.sh
  #{PWD}/out/#{outdir}/merge.sh
"
      """
    end


    puts "Slurm script created at #{PWD}/out/#{outdir}/run-pipeline.slurm"
    FileUtils.chmod('+x', "#{PWD}/out/#{outdir}/make-configs.sh")
    FileUtils.chmod('+x', "#{PWD}/out/#{outdir}/run-parallel.sh")
    FileUtils.chmod('+x', "#{PWD}/out/#{outdir}/merge.sh")
    final_slurm_scripts << "#{PWD}/out/#{outdir}/run-pipeline.slurm"
    puts "--------------------------------------------------------------------"
    # create 
  end # end loop energies
end # end loop campaigns

File.open("hpc/project_scripts/run-#{PROJECT_NAME}.sh", 'w') do |file|
  file.puts "#!/bin/bash"
  file.puts ""
  final_slurm_scripts.each do |script|
      file.puts "sbatch #{script}"
  end
end  

puts "Completed setup ... "
puts "To execute the analysis, run the following **OUTSIDE** the eic-shell environment"
puts "bash #{PWD}/hpc/project_scripts/run-#{PROJECT_NAME}.sh "
