#!/usr/bin/env ruby
require 'fileutils'
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Manually edit these according to your desired simulation

# prefix for output
PROJECT_NAME=""

CAMPAIGNS=["epic.23.10.0",
           "epic.23.11.0"]

DETECTORS = ["epic_craterlake",
             "epic_craterlake"]

ENERGIES=[ ["5x41","10x100","18x275"],
           ["5x41","10x100","18x275"] ]

# Number of Files per Q2 binning
# Set this to very large number for all campaign files
NFILES = 0

NROOT_FILES_PER_JOB = 0

# Points to analysis macro
PATH_TO_ANALYSIS_MACRO = "macro/MY_MACRO.C"

# Path to the directory containing eic-shell
PATH_TO_EIC_SHELL = "#{ENV['EIC_SHELL_PREFIX']}/../"

#################################################################################
#################################################################################
#################################################################################
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
  ENERGIES[index].each do |energy|

    # create output directory
    outdir="#{PROJECT_NAME}___#{campaign}_#{energy}"
    dir_path = "out/#{outdir}"

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
    puts "Creating project --> out/#{outdir}"
    

    # Grab the files from s3
    puts `./s3tools/s3tool.rb -e #{energy} -o #{outdir} -l #{NFILES} -v #{campaign}`
      
      
    # create shell script to count nevents
     File.open("out/#{outdir}/count-nevents.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Counting events for campaign #{campaign} with beam energy #{energy}...(may take a while)...\""
      file.puts "echo \"Results are stored for later access in hpc/nevents_databases for faster computation...\""
      file.puts "python3 ./hpc/src/count_events.py #{campaign} #{detector} #{energy}"
    end

    puts "Shell script created at out/#{outdir}/count-nevents.sh"
      
      
    # create shell script for generating the config files and slurm files    
    File.open("out/#{outdir}/make-configs.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "source environ.sh"
      file.puts "hpc/prepare-multi-roots.rb datarec/#{outdir}/#{energy}/files.config out/#{outdir} #{NROOT_FILES_PER_JOB}"
      file.puts "echo \"y\" | hpc/run-local-slurm.rb #{PATH_TO_ANALYSIS_MACRO} out/#{outdir} #{outdir}"
    end

    puts "Shell script created at out/#{outdir}/make-configs.sh"
      
    # create shell script for running the 'run.slurm' objective created by hpc/run-local-slurm.rb
    # This is a separate script that runs outside the eic-shell environment
    File.open("out/#{outdir}/run-parallel.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Running run.slurm...\""
      file.puts "sbatch --wait out/#{outdir}/scripts/run.slurm"
    end  
    
    puts "Shell script created at out/#{outdir}/run-parallel.sh"
      
    # create shell script for running the 'merge.rb' function to merge the TTrees
    File.open("out/#{outdir}/merge.sh", 'w') do |file|
      file.puts "#!/bin/bash"
      file.puts ""
      file.puts "echo \"Merging TTrees...\""
      file.puts "source environ.sh"
      file.puts "hpc/merge.rb out/#{outdir}/"
    end  
    
    puts "Shell script created at out/#{outdir}/merge.sh"
      
    # create slurm script that will run the above shell scripts
    File.open("out/#{outdir}/run-pipeline.slurm", 'w') do |file|
      file.puts """#!/bin/bash
#SBATCH --job-name=#{outdir}
#SBATCH --account=eic
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --time=24:00:00
#SBATCH --output=out/#{outdir}/pipeline.out
#SBATCH --error=out/#{outdir}/pipeline.err

bash out/#{outdir}/count-nevents.sh
#{PATH_TO_EIC_SHELL}/eic-shell -- out/#{outdir}/make-configs.sh
bash out/#{outdir}/run-parallel.sh
#{PATH_TO_EIC_SHELL}/eic-shell -- out/#{outdir}/merge.sh
      """
    end
    
    
    puts "Slurm script created at out/#{outdir}/run-pipeline.slurm"
    FileUtils.chmod('+x', "out/#{outdir}/make-configs.sh")
    FileUtils.chmod('+x', "out/#{outdir}/run-parallel.sh")
    FileUtils.chmod('+x', "out/#{outdir}/merge.sh")
    final_slurm_scripts << "out/#{outdir}/run-pipeline.slurm"
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
puts "bash hpc/project_scripts/run-#{PROJECT_NAME}.sh "
