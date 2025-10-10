#!/usr/bin/env ruby
#
# run_postprocess_injectExtract.rb
#
# Drives the ROOT macro macros/postprocess_injectExtract.C for each
# (collision,energy) combination, either locally or via SLURM.
# Now supports splitting each energy into multiple CSVs/jobs.

require 'optparse'
require 'fileutils'
require 'time'

options = {
  project:           nil,
  asym_config_path:  nil,
  a_ut:              0.3,
  no_slurm:          false,
  n_splits:          1,
  suffix:            nil
}

OptionParser.new do |opts|
  opts.banner = "Usage: run_postprocess_injectExtract.rb -n <project> --asym-config-path <path> [options]"

  opts.on("-n NAME", "--name NAME", "Project name (required)") do |n|
    options[:project] = n.strip
  end

  opts.on("--asym-config-path PATH", "Directory of scheme_*.csv files (required)") do |p|
    options[:asym_config_path] = p.strip
  end

  opts.on("--A_UT VALUE", Float, "Asymmetry A_UT (default: 0.3)") do |v|
    options[:a_ut] = v
  end

  opts.on("--no-slurm", "Run all jobs locally instead of SLURM") do
    options[:no_slurm] = true
  end

  opts.on("--n-splits N", Integer, "Number of CSV splits per energy (default: 1)") do |n|
    options[:n_splits] = n
  end
    
  opts.on("--suffix STR", "Optional filename suffix (e.g., 'debug' -> base_debug_0000.csv)") do |s|
    options[:suffix] = s.strip
  end
    
  opts.on("-h", "--help", "Show this help message") do
    puts opts
    exit
  end
end.parse!

if options[:project].to_s.empty? || options[:asym_config_path].to_s.empty?
  abort "Error: --name and --asym-config-path are both required."
end

if options[:n_splits] < 1
  abort "Error: --n-splits must be >= 1."
end

# Determine zero-padding width (at least 4 digits)
if options[:n_splits] > 1
  max_idx = options[:n_splits] - 1
  digits  = Math.log10(max_idx).floor + 1
  padding = [4, digits].max
else
  padding = 4
end

COLLISIONS   = %w[ep]
ENERGIES     = %w[5x41 10x100 18x275]
TREE_NAME    = "dihadron_tree"
USE_DEPOL    = true
TARGET_POL   = 0.7

def make_and_run(project, collision, energy, cfg_dir, a_ut, tree_name, use_depol, target_pol, no_slurm, suffix, split_idx, padding)
  cfg_file = File.join(cfg_dir, "scheme_#{energy}.csv")
  unless File.exist?(cfg_file)
    warn "Skipping #{collision} #{energy}: missing #{cfg_file}"
    return
  end

  root_dir = File.join("out", "#{project}_#{energy}")
  FileUtils.mkdir_p(root_dir)
  inj_dir  = File.join(root_dir, "injectionResults")
  FileUtils.mkdir_p(inj_dir)

  base     = File.basename(cfg_file, ".csv")
  idx_str  = split_idx.to_s.rjust(padding, '0')
  suffix_str = (suffix && !suffix.empty?) ? "_#{suffix}" : ""
  out_csv  = File.join(inj_dir, "#{base}#{suffix_str}_#{idx_str}.csv")

  args     = "\"#{root_dir}\",\"#{tree_name}\",\"#{cfg_file}\",\"#{out_csv}\",#{a_ut},#{use_depol},#{target_pol},\"#{collision}\""
  cmd      = "root -l -b -q 'macro/dihadron_IFF/postprocess_injectExtract.C(#{args})'"

  if no_slurm
    puts "[LOCAL] #{collision.upcase} @ #{energy} (split #{idx_str}) → #{out_csv}"
    system(cmd) or warn "  Command failed: #{cmd}"
  else
    yield collision, energy, cmd, out_csv, idx_str
  end
end

if options[:no_slurm]
  puts "=== Running all jobs locally (#{options[:n_splits]} splits per energy) ==="
  COLLISIONS.each do |c|
    ENERGIES.each do |e|
      (0...options[:n_splits]).each do |i|
        make_and_run(
          options[:project],
          c, e,
          options[:asym_config_path],
          options[:a_ut],
          TREE_NAME,
          USE_DEPOL,
          TARGET_POL,
          true,
          i,
          padding
        ) { }  # no block needed for local runs
      end
    end
  end
else
  # SLURM settings
  ACCOUNT       = "eic"
  PARTITION     = "production"
  MEM_PER_CPU   = 4000    # MB
  CPUS_PER_TASK = 2
  TIME_LIMIT    = "24:00:00"

  timestamp = Time.now.strftime("%Y-%m-%d___%H-%M-%S")
  slurm_dir = File.join("hpc","slurm",timestamp)
  log_dir   = File.join(slurm_dir,"log")
  FileUtils.mkdir_p(log_dir)

  run_script = "run_postprocess_jobs.sh"
  File.open(run_script,"w") { |f| f.puts("#!/bin/bash") }

  COLLISIONS.each do |c|
    ENERGIES.each do |e|
      (0...options[:n_splits]).each do |i|
        make_and_run(
          options[:project],
          c, e,
          options[:asym_config_path],
          options[:a_ut],
          TREE_NAME,
          USE_DEPOL,
          TARGET_POL,
          false,
          options[:suffix],
          i,
          padding
        ) do |collision, energy, cmd, out_csv, idx_str|
          job    = "postproc_#{options[:project]}_#{collision}_#{energy}_#{idx_str}"
          slurm  = File.join(slurm_dir, "#{job}.slurm")
          stdout = File.join(log_dir, "#{job}.out")
          stderr = File.join(log_dir, "#{job}.err")

          File.open(slurm,"w") do |f|
            f.puts("#!/bin/bash")
            f.puts("#SBATCH --account=#{ACCOUNT}")
            f.puts("#SBATCH --partition=#{PARTITION}")
            f.puts("#SBATCH --job-name=#{job}")
            f.puts("#SBATCH --cpus-per-task=#{CPUS_PER_TASK}")
            f.puts("#SBATCH --mem-per-cpu=#{MEM_PER_CPU}")
            f.puts("#SBATCH --time=#{TIME_LIMIT}")
            f.puts("#SBATCH --output=#{stdout}")
            f.puts("#SBATCH --error=#{stderr}")
            f.puts
            f.puts("cd #{Dir.pwd}")
            f.puts(cmd)
          end

          FileUtils.chmod("+x", slurm)
          File.open(run_script,"a") { |f| f.puts("sbatch #{slurm}") }
          puts "Prepared SLURM job #{job}; output CSV: #{out_csv}"
        end
      end
    end
  end

  puts "\nAll SLURM scripts written to #{slurm_dir}"
  puts "Submit them via:\n    bash #{run_script}"
end
