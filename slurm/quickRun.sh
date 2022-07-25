#!/bin/bash
#SBATCH --account=eic
#SBATCH --partition=production
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name=quick-eic-analysis
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --chdir=/work/clas12/users/gmat/eic/sidis-eic
#SBATCH --output=/work/clas12/users/gmat/eic/sidis-eic/quickRun.out
#SBATCH --error=/work/clas12/users/gmat/eic/sidis-eic/quickRun.err
root macro/analysisEE_asymmetry.C
