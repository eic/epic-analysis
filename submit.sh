#!/bin/bash

#SBATCH --job-name=athenaSIDIS
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=4000
#SBATCH --gres=disk:1000
#SBATCH --time=24:00:00
##SBATCH --mail-user=$USER@jlab.org

/work/clas12/users/$USER/eic/largex-eic/macro/dis-5x41/job.sh
