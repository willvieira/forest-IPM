#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 3-13:50:00
#SBATCH --mem-per-cpu=1000M
#SBATCH --ntasks=4
#SBATCH --job-name=stan_lm
#SBATCH --mail-user=willian.vieira@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-372

R -f runStan_i.R
