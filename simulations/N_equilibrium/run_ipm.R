# Script to run simulation on cluster

library(tidyverse)

spIds <- unique(
  gsub(
    '.*_|.csv', '', dir(file.path('..', 'data', 'parameters'))
  )
)

# create and export csv sim parameters
sim_pars <- expand_grid(
  species_id = spIds,
  N_init = 1,
  N_inter = 0,
  n_time = 2000,
  deltaTime = 1,
  param_method = 'random',
  replication = 1:50,
  plot_random = 0
) |>
mutate(
  seed = replication
)

sim_pars |>
  write_csv('simulation_pars.csv')


# write slurm bash file
bash_file <- paste0('#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 3-0:00:00
#SBATCH --mem-per-cpu=3000M
#SBATCH --ntasks=1
#SBATCH --job-name=ipm
#SBATCH --mail-user=willian.vieira@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-', nrow(sim_pars), '

module load StdEnv/2020 r/4.1.0

NCORES=\\$SLURM_CPUS_PER_TASK R -f ipm_i.R
')

# save bash file
writeLines(bash_file, 'sub.sh')

# run slurm
system('sbatch sub.sh')
