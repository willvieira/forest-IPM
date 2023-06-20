# Script to run simulation on cluster

library(tidyverse)

spIds <- unique(
  gsub(
    '.*_|.csv', '', dir(file.path('..', 'data', 'parameters'))
  )
)

# create and export csv sim parameters
sim_pars <- expand_grid(
    sp1 = spIds,
    sp2 = spIds,
    Ninit_sp1 = 0.7,
    Ninit_sp2 = 0.7,
    n_time = 2000,
    deltaTime = 1,
    plotSize = 180,
    param_method = 'random',
    replication = 1:50,
    plot_random = 0
  ) |>
  mutate(
    seed = replication
  ) |>
  # remove duplicated species
  filter(sp1 != sp2) |>
  write_csv('simulation_pars.csv')


# write slurm bash file
bash_file <- paste0('#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 1-0:00:00
#SBATCH --mem-per-cpu=1500M
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
