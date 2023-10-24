# Script to run simulation on cluster

set.seed(0.0)

library(tidyverse)

spIds <- read_csv(
  paste0(readLines('_data.path'), 'species_id.csv')
)

# create and export csv sim parameters
sim_pars <- expand_grid(
    species_id = spIds$species_id_old,
    comp = c('low', 'high'),
    clim = c('cold', 'center', 'hot'),
    deltaTime = 1,
    plotSize = 180,
    param_method = 'random',
    replication = 1:500
  ) |>
  mutate(
    seed = replication
  ) |>
  write_csv('simulation_pars.csv')


# Slurm returns an error for array above 10k
# loop over simulations to each batch is limited to 10k

i_max <- floor(nrow(sim_pars)/10000)
batch_max <- rep(10000, i_max)
batch_max[i_max + 1] <- nrow(sim_pars)/10

# write slurm bash file
for(i in 0:i_max)
{
  bash_file <- paste0('#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu=700M
#SBATCH --ntasks=1
#SBATCH --job-name=ipm', i, '
#SBATCH --mail-user=willian.vieira@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-', batch_max[i + 1], '

module load StdEnv/2020 r/4.1.0

BATCH=', i, ' R -f ipm_i.R
')

  # save bash file
  writeLines(bash_file, paste0('sub_', i, '.sh'))

  # run slurm
  system(paste0('sbatch sub_', i, '.sh'))
}
