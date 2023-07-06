# Script to run simulation on cluster

set.seed(0.0)

library(tidyverse)

spIds <- unique(
  gsub(
    '.*_|.csv', '', dir(file.path('..', 'data', 'parameters'))
  )
)

# create and export csv sim parameters
sim_pars <- expand_grid(
    species_id = spIds,
    comp = c('low', 'high'),
    clim = c('bad', 'good'),
    n_time = 2000,
    deltaTime = 1,
    plotSize = 180,
    param_method = 'random',
    replication = 1:1000
  ) |>
  mutate(
    seed = replication
  ) |>
  write_csv('simulation_pars.csv')


# generate focal and competition N vectors (low and high comp)
# so each replication has the same initial state and competition values
# There is an issue that each replicate will have a different Lmax and therefore
# a different N vector. To work around this issue, I generate a init and comp
# vector with the lowest Lmax among spIds and iterations (422) and for each
# sim run I will complete the remaining N vector with zeros
source('R/kernel.R')
pars_fake <- list('growth' = c('Lmax' = 422))
N_init <- init_pop(
  pars_fake,
  L = 127,
  h = 1,
  N = 0.5
)

N_comp_high = init_pop(
  pars_fake,
  L = 127,
  h = 1,
  N = 30
)

saveRDS(
  list(
    'N_init' = N_init,
    'N_comp_high' = N_comp_high
  ),
  'N_init_comp.RDS'
)


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
#SBATCH -t 2:00:00
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
