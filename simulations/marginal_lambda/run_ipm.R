# Script to run simulation on cluster

set.seed(0.0)

library(tidyverse)

spIds <- read_csv(
  paste0(readLines('_data.path'), 'species_id.csv')
)

# create and export csv sim parameters
sim_clim <- expand_grid(
  species_id = spIds$species_id_old,
  clim = seq(0, 1, 0.1),
  var = c('temp', 'prec'),
  comp = c('low', 'high'),
  replication = 1:50
) |>
mutate(arrayID = 1:n())

sim_comp <- expand_grid(
  species_id = spIds$species_id_old,
  var = c('cons', 'het'),
  comp = seq(0, 80, 5),  
  replication = 1:50
) |>
mutate(arrayID = nrow(sim_clim) + (1:n()))

saveRDS(
  list(
    'climate' = sim_clim,
    'competition' = sim_comp
  ),
  'simulation_pars.RDS'
)

sim_pars <- tibble(
    sim = c(rep('climate', nrow(sim_clim)), rep('competition', nrow(sim_comp)))
  ) |>
  mutate(
  deltaTime = 1,
  plotSize = 180,
  param_method = 'random'
  ) |>
  write_csv('simulation_metapars.csv')


dir.create('output')

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
