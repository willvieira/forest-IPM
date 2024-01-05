# Script to run linear model with varying dispersal

library(tidyverse)
library(cmdstanr)


# load batch info
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# load simulation pars for specific run
pars <- readRDS('simulation_pars.RDS')[array_id, ]

# load data
sim_name_i = paste(pars$species_id, pars$border, pars$sim, pars$cond, sep = '_')

data_i <- readRDS(
  paste0(
    'out_species/',
    sim_name_i,
    '.RDS'
  )
)

data_stan <- with(
  data_i,
  list(
    N = nrow(data_i),
    lambda_log = log(lambda),
    MAT = bio_01_mean
  )
)

# Compile stan model
stanModel <- cmdstan_model('model.stan')

# create output dir
dir.create(file.path('output', sim_name_i), recursive = TRUE)

# sample
md_out <- stanModel$sample(
    data = data_stan,
    parallel_chains = 4,
    # output_dir = 'file.path('output', sim_name_i)'
    output_dir = 'output'
)
