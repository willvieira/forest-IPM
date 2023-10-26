###############################################################################
# Layer 1 function to run the IPM in the cluster

# Input
# - Env variables:
#  - array_id (to define which csv simulation to run)
# - csv variables
#  - seed
#  - clim temp and prec [0-1]
#  - N inter for competiton
#  - N intra for initial population size
#  - species_id
#  - time_max of simulation
#  - deltaTime for time iteration
#  - Parameter mean or random pick
#  - Random effects (This could be either a character or a list of values)

# Output
# - List:
# -- matrix of density distribution x time
# -- data.frame with BA_plot_sp and lambda over time
################################################################################

# turn off scientific notation
options(scipen=999)

library(tidyverse)


# Read parameters
######################################################

# Env variables
batch_id <- as.numeric(Sys.getenv('BATCH'))
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + (10000 * batch_id) 
# csv variables
metapars <- read_csv('simulation_metapars.csv')[array_id, ]

pars_ls <- readRDS('simulation_pars.RDS')

pars <- pars_ls[[metapars$sim]] |>
  filter(arrayID == array_id)

# seed
set.seed(pars$arrayID)
# species_id
sp <- pars$species_id

# load IPM funtions
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')


# load parameters (mean or random pick)
pars_sp <- getPars_sp(
  sp = sp,
  method = metapars$param_method,
  model = 'intcpt_plot_comp_clim',
  path = readLines('_data.path')
)


# Load initial conditions
N_init_ls <- readRDS(paste0('N_dist/N_', metapars$sim, '_', sp, '.RDS'))

if(metapars$sim == 'climate')
{
  # deal with competition
  N_init = N_init_ls[[paste0('N_', pars$comp)]]
  N_comp = N_init
  N_comp$Nvec = 0

  # deal with climate
  if(pars$var == 'temp') {
    temp <- pars$clim
    prec <- mean(c(
      pars_sp[['growth']]['optimal_prec'],
      pars_sp[['mort']]['optimal_prec'],
      pars_sp[['rec']]['optimal_prec']
    ))
  }else{
    temp <- mean(c(
      pars_sp[['growth']]['optimal_temp'],
      pars_sp[['mort']]['optimal_temp'],
      pars_sp[['rec']]['optimal_temp']
    ))
    prec <- pars$clim
  }
}else if(metapars$sim == 'competition')
{
  # deal with competition
  if(pars$var == 'cons') {
    N_init = N_init_ls[[as.character(pars$comp)]]
    N_comp = N_init
    N_comp$Nvec = 0
  }else {
    N_init = readRDS(paste0('N_dist/N_climate_', sp, '.RDS'))[['N_low']]
    N_comp = N_init_ls[[as.character(pars$comp)]]
  }

  # deal with climate
  temp <- mean(c(
    pars_sp[['growth']]['optimal_temp'],
    pars_sp[['mort']]['optimal_temp'],
    pars_sp[['rec']]['optimal_temp']
  ))
  prec <- mean(c(
    pars_sp[['growth']]['optimal_prec'],
    pars_sp[['mort']]['optimal_prec'],
    pars_sp[['rec']]['optimal_prec']
  ))
}


# Get kernel projection for initial time
K = mkKernel(
  Nvec_intra = N_init,
  Nvec_inter = N_comp,
  delta_time = metapars$deltaTime,
  plotSize = metapars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp,
  plot_random = rep(0, 3)
)

lambda <- max(Re(eigen(K$K)$values))


# prepare saving
pars_vec <- pars_sp |> unlist()
out_vec <- c(pars_vec, lambda)


writeLines(
  paste0(out_vec, collapse = ','),
  paste0('output/sim_', array_id, '.csv')
)

# create header fo csv file
if(array_id == 1) {
  out <- c(names(pars_vec), 'lambda')
  
  writeLines(
    paste0(out, collapse = ','),
    paste0('output/sim_0.csv')
  )
}
