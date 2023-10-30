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


# prepare parameters
pars |>
  select(contains('.')) |>
  pivot_longer(cols = everything()) |>
  mutate(
    vr = str_replace(name, '\\..*', ''),
    par = str_replace(name, paste0(vr, '.'), '')
  ) |>
  select(!name) |>
  group_split(vr) %>%
  # fuck tidyverse for not wanting to implement an argument to keep a named list
  set_names(map_chr(., ~.x$vr[1])) |>
  map(
    ~.x |>
      select(!vr) |>
      pivot_wider(names_from = par) |>
      as_vector()
  ) ->
pars_sp


# Load initial conditions
N_init_ls <- readRDS(paste0('N_dist/N_', sp, '.RDS'))

# deal with competition
if(metapars$sim == 'pars')
{
  N_init = N_init_ls[['N_init']]
  N_comp = N_init_ls[['N_comp']][['1_1']]

}else if(metapars$sim == 'temp')
{
  N_init = N_init_ls[['N_init']]
  N_comp = N_init_ls[['N_comp']][['1_1']]

}else
{
  N_init = N_init_ls[['N_init']]
  N_comp = N_init_ls[['N_comp']][[paste0(pars$BA_sd, '_', pars$replication)]]
}

# deal with climate
temp <- pars$temp
prec <- mean(c(
  pars_sp[['growth']]['optimal_prec'],
  pars_sp[['mort']]['optimal_prec'],
  pars_sp[['rec']]['optimal_prec']
))



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

writeLines(
  as.character(lambda),
  paste0('output/sim_', array_id, '.csv')
)
