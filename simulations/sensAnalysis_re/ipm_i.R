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

init_time <- Sys.time()

library(tidyverse)


# Read parameters
######################################################

# Env variables
toRun <- readRDS('toRun.RDS')
array_id <- toRun[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
# csv variables
pars <- read_csv('simulation_pars.csv')[array_id, ]
# seed
set.seed(pars$seed)
# species_id
sp1 <- pars$species_id

# load IPM funtions
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')


# load parameters (mean or random pick)
pars_sp1 <- getPars_sp(
  sp = sp1,
  method = pars$param_method,
  path = file.path('..', 'data', 'parameters')
)

# add missing parameter that we do not expore (location of time to reproduce)
pars_sp1[['Rep']] = c(Loc = 0)

# generate random effects
pars_sp1[['random_effects']] <- rnorm(
  n = 3, # growth, mort, and recruitment
  mean = 0,
  sd = c(
    pars_sp1[['growth']]['sigma_plot'],
    pars_sp1[['mort']]['sigma_plot'],
    pars_sp1[['rec']]['sigma_plot']
  )
)

# init pop N with correct mshpoints for focal and competition species
N0_sp1 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = 0
)
N0_sp2 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = 0
)

# Feed with predefined N vectors
N_init_comp <- readRDS('N_init_comp.RDS')
N0_sp1$Nvec[1:length(N_init_comp$N_init$Nvec)] <- N_init_comp$N_init$Nvec
if(pars$comp == 'high')
  N0_sp2$Nvec[1:length(N_init_comp$N_comp_high$Nvec)] <- N_init_comp$N_comp_high$Nvec


# define clim
if(pars$clim == 'good') {
  temp <- mean(c(pars_sp1[['growth']]['optimal_temp'], pars_sp1[['mort']]['optimal_temp']))
  prec <- mean(c(pars_sp1[['growth']]['optimal_prec'], pars_sp1[['mort']]['optimal_prec']))
}else{
  clim_suboptimal <- read_csv('climEffect_suboptimal_clim.csv') |>
    filter(species_id == sp1)
  temp <- clim_suboptimal$temp
  prec <- clim_suboptimal$prec
}

# Get kernel projection for initial time
K0_sp1 = mkKernel(
  Nvec_intra = N0_sp1,
  Nvec_inter = N0_sp2,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp1,
  plot_random = pars_sp1[['random_effects']]
)

lambda <- max(Re(eigen(K0_sp1$K)$values))

# save output
saveRDS(
  list(
    pars_sp1 = pars_sp1,
    lambda = lambda,
    sim_time = Sys.time() - init_time
  ),
  paste0('output/sim_', array_id, '.RDS')
)
