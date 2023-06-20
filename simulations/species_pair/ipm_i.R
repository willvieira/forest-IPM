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
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# csv variables
pars <- read_csv('simulation_pars.csv')[array_id, ]
# seed
set.seed(pars$seed)
# species_id
sp1 <- pars$sp1
sp2 <- pars$sp2

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
pars_sp2 <- getPars_sp(
  sp = sp2,
  method = pars$param_method,
  path = file.path('..', 'data', 'parameters')
)

# add missing parameter that we do not expore (location of time to reproduce)
pars_sp1[['Rep']] = c(Loc = 0)
pars_sp2[['Rep']] = c(Loc = 0)

# initial size distribution of focal species
N_sp1 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = pars$Ninit_sp1
)

# Size distribution of competition species
N_sp2 <- init_pop(
  params = pars_sp2,
  L = 127,
  h = 1,
  N = pars$Ninit_sp2
)

# Prepare initial state and matrix to save size at each time step
ntmat_sp1 <- matrix(0, length(N_sp1$meshpts), pars$n_time)
ntmat_sp1[, 1] <- N_sp1$Nvec

ntmat_sp2 <- matrix(0, length(N_sp2$meshpts), pars$n_time)
ntmat_sp2[, 1] <- N_sp2$Nvec

# define mean optimal temp and prec based on focal species
temp <- mean(c(pars_sp1[['growth']]['optimal_temp'], pars_sp1[['mort']]['optimal_temp']))
prec <- mean(c(pars_sp1[['growth']]['optimal_prec'], pars_sp1[['mort']]['optimal_prec']))

# Get kernel projection for initial time
K0_sp1 = mkKernel(
  Nvec_intra = N_sp1,
  Nvec_inter = N_sp2,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp1,
  plot_random = rep(pars$plot_random, 3)
)

K0_sp2 = mkKernel(
  Nvec_intra = N_sp2,
  Nvec_inter = N_sp1,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp2,
  plot_random = rep(pars$plot_random, 3)
)

nochange = 0

# time dynamic
year_summ = data.frame(
  year = 1,
  N_sp1 = sum(N_sp1$Nvec),
  N_sp2 = sum(N_sp2$Nvec),
  BA_sp1 = size_to_BAplot(N_sp1, plot_size = pars$plotSize),
  BA_sp2 = size_to_BAplot(N_sp2, plot_size = pars$plotSize)#,
  # lambda_sp1 = max(Re(eigen(K0_sp1$K)$values)),
  # lambda_sp2 = max(Re(eigen(K0_sp2$K)$values))
)

for(t in 2:pars$n_time)
{
  # update the state and save output
  ntmat_sp1[, t] <- N_sp1$Nvec <- K0_sp1$K %*% ntmat_sp1[, t - 1]
  ntmat_sp2[, t] <- N_sp2$Nvec <- K0_sp2$K %*% ntmat_sp2[, t - 1]

  # update the kernel
  K0_sp1 <- mkKernel(
    Nvec_intra = N_sp1,
    Nvec_inter = N_sp2,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp1,
    plot_random = rep(pars$plot_random, 3)
  )

  K0_sp2 = mkKernel(
    Nvec_intra = N_sp2,
    Nvec_inter = N_sp1,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp2,
    plot_random = rep(pars$plot_random, 3)
  )

  # time dynamic
  year_summ = rbind(
    year_summ,
    data.frame(
      year = t,
      N_sp1 = sum(N_sp1$Nvec),
      N_sp2 = sum(N_sp2$Nvec),
      BA_sp1 = size_to_BAplot(N_sp1, plot_size = pars$plotSize),
      BA_sp2 = size_to_BAplot(N_sp2, plot_size = pars$plotSize)#,
      # lambda_sp1 = max(Re(eigen(K0_sp1$K)$values)),
      # lambda_sp2 = max(Re(eigen(K0_sp2$K)$values))
    )
  )

  if(
    abs(year_summ$N_sp1[t] - year_summ$N_sp1[t - 1]) < 1e-4 &
    abs(year_summ$N_sp2[t] - year_summ$N_sp2[t - 1]) < 1e-4
  ) nochange = nochange + 1
  if(nochange >= 10) break;

  cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')
}

ntmat_sp1 = ntmat_sp1[, 1:t]
ntmat_sp2 = ntmat_sp2[, 1:t]

# save output
saveRDS(
  list(
    pars_sp1 = pars_sp1,
    pars_sp2 = pars_sp2,
    # ntmat = ntmat,
    year_summ = year_summ,
    sim_time = Sys.time() - init_time
  ),
  paste0('output/sim_', array_id, '.RDS')
)
