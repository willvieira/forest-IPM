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
sp <- pars$species_id

# load IPM funtions
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')


# load parameters (mean or random pick)
pars_sp <- getPars_sp(
  sp = sp,
  method = pars$param_method,
  path = file.path('..', 'data', 'parameters')
)

# add missing parameter that we do not expore (location of time to reproduce)
pars_sp[['Rep']] = c(Loc = 0)


# kernel parameters
mshpts_h <- get_mesh(
  L = 127,                                    # lower value of matrix
  U = round(pars_sp[['growth']]['Lmax'], 0),  # upper value of matrix
  h = 1                                       # bin size
)


# initial size distribution of focal species
N_intra <- init_pop(
  params = pars_sp,
  mshpts = mshpts_h,
  expected_N = pars$N_init
)

# Size distribution of competition species
N_inter <- init_pop(
  params = pars_sp,
  mshpts = mshpts_h,
  expected_N = pars$N_inter
)

# Prepare initial state and matrix to save size at each time step
ntmat <- matrix(0, length(N_intra), pars$n_time)
ntmat[, 1] <- N_intra

# Get kernel projection for initial time
K0 = mkKernel(
  meshpoints = mshpts_h,
  Nvec_intra = N_intra,
  Nvec_inter = N_inter,
  delta_time = pars$deltaTime,
  plot_size = 399,
  Temp = mean(
    c(pars_sp[['growth']]['optimal_temp'], pars_sp[['mort']]['optimal_temp'])
  ),
  Prec = mean(
    c(pars_sp[['growth']]['optimal_prec'], pars_sp[['mort']]['optimal_prec'])
  ),
  pars = pars_sp,
  plot_random = rep(pars$plot_random, 3)
)

# compute pop growth rate
lambda <- max(Re(eigen(K0$K)$values))

nochange = 0

# time dynamic
for(t in 2:pars$n_time)
{
  # update the state
  ntmat[, t] <- K0$K %*% ntmat[, t - 1]

  # update the kernel
  K0 <- mkKernel(
    meshpoints = mshpts_h,
    Nvec_intra = ntmat[, t],
    Nvec_inter = N_inter,
    delta_time = pars$deltaTime,
    plot_size = 399,
    Temp = mean(
      c(pars_sp[['growth']]['optimal_temp'], pars_sp[['mort']]['optimal_temp'])
    ),
    Prec = mean(
      c(pars_sp[['growth']]['optimal_prec'], pars_sp[['mort']]['optimal_prec'])
    ),
    pars = pars_sp,
    plot_random = rep(pars$plot_random, 3)
  )

  # calculate pop growth rate
  lambda[t] = max(Re(eigen(K0$K)$values))

  if(abs(sum(ntmat[, t]) - sum(ntmat[, t-1])) < 1e-7) nochange = nochange + 1
  if(nochange >= 10) break;

  cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')
}

ntmat = ntmat[, 1:t]

# save output
saveRDS(
  list(
    ntmat = ntmat,
    year_summ = ntmat |>
      reshape2::melt() |>
      group_by(Var2) |>
      summarise(
        BA_plot = BAind_to_BAplot(
          size_to_BAind(mshpts_h$meshpts), value, plot_size = 399
        ),
        N = sum(value)
      ) |>
      rename(year = Var2) |>
      bind_cols(
        lambda = lambda,
        BA_inter = BAind_to_BAplot(
          size_to_BAind(mshpts_h$meshpts), N_inter, plot_size = 399
        )
      ),
    sim_time = Sys.time() - init_time
  ),
  paste0('output/sim_', array_id, '.RDS')
)
