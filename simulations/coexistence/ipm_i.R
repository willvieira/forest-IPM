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

# init pop N with correct mshpoints for focal and competition species
N0_sp1 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = 0
)
N0_sp2 <- init_pop(
  params = pars_sp2,
  L = 127,
  h = 1,
  N = 0
)

# Feed with predefined N vectors
N_init_comp <- readRDS('N_init_comp.RDS')
N0_sp1$Nvec[1:length(N_init_comp$N_init$Nvec)] <- N_init_comp$N_init$Nvec

# define mean optimal temp and prec based on resident species
temp <- mean(c(pars_sp2[['growth']]['optimal_temp'], pars_sp2[['mort']]['optimal_temp']))
prec <- mean(c(pars_sp2[['growth']]['optimal_prec'], pars_sp2[['mort']]['optimal_prec']))

# Get kernel projection for initial time
K0_sp1 = mkKernel(
  Nvec_intra = N0_sp1,
  Nvec_inter = N0_sp2,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp1,
  plot_random = rep(0, 3)
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sub simulation 1
# Invader species i alone
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# time dynamic
year_summ = data.frame(
  year = 1,
  N_sp1 = sum(N0_sp1$Nvec),
  BA_sp1 = size_to_BAplot(N0_sp1, plot_size = pars$plotSize)
)

# Save matrix state every 5 steps
mat_sp1 <- matrix(0, nrow = pars$n_time, ncol = length(N0_sp1$Nvec))
mat_row <- 1
mat_sp1[mat_row, ] <- N0_sp1$Nvec

lambdas <- max(Re(eigen(K0_sp1$K)$values))
Ntp1_sp1 <- Nt_sp1 <- N0_sp1
Kt_sp1 <- K0_sp1

nochange = 0

for(t in 2:pars$n_time)
{
  # update the state and save output
  Ntp1_sp1$Nvec <- Kt_sp1$K %*% Nt_sp1$Nvec

  # update the kernel
  Kt_sp1 <- mkKernel(
    Nvec_intra = Ntp1_sp1,
    Nvec_inter = N0_sp2,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp1,
    plot_random = rep(0, 3)
  )

  # time dynamic
  year_summ = rbind(
    year_summ,
    data.frame(
      year = t,
      N_sp1 = sum(Ntp1_sp1$Nvec),
      BA_sp1 = size_to_BAplot(Ntp1_sp1, plot_size = pars$plotSize)
    )
  )

  if(t %% 10000 == 0) {
    mat_row <- mat_row + 1
    mat_sp1[mat_row, ] <- as.vector(Ntp1_sp1$Nvec)
    lambdas[mat_row] <- max(Re(eigen(Kt_sp1$K)$values))
  }

  if(abs(year_summ$N_sp1[t] - year_summ$N_sp1[t - 1]) < 1e-4)
    nochange = nochange + 1
  if(nochange >= 10) break;

  Nt_sp1 <- Ntp1_sp1

  cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')
}

# save last state info
mat_sp1[mat_row + 1, ] <- as.vector(Ntp1_sp1$Nvec)
mat_sp1 <- mat_sp1[1:mat_row + 1, ]
lambdas <- append(
  lambdas,
   max(Re(eigen(Kt_sp1$K)$values))
)
lambda_years <- c(1, 1:mat_row * 10)
lambda_years[mat_row + 1] <- t
year_summ$lambda <- NA
year_summ$lambda[year_summ$year %in% lambda_years] <- lambdas


year_summ_subsim1 <- year_summ
mat_sp1_subsim1  <- mat_sp1





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sub simulation 2
# Species i invading species j
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# First step is to run resident species to equilibrium
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Update resident pop size to the predefined initial state
N0_sp2$Nvec[1:length(N_init_comp$N_init$Nvec)] <- N_init_comp$N_init$Nvec
Ncomp0 <- init_pop(
  params = pars_sp2,
  L = 127,
  h = 1,
  N = 0
)

K0_sp2 = mkKernel(
  Nvec_intra = N0_sp2,
  Nvec_inter = Ncomp0,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp2,
  plot_random = rep(0, 3)
)


# time dynamic
year_summ = data.frame(
  year = 1,
  N_sp2 = sum(N0_sp2$Nvec),
  BA_sp2 = size_to_BAplot(N0_sp2, plot_size = pars$plotSize)
)

# Save matrix state every 5 steps
mat_sp2 <- matrix(0, nrow = pars$n_time, ncol = length(N0_sp2$Nvec))
mat_row <- 1
mat_sp2[mat_row, ] <- N0_sp2$Nvec

lambdas <- max(Re(eigen(K0_sp2$K)$values))
Ntp1_sp2 <- Nt_sp2 <- N0_sp2
Kt_sp2 <- K0_sp2

nochange = 0

for(t in 2:pars$n_time)
{
  # update the state and save output
  Ntp1_sp2$Nvec <- Kt_sp2$K %*% Nt_sp2$Nvec

  # update the kernel
  Kt_sp2 <- mkKernel(
    Nvec_intra = Ntp1_sp2,
    Nvec_inter = Ncomp0,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp2,
    plot_random = rep(0, 3)
  )

  # time dynamic
  year_summ = rbind(
    year_summ,
    data.frame(
      year = t,
      N_sp2 = sum(Ntp1_sp2$Nvec),
      BA_sp2 = size_to_BAplot(Ntp1_sp2, plot_size = pars$plotSize)
    )
  )

  if(t %% 10000 == 0) {
    mat_row <- mat_row + 1
    mat_sp2[mat_row, ] <- as.vector(Ntp1_sp2$Nvec)
    lambdas[mat_row] <- max(Re(eigen(Kt_sp2$K)$values))
  }

  if(abs(year_summ$N_sp2[t] - year_summ$N_sp2[t - 1]) < 1e-4)
    nochange = nochange + 1
  if(nochange >= 10) break;

  Nt_sp2 <- Ntp1_sp2

  cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')
}

# save last state info
mat_sp2[mat_row + 1, ] <- as.vector(Ntp1_sp2$Nvec)
mat_sp2 <- mat_sp2[1:mat_row + 1, ]
lambdas <- append(
  lambdas,
   max(Re(eigen(Kt_sp2$K)$values))
)
lambda_years <- c(1, 1:mat_row * 10)
lambda_years[mat_row + 1] <- t
year_summ$lambda <- NA
year_summ$lambda[year_summ$year %in% lambda_years] <- lambdas


year_summ_sp2_subsim2 <- year_summ
mat_sp2_subsim2  <- mat_sp2






# Second step is to run invader species when resident species is already in equilibrium
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get kernel projection for initial time
K0_sp1 = mkKernel(
  Nvec_intra = N0_sp1,
  Nvec_inter = Ntp1_sp2,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp1,
  plot_random = rep(0, 3)
)

K0_sp2 = mkKernel(
  Nvec_intra = Ntp1_sp2,
  Nvec_inter = N0_sp1,
  delta_time = pars$deltaTime,
  plotSize = pars$plotSize,
  Temp = temp,
  Prec = prec,
  pars = pars_sp2,
  plot_random = rep(0, 3)
)


# time dynamic
year_summ = data.frame(
  year = 1,
  N_sp1 = sum(N0_sp1$Nvec),
  N_sp2 = sum(Ntp1_sp2$Nvec),
  BA_sp1 = size_to_BAplot(N0_sp1, plot_size = pars$plotSize),
  BA_sp2 = size_to_BAplot(Ntp1_sp2, plot_size = pars$plotSize)
)

# Save matrix state every 5 steps
mat_sp1 <- matrix(0, nrow = pars$n_time, ncol = length(N0_sp1$Nvec))
mat_row <- 1
mat_sp1[mat_row, ] <- N0_sp1$Nvec
mat_sp2 <- matrix(0, nrow = pars$n_time, ncol = length(Ntp1_sp2$Nvec))
mat_row <- 1
mat_sp2[mat_row, ] <- Ntp1_sp2$Nvec

lambdas_sp1 <- max(Re(eigen(K0_sp1$K)$values))
lambdas_sp2 <- max(Re(eigen(K0_sp2$K)$values))

Ntp1_sp1 <- Nt_sp1 <- N0_sp1
Kt_sp1 <- K0_sp1
Ntp1_sp2 <- Nt_sp2 <- Ntp1_sp2
Kt_sp2 <- K0_sp2

nochange = 0

for(t in 2:pars$n_time)
{
  # update the state and save output
  Ntp1_sp1$Nvec <- Kt_sp1$K %*% Nt_sp1$Nvec
  Ntp1_sp2$Nvec <- Kt_sp2$K %*% Nt_sp2$Nvec

  # update the kernel
  Kt_sp1 <- mkKernel(
    Nvec_intra = Ntp1_sp1,
    Nvec_inter = Ntp1_sp2,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp1,
    plot_random = rep(0, 3)
  )
  Kt_sp2 <- mkKernel(
    Nvec_intra = Ntp1_sp2,
    Nvec_inter = Ntp1_sp1,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp,
    Prec = prec,
    pars = pars_sp2,
    plot_random = rep(0, 3)
  )

  # time dynamic
  year_summ = rbind(
    year_summ,
    data.frame(
      year = t,
      N_sp1 = sum(Ntp1_sp1$Nvec),
      N_sp2 = sum(Ntp1_sp2$Nvec),
      BA_sp1 = size_to_BAplot(Ntp1_sp1, plot_size = pars$plotSize),
      BA_sp2 = size_to_BAplot(Ntp1_sp2, plot_size = pars$plotSize)
    )
  )

  if(t %% 10000 == 0) {
    mat_row <- mat_row + 1
    mat_sp1[mat_row, ] <- as.vector(Ntp1_sp1$Nvec)
    lambdas_sp1[mat_row] <- max(Re(eigen(Kt_sp1$K)$values))
    mat_sp2[mat_row, ] <- as.vector(Ntp1_sp2$Nvec)
    lambdas_sp2[mat_row] <- max(Re(eigen(Kt_sp2$K)$values))
  }

  if(
    abs(year_summ$N_sp1[t] - year_summ$N_sp1[t - 1]) < 1e-4 &
    abs(year_summ$N_sp2[t] - year_summ$N_sp2[t - 1]) < 1e-4
  ) nochange = nochange + 1
  if(nochange >= 10) break;

  Nt_sp1 <- Ntp1_sp1
  Nt_sp2 <- Ntp1_sp2

  cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')
}

# save last state info
mat_sp1[mat_row + 1, ] <- as.vector(Ntp1_sp1$Nvec)
mat_sp1 <- mat_sp1[1:mat_row + 1, ]
lambdas_sp1 <- append(
  lambdas_sp1,
   max(Re(eigen(Kt_sp1$K)$values))
)
mat_sp2[mat_row + 1, ] <- as.vector(Ntp1_sp2$Nvec)
mat_sp2 <- mat_sp2[1:mat_row + 1, ]
lambdas_sp2 <- append(
  lambdas_sp2,
   max(Re(eigen(Kt_sp2$K)$values))
)
lambda_years <- c(1, 1:mat_row * 10)
lambda_years[mat_row + 1] <- t
year_summ$lambda_sp1 <- year_summ$lambda_sp2 <- NA
year_summ$lambda_sp1[year_summ$year %in% lambda_years] <- lambdas_sp1
year_summ$lambda_sp2[year_summ$year %in% lambda_years] <- lambdas_sp2


year_summ_sp1sp2_subsim2 <- year_summ
mat_sp1sp2_subsim2  <- mat_sp2


# save output
saveRDS(
  list(
    pars_sp1 = pars_sp1,
    pars_sp2 = pars_sp2,
    year_summ_subsim1 = year_summ_subsim1,
    mat_sp1_subsim1 = mat_sp1_subsim1,
    year_summ_sp2_subsim2 = year_summ_sp2_subsim2,
    mat_sp2_subsim2 = mat_sp2_subsim2,
    year_summ_sp1sp2_subsim2 = year_summ_sp1sp2_subsim2,
    mat_sp1sp2_subsim2 = mat_sp1sp2_subsim2,
    sim_time = Sys.time() - init_time
  ),
  paste0('output/sim_', array_id, '.RDS')
)
