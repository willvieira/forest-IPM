##############################################################
# Kernel function with its respective growth*survival kernels
# Will Vieira
# January 20, 2019
# Last updated: March, 18, 2023
##############################################################



######################
 # 1. Growth
 # 2. Survival (1 - mortality)
 # 3. Growth * survival kernel
 # 3. Full kernel
######################



# Probability function for growth
vonBertalanffy_lk = function(
  pars, delta_time, size_t1, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, randomEffects
){
  # get the mean of the prob function for each "ind"
  growth_mean = vonBertalanffy_f(pars, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, randomEffects)

  # likelihood of increment y given defined model and parameters
  growth = dnorm(size_t1, mean = growth_mean, sd = pars['sigma_obs'])

  return( growth )
}



## Survival x growth kernel
P_xEC = function(
  size_t1, size_t0, delta_time, BA_comp_intra, BA_comp_inter, Temp, Prec, parsGrowth, parsMort, randomEffects
){
  pkernel = vonBertalanffy_lk(
    parsGrowth, delta_time, size_t1, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, randomEffects
  ) *
  survival_f(
    parsMort, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, randomEffects
  )

  return( pkernel )
}


# Probability function for ingrowth
ingrowth_lk <- function(
  size_t1, size_t0, delta_time, plot_size, BA_adult_sp, BA_adult, parsIngrowth, parsSizeIngrowth, randomEffects
){
  ingrowth_prob = ingrowth_f(
    parsIngrowth, delta_time, plot_size, BA_adult_sp, BA_adult, randomEffects
  ) *
  truncnorm::dtruncnorm(
    size_t1,
    a = 127,
    b = Inf,
    mean = parsSizeIngrowth['size_int'] + parsSizeIngrowth['phi_time'] * delta_time,
    sd = parsSizeIngrowth['sigma_size']
  )

  return( ingrowth_prob )
}


# Full Kernel
mkKernel = function(
  m, U, L,
  delta_time, plot_size, BA_comp_intra, BA_comp_inter, BA_adult_sp, BA_adult, Temp, Prec,
  pars, randomEffects
){
  # mesh points
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h

  P <- h * outer(
    meshpts, meshpts,
    P_xEC,
    delta_time, BA_comp_intra, BA_comp_inter, Temp, Prec,
    pars[['growth']], pars[['mort']], randomEffects
  )
  F <- h * outer(
    meshpts, meshpts,
    ingrowth_lk,
    delta_time, plot_size, BA_adult_sp, BA_adult, pars[['rec']], pars[['sizeIngrowth']], randomEffects
  )

  K <- P + F

  return(list(K = K, meshpts = meshpts, P = P, F = F))
}



# function to get parameters for specific species_id
getPars_sp <- function(sp)
{
  files_to_load <- dir(
    file.path('data', 'parameters'),
    pattern = sp,
    full.names = TRUE
  )

  map(
    setNames(
      files_to_load,
      gsub(
        paste0('data/parameters/|_|', sp, '|.csv'),
        '',
        files_to_load
      )
    ),
    ~ read_csv(.x, progress = FALSE, show_col_types = FALSE)
  )
}
