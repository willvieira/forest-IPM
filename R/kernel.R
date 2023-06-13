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


# Probability functions for growth, mortality, and recruitment
#################################################################

# Probability function for growth
vonBertalanffy_lk = function(
  pars, delta_time, size_t1, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random
){
  # get the mean of the prob function for each "ind"
  growth_mean = vonBertalanffy_f(pars, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random)

  # likelihood of increment y given defined model and parameters
  growth = dnorm(size_t1, mean = growth_mean, sd = pars['sigma_obs'])

  return( growth )
}



## Survival x growth kernel
P_xEC = function(
  size_t1, size_t0, delta_time, BA_comp_intra, BA_comp_inter, Temp, Prec, parsGrowth, parsMort, plot_random
){
  pkernel = vonBertalanffy_lk(
    parsGrowth, delta_time, size_t1, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random[1]
  ) *
  survival_f(
    parsMort, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random[2]
  )

  return( pkernel )
}


# Probability function for ingrowth
ingrowth_lk <- function(
  size_t1, size_t0, delta_time, plot_size, BA_adult_sp, BA_adult, parsIngrowth, parsSizeIngrowth, parsRep, plot_random
){
  ingrowth_prob = ingrowth_f(
    parsIngrowth, delta_time, plot_size, BA_adult_sp, BA_adult, plot_random
  ) *
  truncnorm::dtruncnorm(
    size_t1,
    a = 127,
    b = Inf,
    mean = parsSizeIngrowth['size_int'] + parsSizeIngrowth['phi_time'] * delta_time,
    sd = parsSizeIngrowth['sigma_size']
  ) *
  (0 + ((1)/(1 + exp(-.5 * (size_t0 - parsRep['Loc'])))^(1)))

  return( ingrowth_prob )
}



# Kernel and functions around building the kernel
#################################################################

# Compute mesh points
get_mesh <- function(L, U, h = 1)
{
  # mesh points
  m <- length(seq(L, U, h))
  meshpts <- L + ((1:m) - 1/2) * h

  return( list(meshpts = meshpts, h = h))
}

# Full Kernel
mkKernel = function(
  meshpoints,
  Nvec_intra, Nvec_inter,
  delta_time, plot_size, Temp, Prec, pars,
  plot_random # vector of [1] growth, [2] mortality, and [3] recruitment ofsets
){
  meshpts = meshpoints$meshpts
  h = meshpoints$h

  # compute competition metrics
  BA_comp_intra <- size_to_BAcomp(meshpts, Nvec_intra, plot_size)
  BA_comp_inter <- size_to_BAcomp(meshpts, Nvec_inter, plot_size)
  BA_adult_sp <- BAind_to_BAplot(size_to_BAind(meshpts), Nvec_intra, plot_size)
  BA_adult <- BAind_to_BAplot(size_to_BAind(meshpts), Nvec_intra + Nvec_inter, plot_size)

  P <- h * outer(
    meshpts, meshpts,
    P_xEC,
    delta_time, BA_comp_intra, BA_comp_inter, Temp, Prec,
    pars[['growth']], pars[['mort']], plot_random
  )
  F <- h * outer(
    meshpts, meshpts,
    ingrowth_lk,
    delta_time, plot_size, BA_adult_sp, BA_adult, pars[['rec']], pars[['sizeIngrowth']], pars[['Rep']], plot_random[3]
  )

  K <- P + F

  return(list(K = K, meshpts = meshpts, P = P, F = F))
}



# function to get parameters for specific species_id
#' sp: species_id
#' method: either `mean` or `random`. `Mean` is the posterior mean and `random` is a single draw from the posterior mean
#' path: default is set to 'data/parameters'
getPars_sp <- function(sp, method, path = NULL)
{
  files_to_load <- dir(
    ifelse(
      is.null(path),
      file.path('data', 'parameters'),
      path
    ),
    pattern = sp,
    full.names = TRUE
  )

  # remove  random effect files
   files_to_load <- files_to_load[
    grep('meanRandomEffect', files_to_load, invert = TRUE)
   ]

  # read pars
  pars <- map(
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

  # export pars
  if(method == 'mean') {
    pars |>
      map(~apply(.x, 2, mean))
  }else if(method == 'random') {
    pars |>
      map(
        ~ .x[sample(1:nrow(.x), 1), ] |>
          pivot_longer(cols = everything()) |>
          deframe()
      )
  }else{
    stop('`param_method` parameter must be either `mean` or `random`')
  }
}
