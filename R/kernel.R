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

# Full Kernel
#' @export
mkKernel = function(
  Nvec_intra, Nvec_inter,
  delta_time, plotSize, Temp, Prec, pars,
  plot_random # vector of [1] growth, [2] mortality, and [3] recruitment ofsets
){
  meshpts = Nvec_intra$meshpts
  h = Nvec_intra$h

  # compute competition metrics
  BA_comp_intra <- size_to_BAcomp(
    N_intra = Nvec_intra,
    plot_size = plotSize
  )
  BA_comp_inter <- size_to_BAcomp(
    N_intra = Nvec_intra,
    N_inter = Nvec_inter,
    plot_size = plotSize
  )
  BAplot_intra <- size_to_BAplot(
    N = Nvec_intra,
    plot_size = plotSize
  )
  BAplot_inter <- size_to_BAplot(
    N = Nvec_inter,
    plot_size = plotSize
  )
  BAplot_total <- BAplot_intra + BAplot_inter

  # kernel
  P <- h * outer(
    meshpts, meshpts,
    P_xEC,
    delta_time, BA_comp_intra, BA_comp_inter, Temp, Prec,
    pars[['growth']], pars[['mort']], plot_random
  )
  F <- h * outer(
    meshpts, meshpts,
    ingrowth_lk,
    delta_time, plotSize, BAplot_intra, BAplot_total, pars[['rec']], pars[['sizeIngrowth']], pars[['Rep']], plot_random[3]
  )

  K <- P + F

  return(list(K = K, P = P, F = F))
}



# function to get parameters for specific species_id
#' sp: species_id
#' method: either `mean` or `random`. `Mean` is the posterior mean and `random` is a single draw from the posterior mean
#' path: default is set to 'data/parameters'
#' @export
getPars_sp <- function(sp, method, path = NULL)
{
  pars_path <- ifelse(
    is.null(path),
    file.path('data', 'parameters'),
    path
  )

  files_to_load <- dir(
    pars_path,
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
          paste0(pars_path, '/|_|', sp, '|.csv'),
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



# Function to generate meshpoints and a smooth initial size distribution in
# function of expected population size N
#' params: species specific parameters
#' L: lower boundary of kernel
#' h: intgration bin
#' N: approximated population size N expected for the initial size dist
#' accuracy: minimum accepted error difference between N and expected N
#' meanSize: mean of lognormal distribution for size distribution
#' sdSize: standard deviation of lognormal distribution for size distribution
#' Both meanSize and sdSize parameters are in natural scale
#' @export
init_pop <- function(
  params,
  L,
  h,
  N,
  accuracy = 0.001,
  meanSize = 130,
  sdSize = 1.8
){
  # Compute mesh points
  Lmax <- round(params[['growth']]['Lmax'], 0)
  m <- length(seq(L, Lmax, h))
  msh <- L + ((1:m) - 1/2) * h

  if(N < 0)
  {
    stop('Argument `N` must be larger or equal than zero.')
  }
  # In case `N` equal zero, return empty dist vector
  else if(N == 0)
  {
    N_out <- rep(0, length(msh))
  }
  # Generate smooth size dist in function of N
  else
  {
    # generate random individuals from the lognorm distribution
    dbh <- qlnorm(
      runif(
        n = 1e4,
        min = plnorm(min(msh), log(meanSize), log(sdSize)),
        max = plnorm(
          Lmax, log(meanSize), log(sdSize)
        )
      ),
      log(meanSize), log(sdSize)
    )

    # get density distribution from genereted individual sizes
    dbh_den <- density(dbh, n = length(msh))$y

    # transform density distribution to approximate total pop size to
    # the expected N argument
    diff_N <- N - sum(dbh_den)

    Min = 0; Max = N * 5
    while(abs(diff_N) > accuracy) {
      fct = runif(1, Min, Max)
      new_dbh_den <- dbh_den * fct
      diff_N <- N - sum(new_dbh_den)
      if(diff_N > 0) {
        Min = fct
      }else{
        Max = fct
      }
    }

    # In case the generated size dist is already close to the expected N
    if(exists('fct')) {
      N_out <- dbh_den * fct
    }else{
      N_out <- dbh_den
    }
  }

  return(
    list(
      meshpts = msh,
      Nvec = N_out,
      h = h
    )
  )
}
