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
  size_t1, size_t0, delta_time, plot_size, BA_adult_sp, BA_adult, Temp, Prec, parsIngrowth, parsSizeIngrowth, plot_random
){
  ingrowth_prob = ingrowth_f(
    parsIngrowth, delta_time, plot_size, BA_adult_sp, BA_adult, Temp, Prec, plot_random
  ) *
  truncnorm::dtruncnorm(
    size_t1,
    a = 127,
    b = Inf,
    mean = parsSizeIngrowth['size_int'] + parsSizeIngrowth['phi_time'] * delta_time,
    sd = parsSizeIngrowth['sigma_size']
  ) #*
  # 
  # (0 + ((1)/(1 + exp(-.5 * (size_t0 - parsRep['Loc'])))^(1)))

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
    delta_time, plotSize, BAplot_intra, BAplot_total, Temp, Prec, pars[['rec']], pars[['sizeIngrowth']], plot_random[3]
  )

  K <- P + F

  return(list(K = K, P = P, F = F))
}



# function to get parameters for specific species_id
#' sp: species_id
#' method: either `mean` or `random`. `Mean` is the posterior mean and `random` is a single draw from the posterior mean
#' model: name of specific demographic model. Can be either a character or a vector of size 3 with the name of the model for each of the [1] growth, [2] mort, or [3] recruitment model. Accepted values are `intcpt`, `plot`, `comp`, and `clim` that can be used in incremental form. Defaut set to the complete model: `intcpt_plot_comp_clim`.
#' path: default is set to 'data/output_sim_processed'
#' @export
getPars_sp <- function(
  sp,
  method,
  model = 'intcpt_plot_comp_clim',
  path = file.path('data', 'output_sim_processed')
){
  # deal for `model` parameter
  if(length(model) == 1) {
    model = rep(model, 3)
  }else if(length(model) != 3) {
    stop('`method` parameter must be a character of length 1 or 3.')
  }

  for(MD in model)
    if(!MD %in% c('intcpt', 'intcpt_plot', 'intcpt_plot_comp', 'intcpt_plot_comp_clim')) {
      stop(
        paste0('model parameter should be one of the following character:\n', paste0(paste0('-', c('intcpt', 'intcpt_plot', 'intcpt_plot_comp', 'intcpt_plot_comp_clim')), collapse='\n'))
      )
    }
  
  # define full path to load parameters for each vital rate
  pars_dir <- setNames(
    paste0(
      path, '/', c('growth', 'mort', 'recruit', 'sizeIngrowth'), '/', c(model, 'time_truc'), '/posterior_pop_', sp, '.RDS'
    ),
    c('growth', 'mort', 'rec', 'sizeIngrowth')
  )

  # export pars
  if(method == 'mean') {
    pars_dir |>
      map(
        ~ readRDS(.x) |>
          filter(par != 'lp__') |>
          group_by(par) |>
          reframe(value = mean(value)) |>
          pivot_wider(names_from = 'par') |>
          as_vector()
      )
  }else if(method == 'random') {
    pars_dir |>
      map(
        ~ readRDS(.x) |>
          filter(par != 'lp__') |>
          filter(iter == sample(1:4000, 1)) |>
          select(!iter) |>
          pivot_wider(names_from = 'par') |>
          as_vector()      
      )
  }else{
    stop('`param_method` parameter must be either `mean` or `random`')
  }
}

# Function to tranform rowwise tibble of parameters in list format
# Which is the current way of passing parameters to the IPM
pars_to_list <- function(pars)
{
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
    )
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
