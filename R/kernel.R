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

#' Build the IPM Kernel
#'
#' Assembles the full Integral Projection Model kernel (K = P + F) for a
#' focal species given its size distribution vectors and environmental conditions.
#'
#' @param Nvec_intra List. Intraspecific size distribution; output of \code{init_pop()}.
#' @param Nvec_inter List. Interspecific size distribution; output of \code{init_pop()}.
#' @param delta_time Numeric. Time step in years.
#' @param plotSize Numeric. Plot area in square meters.
#' @param Temp Numeric. Mean annual temperature (degrees Celsius).
#' @param Prec Numeric. Mean annual precipitation (mm).
#' @param pars Named list. Species-specific parameters with elements
#'   \code{growth}, \code{mort}, \code{rec}, and \code{sizeIngrowth}.
#' @param plot_random Numeric vector of length 3. Plot-level random effects
#'   for growth (1), mortality (2), and recruitment (3).
#'
#' @return A named list with elements \code{K} (full kernel), \code{P} (growth
#'   x survival kernel), and \code{F} (recruitment kernel), each a square matrix
#'   of dimension equal to the number of mesh points.
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



#' Generate Initial Population Size Distribution
#'
#' Generates mesh points and a smooth initial size distribution scaled to an
#' expected total population size \code{N}, using a log-normal distribution.
#'
#' @param params Named list. Species-specific parameters; must contain
#'   \code{params[["growth"]]["Lmax"]}.
#' @param L Numeric. Lower boundary of the kernel (minimum size in mm).
#' @param h Numeric. Integration bin width in mm.
#' @param N Numeric. Approximated target total population size (number of individuals).
#' @param accuracy Numeric. Minimum accepted absolute error between realized
#'   and target \code{N}. Default is \code{0.001}.
#' @param meanSize Numeric. Mean of the log-normal size distribution in natural
#'   scale (mm). Default is \code{130}.
#' @param sdSize Numeric. Standard deviation of the log-normal size distribution
#'   in natural scale. Default is \code{1.8}.
#'
#' @return A named list with elements \code{meshpts} (numeric vector of mesh
#'   points), \code{Nvec} (numeric vector of individual counts per mesh point),
#'   and \code{h} (integration bin width).
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

    fct <- 1
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

    N_out <- dbh_den * fct
  }

  return(
    list(
      meshpts = msh,
      Nvec = N_out,
      h = h
    )
  )
}
