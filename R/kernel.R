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
#' @param Nvec_intra List. Intraspecific size distribution.
#' @param Nvec_inter List. Interspecific size distribution.
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


