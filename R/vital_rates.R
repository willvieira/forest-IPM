##############################################################
# Vital rate functions for growth, mortality and reproduction
# Will Vieira
# January 20, 2019
# Last updated: March, 18, 2023
##############################################################



######################
 # 1. Growth
 # 2. Mortality
 # 3. Ingrowth
 ######################



# 1. Von Bertalanffy model to predict future size_t1 as a function of:
# - time interval between size_t0 and size_t1
# - Size_t0 (dbh in mm)
# - Competition (intra and interspecific basal area of indivudal higher than the focal individual)
# - Mean annual temperature
# - Mean annual precipitation
# - plot random effects
# - individual random effects

vonBertalanffy_f <- function(
  pars, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random
){
  # Compute r
  rPlotInd = exp(
    pars['r'] + # intercept
    plot_random +
    pars['Beta'] * (BA_comp_intra + pars['theta'] * BA_comp_inter) + # Comp
    -pars['tau_temp'] * (Temp - pars['optimal_temp'])^2 + # temp effect
    -pars['tau_prec'] * (Prec - pars['optimal_prec'])^2 # prec effect
  )

# pre calculate component of the model
rPlotTime = exp(-rPlotInd * delta_time)

# mean
mu_obs = size_t0 *
  rPlotTime +
  pars['Lmax'] * (1 - rPlotTime)

return( mu_obs )
}
 


# 2. Probability of mortality (Surv = 1 - mort) as a function of:
# - time interval between size_t0 and size_t1
# - Size_t0 (dbh in mm) - /!\ deprecated /!\
# - Competition (intra and interspecific basal area of indivudal higher than the focal individual)
# - Mean annual temperature
# - Mean annual precipitation
# - plot random effects
# - year random effects
survival_f = function(
  pars, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random
){
  # longevity rate
  longev_log <- 1/(1 + exp(
      -(
        pars['psi'] +
        plot_random +
        # -(log(size_t0/pars['size_opt'])/pars['size_var'])^2 +
        pars['Beta'] * (BA_comp_intra + pars['theta'] * BA_comp_inter) +
        -pars['tau_temp'] * (Temp - pars['optimal_temp'])^2 +
        -pars['tau_prec'] * (Prec - pars['optimal_prec'])^2
      )
    )
  )

	# account for the time interval between sensus
	survival_prob = longev_log^delta_time;
	
  return( survival_prob )
}




# 3. Ingrowth rate as a function of:
# - plot size
# - time interval
# - Basal area of adults from conspecific species
# - Basal area from all adult species

ingrowth_f <- function(
  pars, delta_time, plot_size, BA_adult_sp, BA_adult, Temp, Prec, plot_random
){
  mPlot <- exp(
    pars['mPop_log'] +
    plot_random +
    (-1/pars['sigma_BA']^2) * (BA_adult_sp - pars['optimal_BA'])^2 +
    -pars['tau_temp'] * (Temp - pars['optimal_temp'])^2 +
    -pars['tau_prec'] * (Prec - pars['optimal_prec'])^2
  )

  p <- exp(
    -exp(
      pars['p_log']
    ) +
    BA_adult * -pars['beta_p']
  )

  # mean
  ingrowth <- mPlot * plot_size * (1 - p^delta_time) / (1 - p)

  return( ingrowth )
}
