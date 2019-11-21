##############################################################
# Vital rate functions for growth, mortality and reproduction
# Will Vieira
# January 20, 2019
##############################################################



######################
 # 1. Growth
 # 2. Mortality
 # 3. Reproduction
 ######################



# 2. growth function in function of size (x), temperature and precipitaion (E[1 or 2]) and competition (C)
 growth_XEC <- function(x, E, C, parsGrowth)
 {

   # annual growth rate
   growth = parsGrowth$pdg *
            (C + (1 - C) * parsGrowth$beta) * # competition effect
            exp(-(E[1] - parsGrowth$A)/parsGrowth$B) * exp((-1 - parsGrowth$C) * log(1 + exp(-(E[1] - parsGrowth$A)/parsGrowth$B))) * parsGrowth$C *(1 + 1/parsGrowth$C)^(1 + parsGrowth$C) * # temp effect
            exp(-(E[2] - parsGrowth$P_opt) * (E[2] - parsGrowth$P_opt)/parsGrowth$sigmaP_opt^2) * # prec effet
            exp(-log(x/parsGrowth$Phi_opt)*log(x/parsGrowth$Phi_opt)/parsGrowth$sigmaPhi_opt^2) # size effect

   return(growth)
 }
#



# 2. Probability of mortality (Surv = 1 - mort) in function of size (x), temperature and precipitation (E[1 or 2]) and competition (C)
  mort_XEC = function(x, E, C, parsMort)
  {

    lambda = log(99)/(parsMort$DBH001 * (1 - parsMort$theta))

    firstPart = parsMort$psi * # optimal survival
                (C + (1 - C) * parsMort$beta) * # competition effect
                exp(-0.5 * (E[1] - parsMort$T_opt)^2/parsMort$sigmaT_opt^2) * # temp effect
                exp(-0.5 * (E[2] - parsMort$P_opt)^2/parsMort$sigmaP_opt^2) * # prec effect
                (x/100)^parsMort$phi/(1 + exp(lambda * ((x/10) - parsMort$theta * parsMort$DBH001))) # size effect

    return(1/(1 + firstPart)) # logit
  }
#



# 3. Fecundity function (#TODO to be defined)

  gLogit <- function(x, L, U, B, C, v, Q, M)
  {
    # https://en.wikipedia.org/wiki/Generalised_logistic_function
    return( L + ((U - L)/(C + exp(-B * (x - M)))^(1/v)) )
  }

  # this is a constant times a normal distribution times an exponential function
  fec_X = function(y, x, parsFec)
  {

    # establishment proportion
    parsFec$establishment.prob *
    # reproduction probability
    gLogit(x,
      L = 0, U = 1, B = 0.05, C = 1, v = 1, Q = 1, M = parsFec$location) *
    # Recruitment size distribition
    dnorm(y,
      mean = parsFec$recruit.size.mean,
      sd = parsFec$recruit.size.sd) *
    # Seed production
    exp(parsFec$seed.int + parsFec$seed.slope * x)
  }
#
