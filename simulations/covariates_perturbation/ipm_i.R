#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Layer 1 function to run the IPM in the cluster
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# turn off scientific notation
options(scipen = 999)

library(tidyverse)



# Simulation meta parameters (similar for all rows in simulation_pars.RDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replications = 100
perturbation_size = 0.01

# Read parameters and functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Env variables
batch_id <- as.numeric(Sys.getenv('BATCH'))
array_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) + (10000 * batch_id) 

set.seed(array_id)

# simulation variables
sim_pars <- readRDS('simulation_pars.RDS')[array_id, ]

Sp = sim_pars$species_id
Plot = sim_pars$plot_id
Year = sim_pars$year_measured
plot_size = sim_pars$plot_size
temp = sim_pars$bio_01_mean
prec = sim_pars$bio_12_mean

# load parameters
pop_pars <- readRDS('pop_pars.RDS') |>
  filter(species_id == Sp) |>
  select(!species_id) |>
  slice_sample(n = replications)

plot_pars <- readRDS(paste0('plot_parameters/', Sp, '.RDS')) |>
  filter(plot_id == Plot) |>
  filter(year_measured == Year) |>
  pull(plot_pars) %>%
  # jeeez
  .[[1]]

# load IPM funtions
source('R/vital_rates.R')
source('R/kernel.R')
source('R/BasalArea_competition.R')

# Variables range for scaling
vars_rg <- append(
  setNames(
    readRDS('../data/climate_scaleRange.RDS'),
    c('temp', 'prec')
  ),
  list('dbh' = readRDS('dbh_range.RDS'))
)

# functions to scale OR unscale climate variables
scale_vars <- function(
  value, # covariate value
  cov, # either 'temp', 'prec', 'dbh'
  direction, # either 'scale' or 'unscale'
  range_dt = vars_rg
){
  if(cov %in% c('temp', 'prec', 'dbh')){
    cov_rg = range_dt[[cov]]
  }else{
    stop("`cov` must be one of 'temp', 'prec', or 'dbh' character.")
  }

  min_v = cov_rg[1]
  max_v = cov_rg[2]

  if(direction == 'scale') {
    return( (value - min_v)/(max_v - min_v) )
  }else if(direction == 'unscale') {
    return( value * (max_v - min_v) + min_v )
  }else {
    stop('`direction` argument must be either `scale` or `unscale`.')
  }
}

# scale climate variables
temp_scl = scale_vars(temp, 'temp', 'scale')
prec_scl = scale_vars(prec, 'prec', 'scale')

# output table
out_tb <- tibble()

for(i in 1:replications)
{
  # get parameters
  pop_pars_i <- pars_to_list(pop_pars[i, ])
  re_pars_i <- plot_pars[i, ] |> unlist() |> unname()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LAMBDA BASE
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # competition
  N_ref <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0
  )
  
  N_con <- dbh_to_sizeDist(
    dbh = unlist(sim_pars$con_dbh),
    N_intra = N_ref
  )

  # some plots do not have heterospecific individuals
  if(is.null(unlist(sim_pars$het_dbh))) {
    # then use the N_ref object that have N=0 individuals
    N_het <- N_ref
  }else {
    N_het <- dbh_to_sizeDist(
      dbh = unlist(sim_pars$het_dbh),
      N_intra = N_ref
    )
  }

  # get lambda base as reference
  lambda_base <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = temp_scl,
          Prec = prec_scl,
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LAMBDA BA_con
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # scale and perturb individual dbh
  dbh_con_scl <- scale_vars(unlist(sim_pars$con_dbh), 'dbh', 'scale')
  dbh_con_scl_pertb <- dbh_con_scl + perturbation_size
  dbh_con_pertb <- scale_vars(dbh_con_scl_pertb, 'dbh', 'unscale')
  
  N_con_pertb <- dbh_to_sizeDist(
    dbh = dbh_con_pertb,
    N_intra = N_ref
  )

  # get lambda base as reference
  lambda_bacon <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con_pertb,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = temp_scl,
          Prec = prec_scl,
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LAMBDA BA_het
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  # some plots do not have heterospecific individuals
  if(is.null(unlist(sim_pars$het_dbh))) {
    # if so, just skip this
    N_het_pertb = N_het
    lambda_bahet = lambda_base
  }else {
    # scale and perturb individual dbh
    dbh_het_scl <- scale_vars(unlist(sim_pars$het_dbh), 'dbh', 'scale')
    dbh_het_scl_pertb <- dbh_het_scl + perturbation_size
    dbh_het_pertb <- scale_vars(dbh_het_scl_pertb, 'dbh', 'unscale')
    
    N_het_pertb <- dbh_to_sizeDist(
      dbh = dbh_het_pertb,
      N_intra = N_ref
    )

    # get lambda base as reference
    lambda_bahet <- max(
      Re(
        eigen(
          mkKernel(
            Nvec_intra = N_con,
            Nvec_inter = N_het_pertb,
            delta_time = 1,
            plotSize = plot_size,
            Temp = temp_scl,
            Prec = prec_scl,
            pars = pop_pars_i,
            plot_random = re_pars_i
          )$K
        )$values
      )
    )
  }


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LAMBDA TEMP
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # get lambda base as reference
  lambda_temp <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = temp_scl + perturbation_size,
          Prec = prec_scl,
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LAMBDA PREC
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # get lambda base as reference
  lambda_prec <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = temp_scl,
          Prec = prec_scl + perturbation_size,
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # prepare output
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  # assign to out table
  # `par.` refers to partial derivative
  out_tb <- rbind(
    out_tb,
    tibble(
      rep = i,
      temp = temp,
      prec = prec,
      BA_con = size_to_BAplot(N = N_con, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda_base = lambda_base,
      # BA_con perturbation
      BA_con_pertb = size_to_BAplot(N = N_con_pertb, plot_size),
      par.BA_con = abs(lambda_bacon - lambda_base)/(size_to_BAplot(N = N_con_pertb, plot_size) - size_to_BAplot(N = N_con, plot_size)),
      # BA_het perturbation
      BA_het_pertb = size_to_BAplot(N = N_het_pertb, plot_size),
      par.BA_het = abs(lambda_bahet - lambda_base)/(size_to_BAplot(N = N_het_pertb, plot_size) - size_to_BAplot(N = N_het, plot_size)),
      # temp perturbation
      temp_pertb = scale_vars(temp_scl + perturbation_size, 'temp', 'unscale'),
      par.temp = abs(lambda_temp - lambda_base)/(perturbation_size),
      # prec perturbation
      prec_pertb = scale_vars(prec_scl + perturbation_size, 'prec', 'unscale'),
      par.prec = abs(lambda_prec - lambda_base)/(perturbation_size)
    )
  )
}



# Save output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveRDS(
  out_tb,
  paste0('output/sim_', array_id, '.RDS')
)
