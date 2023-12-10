#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Layer 1 function to run the IPM in the cluster
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# turn off scientific notation
options(scipen = 999)

library(tidyverse)



# Simulation meta parameters (similar for all rows in simulation_pars.RDS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
replications = 100



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

# climate scaling
clim_rg <- readRDS('../data/climate_scaleRange.RDS')
scale_clim <- function(value, var, climRange = clim_rg) {
  if(var == 'temp') {
    return((value - climRange[[1]][1])/(climRange[[1]][2] - climRange[[1]][1]))
  }else if(var == 'prec') {
    return((value - climRange[[2]][1])/(climRange[[2]][2] - climRange[[2]][1]))
  }else{
    stop('`var` must be either `temp` or `prec` character.')
  }
}


# output table
out_tb <- tibble()

# Below there will be four set of simulations
# Sim0: par mean, env mean, space mean
# Sim1: par var,  env mean, space mean
# Sim2: par var,  env var,  space mean
# Sim3: par var,  env var,  space var

# For each simulation set, we will compute lambda for three conditions
# cond1: lambda without competition
# cond2: lambda with heterospecific competition only
# cond3: lambda with heterospecific and conspecific competition



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sim 0: par mean, condition man, space mean

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pop_pars_i <- pop_pars |>
  reframe(across(everything(), mean)) |>
  pars_to_list()
re_pars_i <- rep(0, 3)

# competition
N_ref <- init_pop(
  pop_pars_i,
  L = 127,
  h = 1,
  N = 0.1
)

N_comp0 <- init_pop(
  pop_pars_i,
  L = 127,
  h = 1,
  N = 0
)


# CONDITION 1:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get lambda
lbd <- max(
  Re(
    eigen(
      mkKernel(
        Nvec_intra = N_ref,
        Nvec_inter = N_comp0,
        delta_time = 1,
        plotSize = plot_size,
        Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
        Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
        pars = pop_pars_i,
        plot_random = re_pars_i
      )$K
    )$values
  )
)

# assign to out table
out_tb <- rbind(
  out_tb,
  tibble(
    sim = 0,
    cond = 1,
    rep = 0,
    temp = sim_pars$bio_01_mean,
    temp_sd = 0,
    prec = sim_pars$bio_12_mean,
    prec_sd = 0,
    BA_con = size_to_BAplot(N = N_ref, plot_size),
    BA_het = 0,
    lambda = lbd
  )
)

# CONDITION 2:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# some plots does not have heterospecific individuals
if(is.null(unlist(sim_pars$het_dbh))) {
  # then use the N_ref object that have N=0 individuals
  N_het <- N_comp0
}else {
  N_het <- dbh_to_sizeDist(
    dbh = unlist(sim_pars$het_dbh),
    N_intra = N_ref
  )
}

# get lambda
lbd <- max(
  Re(
    eigen(
      mkKernel(
        Nvec_intra = N_ref,
        Nvec_inter = N_het,
        delta_time = 1,
        plotSize = plot_size,
        Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
        Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
        pars = pop_pars_i,
        plot_random = re_pars_i
      )$K
    )$values
  )
)

# assign to out table
out_tb <- rbind(
  out_tb,
  tibble(
    sim = 0,
    cond = 2,
    rep = 0,
    temp = sim_pars$bio_01_mean,
    temp_sd = 0,
    prec = sim_pars$bio_12_mean,
    prec_sd = 0,
    BA_con = size_to_BAplot(N = N_ref, plot_size),
    BA_het = size_to_BAplot(N = N_het, plot_size),
    lambda = lbd
  )
)

# CONDITION 3:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N_con <- dbh_to_sizeDist(
  dbh = unlist(sim_pars$con_dbh),
  N_intra = N_ref
)

# get lambda
lbd <- max(
  Re(
    eigen(
      mkKernel(
        Nvec_intra = N_con,
        Nvec_inter = N_het,
        delta_time = 1,
        plotSize = plot_size,
        Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
        Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
        pars = pop_pars_i,
        plot_random = re_pars_i
      )$K
    )$values
  )
)

# assign to out table
out_tb <- rbind(
  out_tb,
  tibble(
    sim = 0,
    cond = 3,
    rep = 0,
    temp = sim_pars$bio_01_mean,
    temp_sd = 0,
    prec = sim_pars$bio_12_mean,
    prec_sd = 0,
    BA_con = size_to_BAplot(N = N_con, plot_size),
    BA_het = size_to_BAplot(N = N_het, plot_size),
    lambda = lbd
  )
)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sim 1: par var, condition man, space mean

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in 1:replications)
{
  # get parameters
  pop_pars_i <- pars_to_list(pop_pars[i, ])
  re_pars_i <- rep(0, 3)
  
  # competition
  N_ref <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0.1
  )
  
  N_comp0 <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0
  )


  # CONDITION 1:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_comp0,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
          Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 1,
      cond = 1,
      rep = i,
      temp = sim_pars$bio_01_mean,
      temp_sd = 0,
      prec = sim_pars$bio_12_mean,
      prec_sd = 0,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = 0,
      lambda = lbd
    )
  )


  # CONDITION 2:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # some plots does not have heterospecific individuals
  if(is.null(unlist(sim_pars$het_dbh))) {
    # then use the N_ref object that have N=0 individuals
    N_het <- N_comp0
  }else {
    N_het <- dbh_to_sizeDist(
      dbh = unlist(sim_pars$het_dbh),
      N_intra = N_ref
    )
  }

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
          Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 1,
      cond = 2,
      rep = i,
      temp = sim_pars$bio_01_mean,
      temp_sd = 0,
      prec = sim_pars$bio_12_mean,
      prec_sd = 0,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )


  # CONDITION 3:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  N_con <- dbh_to_sizeDist(
    dbh = unlist(sim_pars$con_dbh),
    N_intra = N_ref
  )

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(sim_pars$bio_01_mean, var = 'temp'),
          Prec = scale_clim(sim_pars$bio_12_mean, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 1,
      cond = 3,
      rep = i,
      temp = sim_pars$bio_01_mean,
      temp_sd = 0,
      prec = sim_pars$bio_12_mean,
      prec_sd = 0,
      BA_con = size_to_BAplot(N = N_con, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )

}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sim 2: par var, condition var, space mean

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in 1:replications)
{
  # get parameters
  pop_pars_i <- pars_to_list(pop_pars[i, ])
  re_pars_i <- rep(0, 3)
  
  # competition
  N_ref <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0.1
  )
  
  N_comp0 <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0
  )

  # generate random climate
  temp_random <- rnorm(1, sim_pars$bio_01_mean, sim_pars$bio_01_sd)
  prec_random <- rnorm(1, sim_pars$bio_12_mean, sim_pars$bio_12_sd)


  # CONDITION 1:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_comp0,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 2,
      cond = 1,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = 0,
      lambda = lbd
    )
  )


  # CONDITION 2:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # some plots does not have heterospecific individuals
  if(is.null(unlist(sim_pars$het_dbh))) {
    # then use the N_ref object that have N=0 individuals
    N_het <- N_comp0
  }else {
    N_het <- dbh_to_sizeDist(
      dbh = unlist(sim_pars$het_dbh),
      N_intra = N_ref
    )
  }

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 2,
      cond = 2,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )


  # CONDITION 3:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  N_con <- dbh_to_sizeDist(
    dbh = unlist(sim_pars$con_dbh),
    N_intra = N_ref
  )

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 2,
      cond = 3,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_con, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Sim 3: par var, condition var, space var

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(i in 1:replications)
{
  # get parameters
  pop_pars_i <- pars_to_list(pop_pars[i, ])
  re_pars_i <- plot_pars[i, ] |> unlist() |> unname()
  
  # competition
  N_ref <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0.1
  )
  
  N_comp0 <- init_pop(
    pop_pars_i,
    L = 127,
    h = 1,
    N = 0
  )

  # generate random climate
  temp_random <- rnorm(1, sim_pars$bio_01_mean, sim_pars$bio_01_sd)
  prec_random <- rnorm(1, sim_pars$bio_12_mean, sim_pars$bio_12_sd)


  # CONDITION 1:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_comp0,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 3,
      cond = 1,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = 0,
      lambda = lbd
    )
  )


  # CONDITION 2:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # some plots does not have heterospecific individuals
  if(is.null(unlist(sim_pars$het_dbh))) {
    # then use the N_ref object that have N=0 individuals
    N_het <- N_comp0
  }else {
    N_het <- dbh_to_sizeDist(
      dbh = unlist(sim_pars$het_dbh),
      N_intra = N_ref
    )
  }

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_ref,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 3,
      cond = 2,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_ref, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )


  # CONDITION 3:
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  N_con <- dbh_to_sizeDist(
    dbh = unlist(sim_pars$con_dbh),
    N_intra = N_ref
  )

  # get lambda
  lbd <- max(
    Re(
      eigen(
        mkKernel(
          Nvec_intra = N_con,
          Nvec_inter = N_het,
          delta_time = 1,
          plotSize = plot_size,
          Temp = scale_clim(temp_random, var = 'temp'),
          Prec = scale_clim(prec_random, var = 'prec'),
          pars = pop_pars_i,
          plot_random = re_pars_i
        )$K
      )$values
    )
  )

  # assign to out table
  out_tb <- rbind(
    out_tb,
    tibble(
      sim = 3,
      cond = 3,
      rep = i,
      temp = temp_random,
      temp_sd = sim_pars$bio_01_sd,
      prec = prec_random,
      prec_sd = sim_pars$bio_12_sd,
      BA_con = size_to_BAplot(N = N_con, plot_size),
      BA_het = size_to_BAplot(N = N_het, plot_size),
      lambda = lbd
    )
  )
}



# Save output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

saveRDS(
  out_tb,
  paste0('output/sim_', array_id, '.RDS')
)
