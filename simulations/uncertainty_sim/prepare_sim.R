# script to set climate, competition, and initial conditions for sensitivity analysis
# Will Vieira
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# /!\ This should be ran locally (pre server run) /!\

library(tidyverse)

set.seed(0.0)


# Competition and climate range
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

treeData <- readRDS(
  paste0(readLines('_data.path'), 'treeData.RDS')
)

spIds <- read_csv(
  paste0(readLines('_data.path'), 'species_id.csv')
)

sp_analysis <- c('NAQUEPRI', '183302PICMAR', '183319PINBAN', '19447QUEVEL', '32931FRAAME')

temp_mean  <- c(7.8, -3.8, 7.7, 17.4, 19.0)
temp_sd <- seq(0.1, 1.5, 0.2)
BA_comp <- c(51, 41.8, 31.5, 44.3, 41.1)
BA_sd <- seq(1, 8, 1)



# using the 500 draws from the sensAnalyis_v3 to compute demo stochasticity
readRDS('simulations/sensAnalysis_v3/output_processed/lambdas.RDS') |>
  filter(species_id %in% sp_analysis) |>
  filter(comp == 'high' & clim == 'cold') |>
  select(!c(comp, clim)) ->
pars

# define temperature variability
tibble(
    species_id = sp_analysis,
    temp_mean = temp_mean
  ) |>
  group_by(species_id) |>
  expand_grid(
    temp_sd
  ) |>
  group_by(species_id, temp_sd) |>
  expand_grid(
    replication = 1:500
  ) |>
  mutate(
    temp = rnorm(n(), temp_mean, temp_sd)
  ) ->
temp_sim

# define competition variability
tibble(
    species_id = sp_analysis,
    BA_mean = BA_comp
  ) |>
  group_by(species_id) |>
  expand_grid(
    BA_sd
  ) |>
  group_by(species_id, BA_sd) |>
  expand_grid(
    replication = 1:500
  ) |>
  mutate(
    BA = rnorm(n(), BA_mean, BA_sd)
  ) ->
BA_sim


# Initial condition and competition conditions for each species
# So each replication uses the exact the same conditions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get lower `Lmax` value per species
spIds$species_id_old |>
  set_names() |>
  map(
    ~readRDS(
      paste0(readLines('_data.path'), 'output_sim_processed/growth/intcpt_plot_comp_clim/posterior_pop_', .x, '.RDS')
    ) |>
    group_by(par) |>
    reframe(min_value = min(value)) |>
    filter(par == 'Lmax') |>
    pull(min_value)
  ) |>
  as_vector() ->
Lmax_sp


# generate N_sp and N_competition for each species
source('R/kernel.R')
source('R/BasalArea_competition.R')

N_dist_folder <- file.path('simulations', 'uncertainty_sim', 'N_dist')
dir.create(N_dist_folder)

new_BA_sim = list()
for(Sp in sp_analysis)
{
  # N_init with small pop size
  pars_fake <- list('growth' = c('Lmax' = unname(Lmax_sp[Sp])))

  # Initial state to be used for all replications
  N_init <- init_pop(
    params = pars_fake,
    L = 127,
    h = 1,
    N = 0.1
  )

  # Competitive state to be used in fixed competition state
  # expected 
  BA_sim |>
    filter(species_id == Sp) ->
  BA_sim_sp
  
  BA_ls <- list()
  for(i  in 1:nrow(BA_sim_sp))
  {
    BA_expected  = BA_sim_sp$BA[i]
    
    fact = 3
    N_comp <- init_pop(
      params = pars_fake,
      L = 127,
      h = 1,
      N = BA_expected/fact
    )
    
    nb_try = 0

    diff_BA = BA_expected - size_to_BAplot(N_comp, plot_size = 180)

    while(abs(diff_BA) > 0.1 & nb_try <= 10)
    {
      if(abs(diff_BA) > 0.1) {
        if(diff_BA < 0) {
          fact = fact + 0.02
        }else{
          fact = fact - 0.02
        }
      }else{
        nb_try = nb_try + 1
      }
      N_comp <- init_pop(
        params = pars_fake,
        L = 127,
        h = 1,
        N = BA_expected/fact
      )

      BA_estimated = size_to_BAplot(N_comp, plot_size = 180)

      diff_BA = BA_expected - BA_estimated
    }

    BA_sim_sp$BA[i] = BA_estimated
    BA_ls[[paste0(BA_sim_sp$BA_sd[i], '_', BA_sim_sp$replication[i])]] = N_comp

    cat(round(i/nrow(BA_sim_sp) * 100, 1), '%\r')
  }

  saveRDS(
    list(
      'N_init' = N_init,
      'N_comp' = BA_ls
    ),
    paste0(N_dist_folder, '/N_', Sp, '.RDS')
  )

  new_BA_sim[[Sp]] <- BA_sim_sp
}


# Scale temperature variables
temp_rg <- readRDS(file.path(readLines('_data.path'), 'climate_scaleRange.RDS'))$bio_01_mean

# Save sim and metasim file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Parameters only
pars |>
  left_join(
    tibble(
      species_id = sp_analysis,
      temp = (temp_mean - temp_rg[1])/(temp_rg[2] - temp_rg[1])
    )
  ) |>
  mutate(arrayID = 1:n()) |>
  select(!lambda) ->
pars_sim

# temperature + parameters
temp_sim |>
  left_join(
    pars |>
      select(!lambda)
  ) |>
  # scale temp
  mutate(
    temp = (temp - temp_rg[1])/(temp_rg[2] - temp_rg[1])
  ) |>
  mutate(arrayID = (1:n()) + nrow(pars_sim)) ->
temp_sim

# BA + parameters
new_BA_sim |>
  bind_rows() |>
  left_join(
    pars |>
      select(!lambda)
  ) |>
  mutate(arrayID = (1:n()) + (nrow(pars_sim) + nrow(temp_sim))) ->
BA_sim

# BA + temperature + parameters
new_BA_sim |>
  bind_rows() |>
  # this is necessary to merge because it cannot identify BA_sp and temp_sd due to diff values
  mutate(arrayID = (1:n()) + nrow(pars_sim)) |>
  left_join(
    temp_sim,
    by = c('species_id', 'replication', 'arrayID')
  ) |>
  mutate(arrayID = (1:n()) + (nrow(pars_sim) + nrow(temp_sim) + BA_sim)) ->
BAtemp_sim



tibble(
    sim = c(rep('pars', nrow(pars_sim)), rep('temp', nrow(temp_sim)), rep('BA', nrow(BA_sim)), rep('BAtemp', nrow(BAtemp_sim)))
  ) |>
  mutate(
  deltaTime = 1,
  plotSize = 180,
  param_method = 'random'
  ) |>
  write_csv(file.path('simulations', 'uncertainty_sim', 'simulation_metapars.csv'))

saveRDS(
  list(
    'pars' = pars_sim,
    'temp' = temp_sim,
    'BA' = BA_sim,
    'BAtemp' = BAtemp_sim
  ),
  file.path('simulations', 'uncertainty_sim', 'simulation_pars.RDS')
)
