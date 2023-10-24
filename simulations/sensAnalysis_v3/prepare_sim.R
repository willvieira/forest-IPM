# script to set climate, competition, and initial conditions for sensitivity analysis
# Will Vieira
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# /!\ This should be ran locally (pre server run) /!\

library(tidyverse)

set.seed(0.0)


# Competition and climate range
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load treeData to set climate and competition limits
# the threshold for climate and competition will be set to the 10% quantile dist
# range hot: lat <- quantile(Mean Annual Temperature, 1%)
# range cold: lat <- quantile(Mean Annual Temperature, 99%)
# range center: lat <- quantile(Mean Annual Temperature, 50%)
# comp low: BA_comp_sp = 0.5
# comp high: BA_comp_sp <- quantile(BA_comp_sp, 99%)

treeData <- readRDS(
  paste0(readLines('_data.path'), 'treeData.RDS')
)

spIds <- read_csv(
  paste0(readLines('_data.path'), 'species_id.csv')
)

treeData |>
  filter(species_id %in% spIds$species_id_old) |>
  group_by(species_id) |>
  reframe(
    temp_hot = quantile(bio_01_mean_scl, 0.01, na.rm = TRUE),
    temp_center = quantile(bio_01_mean_scl, 0.5, na.rm = TRUE),
    temp_cold = quantile(bio_01_mean_scl, 0.99, na.rm = TRUE),
    comp_high = quantile(BA_sp, 0.99, na.rm = TRUE)
  ) |>
  write_csv(
    file = file.path('simulations', 'sensAnalysis_v3', 'comp_clim_range.csv')
  )




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

N_dist_folder <- file.path('simulations', 'sensAnalysis_v3', 'N_dist')
dir.create(N_dist_folder)

comp_rg <- read_csv(file.path('simulations', 'sensAnalysis_v3', 'comp_clim_range.csv'))

# save estimated BA_plot to compare with the expected from the quantile dist
# This is because the init function requires a total population size, and not
# plot basal area. To get the expected plot basal area from the population size
# there is a litle algorithm to brute force the best N for each species to get
# the closest expected plot basal area
comp_rg$comp_estimated = 0
comp_rg$comp_init = 0
for(Sp in spIds$species_id_old)
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
  comp_rg$comp_init[comp_rg$species_id == Sp] = size_to_BAplot(N_init, plot_size = 180)

  # expected 
  BA_expected <- comp_rg |> filter(species_id == Sp) |> pull(comp_high)
  
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

  comp_rg$comp_estimated[comp_rg$species_id == Sp] = BA_estimated

  saveRDS(
    list(
      'N_init' = N_init,
      'N_comp' = N_comp
    ),
    paste0(N_dist_folder, '/N_', Sp, '.RDS')
  )
}

  # comp_rg |>
  #   ggplot() + 
  #   aes(comp_high, comp_estimated) +
  #   geom_point() +
  #   geom_abline(slope = 1, intercept = 0) +
  #   theme_classic()
