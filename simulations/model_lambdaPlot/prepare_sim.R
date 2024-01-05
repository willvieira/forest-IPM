# script to set data parameters to run linear model - lambda against temperature
# Will Vieira
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# /!\ This should be ran locally (pre server run) /!\

library(tidyverse)

set.seed(0.0)

sim_name <- 'simulations/model_lambdaPlot'


# For each species-border condition, there are 6 different simulations
# This yields a total of 372 stan models to be run
pars <- readRDS('simulations/lambda_plot/simulation_pars.RDS')

pars |>
  ungroup() |>
  select(species_id) |>
  distinct() |>
  group_by(species_id) |>
  expand_grid(
    border = c('cold', 'hot'),
    sim = 1:2,
    cond = 1:3
  ) |>
  mutate(array_id = row_number()) |>
  saveRDS(file = paste0(sim_name, '/simulation_pars.RDS'))
