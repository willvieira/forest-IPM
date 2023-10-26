# script to merge output simulations into a single file
# Will Vieira
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)

# load sim and metasim files
sims <- readRDS('simulation_pars.RDS')
sim_clim <- sims[['climate']]
sim_comp <- sims[['competition']]
metasim <- read_csv('simulation_metapars.csv')


# load parameters into data.frame
1:nrow(metasim) |>
  map(~ readLines(paste0('output/sim_', .x, '.csv'))) |>
  map(~ str_split(.x, pattern = ',', simplify = TRUE)) |>
  map(as.numeric) ->
vals_ls

vals_df <- as_tibble(do.call(rbind, vals_ls))
names(vals_df) = str_split(readLines('output/sim_0.csv'), ',', simplify = TRUE)


# merge with sim info
sim_clim |>
  mutate(sim = 'climate') |>
  bind_rows(
    sim_comp |>
      mutate(
        sim = 'competition',
        clim = NA,
        comp = as.character(comp)
      )
  ) |>
  mutate(array_id = 1:n()) |>
  left_join(
    vals_df |>
      mutate(array_id = 1:n())
  ) |>
  saveRDS('output_complete.RDS')
