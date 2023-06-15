# N equilibrium figure

library(tidyverse)
library(ggdist)

# load parameter simulations
sim_pars <- read_csv('simulations/N_equilibrium/output/simulation_pars.csv')

out <- parse_number(
    dir(
      'simulations/N_equilibrium/output/', pattern = '.RDS'
    )
  ) |>
  sort() |>
  map_dfr(
    ~readRDS(
      paste0('simulations/N_equilibrium/output/sim_', .x, '.RDS')
    )[['year_summ']] |>
    bind_cols(sim = .x)
  ) |>
  # add parameter simulations
  left_join(
    sim_pars |>
      select(!c(param_method, seed, contains('random'))) |>
      mutate(sim = row_number())
  )

# distribution of N equilibrium across species
out |>
  group_by(species_id) |>
  filter(year == max(year)) |>
  ggplot(aes(N, species_id)) +
    stat_pointinterval() +
    theme_minimal() +
    xlab('Population size') +
    ylab('')

# is population distribution at equilibirum stable?
out |>
  group_by(species_id) |>
  mutate(q80 = quantile(year, prob = 0.8)) |>
  filter(year > q80) |>
  ggplot(aes(year, N, color = as.factor(replication))) +
    geom_line(alpha = 0.5) +
    facet_wrap(~species_id, scales = 'free') +
    theme_minimal() +
    theme(legend.position = 'none')


# # special case for 18086LIRTUL and 18048JUNVIR
# out |>
#   filter(species_id == '18086LIRTUL' & year %in% 1500:2000) |>
#   ggplot(aes(year, N, color = as.factor(replication))) +
#     geom_line() +
#     theme(legend.position = 'none')

# out |>
#   filter(species_id == '18048JUNVIR' & replication %in% c(23, 47)) |>
#   ggplot(aes(year, N, color = as.factor(replication))) +
#     geom_line()

# library(plotly)
# x = readRDS('simulations//N_equilibrium/output/sim_123.RDS')[[1]]
# ggplotly(
#   x[, 1300:2000] |>
#     reshape2::melt() |>
#     ggplot(aes(Var1, value)) +
#       geom_line(aes(frame = Var2))
# )
