library(tidyverse)
library(ggtext)

# load simulation info parameters
pars_sim <- read_csv('simulations/coexistence/output/simulation_pars.csv')

# function to retrieve info from each simulation
get_sim <- function(array_id) {
  # read sim
  sim_i <- readRDS(
    paste0('simulations/coexistence/output/sim_', array_id, '.RDS')
  )

  # species i alone
  sp1_subsim1 <- c(
    head(sim_i[['year_summ_subsim1']], 1)$lambda,
    tail(sim_i[['year_summ_subsim1']], 1)$N_sp1,
    tail(sim_i[['year_summ_subsim1']], 1)$BA_sp1
  )

  # species j alone
  sp2_subsim1 <- c(
    head(sim_i[['year_summ_sp2_subsim2']], 1)$lambda,
    tail(sim_i[['year_summ_sp2_subsim2']], 1)$N_sp2,
    tail(sim_i[['year_summ_sp2_subsim2']], 1)$BA_sp2
  )

  # species i and j together
  sp1_subsim2 <- c(
    head(sim_i[['year_summ_sp1sp2_subsim2']], 1)$lambda_sp1,
    tail(sim_i[['year_summ_sp1sp2_subsim2']], 1)$N_sp1,
    tail(sim_i[['year_summ_sp1sp2_subsim2']], 1)$BA_sp1
  )
  sp2_subsim2 <- c(
    head(sim_i[['year_summ_sp1sp2_subsim2']], 1)$lambda_sp2,
    tail(sim_i[['year_summ_sp1sp2_subsim2']], 1)$N_sp2,
    tail(sim_i[['year_summ_sp1sp2_subsim2']], 1)$BA_sp2
  )

  as.data.frame(
    t(
      as.matrix(
        tibble(
          sp1_subsim1, sp2_subsim1, sp1_subsim2, sp2_subsim2
        )
      )
    )
  ) |>
  rownames_to_column() |>
  rename(
    lambda = V1,
    N_eq = V2,
    BA_eq = V3
  ) |>
  mutate(
    sp = gsub('_.*', '', rowname),
    subsim = gsub('.*_', '', rowname)
  ) |>
  select(!rowname)
}

out <- parse_number(dir('simulations/coexistence/output/', pattern = 'RDS')) |>
  sort() |>
  map_dfr(
    ~ get_sim(.x) |>
    bind_cols(sim = .x)
  )


# compute sensitivity
sens <- out |>
  filter(sp == 'sp1') |>
  group_by(sim) |>
  mutate(
    sens_lambda = (lambda[subsim == 'subsim1'] - lambda[subsim == 'subsim2'])/lambda[subsim == 'subsim1'],
    sens_N = (N_eq[subsim == 'subsim1'] - N_eq[subsim == 'subsim2'])/N_eq[subsim == 'subsim1'],
    sens_BA = (BA_eq[subsim == 'subsim1'] - BA_eq[subsim == 'subsim2'])/BA_eq[subsim == 'subsim1']
  ) |>
  slice_head(n = 1) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, seed)
      )
  ) |>
  as_tibble()


# Matrix of sensitivity for pair of species
sens |>
  group_by(sp_pair) |>
  mutate(
    mean_lambda = mean(sens_lambda),
    mean_N = mean(sens_N),
    mean_BA = mean(sens_BA)
  ) |>
  slice_head(n = 1) |>
  ungroup() |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(species_name_sp1 = species_name) |>
      select(species_id_old, species_name_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(species_name_sp2 = species_name) |>
      select(species_id_old, species_name_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  ggplot(aes(species_name_sp2, species_name_sp1)) + 
    geom_raster(aes(fill = mean_BA)) + 
    scale_fill_gradient2() +
    theme_classic() +
    labs(
      fill = 'Sensitivity',
      title = 'Mean basal area sensitivity of invader species *i* to the resident species *j* (*S<sub>ij</sub>*)'
    ) +
    theme(
      axis.text.x = element_text(
        size = 9, face = 'italic', angle = 45, hjust = 1, vjust = 1
      ),
      axis.text.y = element_text(size = 9, face = 'italic'),
      plot.title = element_markdown()
    ) +
    xlab('Resident species') +
    ylab('Invader species')

