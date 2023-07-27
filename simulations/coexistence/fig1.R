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


# Sensitivity distribution of species i across all other species
p1 <- sens |>
  select(sp1, sp2, sens_BA) |>
  filter(
    sens_BA > quantile(sens_BA, probs = 0.001) &
    sens_BA < quantile(sens_BA, probs = 0.999)
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(!species_id_new),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp2 = shade) |>
      select(species_id_old, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    species_name = fct_reorder(species_name, sens_BA),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(sens_BA, species_name, fill = shade_sp2)) +
    ggridges::geom_density_ridges2(
      color = 'transparent',
      alpha = 0.7
    ) +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab('Basal area sensitivity') +
    theme(
      axis.text.y = element_text(face = "italic"),
      legend.position = 'top'
    ) +
    labs(
      fill = 'Shade tolerance\nof resident species'
    )

p2 <- sens |>
  select(sp1, sens_BA) |>
  filter(
    sens_BA > quantile(sens_BA, probs = 0.001) &
    sens_BA < quantile(sens_BA, probs = 0.999)
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(!species_id_new),
    by = c('sp1' = 'species_id_old')
  ) |>
  mutate(
    shade = factor(shade, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade, sens_BA)) +
    geom_boxplot(alpha = 0.8, fill = c("#87bc45", "#edbf33", "#ea5545")) +
    theme_minimal() +
    xlab('Shade tolerance of invader species') +
    ylab('Basal area sensitivity')

ggpubr::ggarrange(p1, p2, ncol = 2)



# correlation between BA sensitivity and a_ij
out |>
  filter(sp == 'sp1') |>
  select(sim, subsim, BA_eq) |>
  pivot_wider(
    names_from = subsim,
    values_from = BA_eq
  ) |>
  mutate(
    # first deal with 0 BA or values too close to it
    subsim2 = ifelse(subsim2 < 0.1, 0.1, subsim2),
    subsim1 = ifelse(subsim1 < 0.1, 0.1, subsim1),
    a_ii = 1/subsim1,
    a_ij = log(subsim2/subsim1),
    sens = (subsim1 - subsim2)/subsim1
  ) |>
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
  ggplot(aes(sens, a_ij)) +
    geom_point() +
    theme_minimal() +
    xlab('BA sensitivity') +
    ylab(expression(alpha[ij])) +
    lab()


## Compute alpha_ii with the equation 1/BA*
p1 <-out |>
  filter(sp == 'sp1' & subsim == 'subsim1') |>
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
  as_tibble() |>
  mutate(
    K = 1/BA_eq
  ) |>
  filter(K < quantile(K, probs = 0.925)) |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(species_id_old, species_name),
    by = c('sp1' = 'species_id_old')
  ) |>
  ggplot(aes(K, fct_reorder(species_name, K))) +
    ggridges::geom_density_ridges2(color = rgb(0.3,0.5,0.4,0.6), fill = rgb(0.3,0.5,0.4,0.6), alpha = 0.8) +
    theme_classic() +
    ylab('Invader species') +
    xlab(expression(alpha[ii])) +
    theme(
      axis.text.y = element_text(face = "italic")
    )

p2 <- out |>
  filter(sp == 'sp1' & subsim == 'subsim1') |>
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
  as_tibble() |>
  mutate(
    K = 1/BA_eq
  ) |>
  filter(K < quantile(K, probs = 0.9)) |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(species_id_old, shade),
    by = c('sp1' = 'species_id_old')
  ) |>
  mutate(
    shade = factor(shade, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade, K)) +
    geom_boxplot(
      alpha = 0.8,
      fill = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_minimal() +
    xlab('Shade tolerance') +
    ylab(expression(alpha[ii]))

ggpubr::ggarrange(p1, p2, ncol = 2)

