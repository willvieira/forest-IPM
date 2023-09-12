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



# Compute coexistence metrics from the sensisitivy
## Niche difference (ND)
## Relative fitness difference (RFD)

coex <- sens |>
  select(replication, sp1, sp2, contains('sens_')) |>
  pivot_longer(
    cols = contains('sens_'),
    values_to = 'sens_sp1'
  ) |>
  mutate(
    sens_sp1 = ifelse(sens_sp1 < 0.001, 0.001, sens_sp1),
    name = gsub('sens_', '', name)
  ) |>
  full_join(
   sens |>
    select(replication, sp1, sp2, contains('sens_')) |>
    pivot_longer(
      cols = contains('sens_'),
      values_to = 'sens_sp2'
    ) |>
    mutate(
      sens_sp2 = ifelse(sens_sp2 < 0.001, 0.001, sens_sp2),
      name = gsub('sens_', '', name)
    ),
    by = c(
      'replication' = 'replication',
      'sp2' = 'sp1', 'sp1' = 'sp2',
      'name' = 'name'
    )
  ) |>
  drop_na() |>
  mutate(
    ND = 1 - sqrt(sens_sp1 * sens_sp2),
    RFD = sqrt(sens_sp1/sens_sp2)
  )
  
  
coex |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(species_id_old, shade),
    by = c('sp2' = 'species_id_old')
  ) |>
  ggplot(aes(ND, RFD, color = shade)) +
    geom_point(alpha = 0.5, size = 0.5) +
    facet_grid(~name, scales = 'free')


coex |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp1 = shade,
        species_name_sp1 = species_name
      ) |>
      select(species_id_old, species_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp2 = shade,
        species_name_sp2 = species_name
      ) |>
      select(species_id_old, species_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  filter(name == 'BA') |>
  select(species_name_sp1, ND, RFD, shade_sp2) |>
  pivot_longer(
    cols = c(ND, RFD)
  ) |>
  mutate(
    name = case_match(
      name,
      'ND' ~ 'Niche difference',
      'RFD'~ 'Relative fitness difference'
    )
  ) |>
  ggplot(aes(value, species_name_sp1, fill = shade_sp2)) +
    ggridges::geom_density_ridges2(
      # color = rgb(0.3,0.5,0.4,0.6),
      # fill = rgb(0.3,0.5,0.4,0.6),
      color = 'transparent',
      alpha = 0.8
    ) +
    facet_grid(~name, scales = 'free') +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab('') +
    theme(
      axis.text.y = element_text(face = "italic")
    ) +
    labs(fill = 'shade tolerance\nof resident species')

coex |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp1 = shade) |>
      select(species_id_old, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp2 = shade) |>
      select(species_id_old, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    shade_sp1 = factor(shade_sp1, levels = c('tolerant', 'intermediate', 'intolerant')),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade_sp1, RFD, fill = shade_sp2)) +
    geom_boxplot() +
    theme_minimal() +
    xlab('Shade tolerance of invader species') +
    ylab('Relative fitness difference') +
    labs(fill = 'Shade tolerance\nof resident species')


# How does RFD relates to lambda_sp1/lambda_sp2
coex |>
  filter(name == 'BA') |>
  full_join(
    out |>
      filter(subsim == 'subsim1') |>
      select(sim, lambda, sp) |>
      pivot_wider(
        names_from = sp,
        values_from = lambda
      ) |>
      mutate(
        lambda_diff = sp1/sp2
      ) |>
      select(!contains('sp')) |>
      left_join(
        pars_sim |>
          mutate(
            sim = row_number(),
            sp_pair = paste0(sp1, '_', sp2)
          ) |>
          select(
            !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
          )
      ) |>
      select(!sim)
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp1 = shade) |>
      select(species_id_old, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp2 = shade) |>
      select(species_id_old, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  ggplot(aes(RFD, lambda_diff, color = shade_sp1)) +
    geom_point(alpha = 0.5, size = 0.3) +
    theme_minimal()





# Coexistence using the formulas from Lyu 2023 EcoLet
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# source: https://github.com/ShengmanLyu/Compensatory_Responses_and_species_coexistence/blob/b5c4010c7dff7b6f3c696d2ab8948afba2da338a/Functions.R#L254
sensitivity.igr <- function(pgr.intrinsic, pgr.invasion) {
  # pgr.intrinsic: intrinsic growth rate of the focal species
  # pgr.invasion: invasion growth rate of the focal species

  # log
  pgr.intrinsic.log <- log(pgr.intrinsic)
  pgr.invasion.log <- log(pgr.invasion)
  
  # calculate sensitivity
  if(pgr.intrinsic.log > 0) {
    Sens <- 1- pgr.invasion.log/pgr.intrinsic.log
  }
  else if(pgr.intrinsic.log < 0) {
    Sens <- pgr.invasion.log/pgr.intrinsic.log - 1
  }
  # output
  return(Sens)
}

sens_l <- out |>
  filter(sp == 'sp1') |>
  select(sim, subsim, lambda) |>
  pivot_wider(
    names_from = subsim,
    values_from = lambda
  ) |>
  rowwise() |>
  mutate(
    sens_sp1 = sensitivity.igr(subsim1, subsim2)
  ) |>
  ungroup() |>
  #/!\ Important filter here /!\
  filter(
    sens_sp1 >= 0 & # this happens when theta is positive
    sens_sp1 < quantile(sens_sp1, probs = 0.95)
  ) |>
  select(!contains('subsim')) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  )
# modified from: https://github.com/ShengmanLyu/Compensatory_Responses_and_species_coexistence/blob/b5c4010c7dff7b6f3c696d2ab8948afba2da338a/Functions.R#L330
coex.ndfd <- function(s12, s21) {
  # s12: sensitivity of species 1 invading species 2
  s <- c(s12, s21)
  # ND
  nd <- 1 - sqrt(s12*s21)
  names(nd) <- "ND"
  
  # FD
  # afd <- exp( sqrt( (log(s12/gm) + log(s21/gm))/2 ) )
  # afd <- exp(sqrt(mean(log(s)**2) - (mean(log(s)))**2))
  rfd <- sqrt(s21/s12)
  names(rfd)  <- "RFD"
  
  # superior
  if(rfd > 1) superior <- "sps1"
  else superior <- "sps2"
  
  # Coexistence outcomes
  # make sure FD always > 1
  if(rfd < 1) rfd2 <- 1/rfd
  else rfd2 <- rfd
  
  # coexistence metrics
  coex.metric <- 1/ (rfd2*(1-nd))
  names(coex.metric) <- "Coexistence.metric"
  
  if(rfd2 <= 1/(1-nd)) { 
    outcome <- "Coexistence"
    winner <- "sp1_sps2"
    }
  else if(rfd > 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps1"
    }
  else if(rfd < 1) {
    outcome <- "Competitive exclusion"
    winner <- "sps2"
  }else if(rfd == 1) {
    outcome <- "Coexistence"
    winner <- "sp1_sps2"
  }
  return(c(nd, rfd, coex.metric, outcome = outcome, winner = winner, superior = superior) )
}


coex_l <- sens_l |>
  full_join(
    sens_l |>
      rename(sens_sp2 = sens_sp1),
    by = c(
      'replication' = 'replication',
      'sp2' = 'sp1', 'sp1' = 'sp2'
    )
  ) |>
  # remove NAs because I removed negative sensitivities
  drop_na() |>
  rowwise() |>
  mutate(
    coex = list(coex.ndfd(sens_sp1, sens_sp2))
  ) |>
  mutate(var = list(names(coex))) |>
  unnest(c(var, coex)) |>
  pivot_wider(
    names_from = var,
    values_from = coex
  ) |>
  mutate_at(c('ND', 'RFD', 'Coexistence.metric'), as.numeric) |>
  ungroup() |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        sp_name_sp1 = species_name,
        shade_sp1 = shade
      ) |>
      select(species_id_old, sp_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        sp_name_sp2 = species_name,
        shade_sp2 = shade
      ) |>
      select(species_id_old, sp_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  select(!c(replication, sp1, sp2))



# Figures
#################################

# first figure to help have an intuition of sensitivity and niche overlap
expand_grid(
  s12 = seq(0.05, 3, 0.05),
  s21 = seq(0.05, 3, 0.05)
) |>
rowwise() |>
mutate(
  coex = list(coex.ndfd(s12, s21))
) |>
mutate(var = list(names(coex))) |>
unnest(c(var, coex)) |>
pivot_wider(
  names_from = var,
  values_from = coex
) |>
mutate_at(c('ND', 'RFD', 'Coexistence.metric'), as.numeric) |>
mutate(
  RFD2 = ifelse(RFD < 1, 1/RFD, RFD)
) |>
ggplot(aes(s12, s21)) +
  geom_tile(aes(fill = log(1-ND))) +
  # geom_tile(aes(fill = log(RFD))) +
  # geom_tile(aes(fill = log(1/(RFD2 * (1-ND))))) +
  # geom_tile(aes(fill = log(Coexistence.metric)))) +
  scale_fill_gradient2()


# Functions to create coexistence polygon
# polygan of coexistence area
polygan.coex.no = function(x1 = -1, x2=0) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = seq(x2,x1,length.out = 100)
  xx3 = rep(x1,100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(-xx1, xx2, seq(x1,-x1, length.out=100))
  return(data.frame(x=xx,y=yy))
}

polygan.prio.no = function(x1 = 0, x2=1) {
  xx1 = seq(x1,x2, length.out = 100)
  xx2 = rep(x2,100)
  xx3 = seq(x2,x1,length.out = 100)
  xx <- c(xx1,xx2,xx3)
  yy <- c(xx1,seq(x2,-x2, length.out=100), -xx3)
  return(data.frame(x=xx,y=yy))
}

ggplot(polygan.coex.no(x1 = -6), aes(x=x,y=y)) +
  geom_polygon(alpha = 0.8) +
  geom_polygon(
    data = polygan.prio.no(),
    aes(x=x,y=y),
    alpha = 0.8
  ) +
  geom_point(
    data = coex_l |>
      mutate(sp_pair = paste0(sp_name_sp1, sp_name_sp2)) |>
      group_by(sp_pair) |>
      mutate(
        RFD2 = ifelse(RFD < 1, 1/RFD, RFD),
        NO = mean(log(1-ND)),
        RFD = mean(log(RFD))
      ) |>
      slice_head(n = 1),
    aes(NO, RFD, color = shade_sp1),
    alpha = 0.6,
    size = 0.9
  ) +
  theme_minimal() +
  xlab('log(Niche overlap)') +
  ylab('log(Relative fitness difference)') +
  labs(color = 'Shade tolerance of\ninvasive species')

coex_l |>
  ggplot(aes(sp_name_sp2, sp_name_sp1, fill = sens_sp1)) +
    geom_raster() +
    scale_fill_gradient2() +
    theme_classic() +
    labs(
      fill = ''
    ) +
    theme(
      axis.text.x = element_text(
        size = 9, face = 'italic', angle = 45, hjust = 1, vjust = 1
      ),
      axis.text.y = element_text(size = 9, face = 'italic'),
      plot.title = element_markdown()
    ) +
    labs(subtitle = 'Intensity of coexistence (blue) vs competitive exclusion (red) between pair of species') +
    xlab('Resident species') +
    ylab('Invader species')

coex_l |>
  ggplot(
    aes(log(1-ND), fct_reorder(sp_name_sp1, log(1-ND)), fill = shade_sp2)
  ) +
  ggridges::geom_density_ridges2(
      # color = rgb(0.3,0.5,0.4,0.6),
      # fill = rgb(0.3,0.5,0.4,0.6),
      color = 'transparent',
      alpha = 0.8
    ) +
    #facet_grid(~name, scales = 'free') +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab('log(Niche overlap)') +
    theme(
      axis.text.y = element_text(face = "italic")
    ) +
    labs(fill = 'shade tolerance\nof resident species')

coex_l |>
  select(c(contains('sp_'), shade_sp2, ND, RFD)) |>
  mutate(
    ND = 1 - ND
  ) |>
  mutate(across(c('ND', 'RFD'), log)) |>
  rename(
    `Niche overlap` = ND,
    `Relative fitness difference` = RFD
  ) |>
  pivot_longer(
    cols = c(`Niche overlap`, `Relative fitness difference`)
  ) |>
  ggplot(
    aes(value, fct_reorder(sp_name_sp1, value), fill = shade_sp2)
  ) +
  ggridges::geom_density_ridges2(
      # color = rgb(0.3,0.5,0.4,0.6),
      # fill = rgb(0.3,0.5,0.4,0.6),
      color = 'transparent',
      alpha = 0.8
    ) +
    geom_vline(xintercept = 0, alpha = 0.9, linetype = 2) +
    facet_grid(~name, scales = 'free') +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab('log(value)') +
    theme(
      axis.text.y = element_text(face = "italic")
    ) +
    labs(fill = 'shade tolerance\nof resident species')


coex_l |>
  mutate(
    shade_sp1 = factor(shade_sp1, levels = c('tolerant', 'intermediate', 'intolerant')),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade_sp1, log(RFD), fill = shade_sp2)) +
    geom_boxplot() +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_minimal() +
    xlab('Shade tolerance of invader species') +
    ylab('log(Relative fitness difference)') +
    labs(fill = 'Shade tolerance of\nresident species')

coex_l |>
  mutate(
    shade_sp1 = factor(shade_sp1, levels = c('tolerant', 'intermediate', 'intolerant')),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade_sp1, log(1 - ND), fill = shade_sp2)) +
    geom_boxplot() +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_minimal() +
    xlab('Shade tolerance of invader species') +
    ylab('log(Niche overlap)') +
    labs(fill = 'Shade tolerance of\nresident species')

# From the above figure:
# - RFD is more sensitive to shade tollerance than ND
# - RFD shows quite well the importance of tolerance to shade betwen the 
# invader and the resident species to define the fitness ratio
# - Niche Overlap (1-ND) is quite constant among shade tolerance combinations,
# exept for the intolerant vs intolerant pair in which niche overlap increases
# The question is what is driving then the difference in ND?

library(ranger)

# extract parameters used in each simulation
pars_sp1 <- parse_number(
  dir(
    'simulations/coexistence/output/',
    pattern = 'RDS'
  )
) |>
sort() |>
map_dfr(
  ~ readRDS(
    paste0(
      'simulations/coexistence/output/sim_', .x, '.RDS'
    )
  )[['pars_sp1']] |>
  unlist() |>
  enframe() |>
  bind_cols(sim = .x) #|>
  #filter(name %in% c('growth.Beta', 'growth.theta', 'mort.Beta', 'mort.theta', 'rec.beta_p', 'rec.optimal_BA', 'rec.sigma_BA'))
)

rf_dat <- coex_l |>
  select(!c(sens_sp2, outcome, winner, superior, shade_sp1, shade_sp2)) |>
  left_join(
    pars_sp1 |>
      rename(sim.x = sim) |>
      pivot_wider()      
      # mutate(
      #   growth.BetaInter = growth.Beta * growth.theta,
      #   mort.BetaInter = mort.Beta * mort.theta,
      #   growth.theta = NULL,
      #   mort.theta = NULL
      # )
  ) |>
  select(!c(sim.x, sim.y))

rf_sensSp1 <- ranger(
  sens_sp1 ~ ., 
  data = rf_dat |>
    # select(!c(ND, RFD, Coexistence.metric)) |>
    select(c(sens_sp1, sp_name_sp1, sp_name_sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)),
  importance = 'impurity_corrected'
)

knitr::kable(sort(importance(rf_sensSp1)/max(importance(rf_sensSp1))))

rf_ND <- ranger(
  log(1-ND)~ ., 
  data = rf_dat |>
    # select(!c(sens_sp1, RFD, Coexistence.metric)),
    select(c(ND, sp_name_sp1, sp_name_sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)),
  importance = 'impurity_corrected'
)

knitr::kable(sort(importance(rf_ND)/max(importance(rf_ND))))

rf_RFD <- ranger(
  log(RFD) ~ ., 
  data = rf_dat |>
    # select(!c(sens_sp1, ND, Coexistence.metric)),
    select(c(RFD, sp_name_sp1, sp_name_sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)),
  importance = 'impurity_corrected'
)

knitr::kable(sort(importance(rf_RFD)/max(importance(rf_RFD))))


rf_coex <- ranger(
  log(Coexistence.metric) ~ ., 
  data = rf_dat |>
    # select(!c(sens_sp1, ND, Coexistence.metric)),
    select(c(Coexistence.metric, sp_name_sp1, sp_name_sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)),
  importance = 'impurity_corrected'
)

knitr::kable(sort(importance(rf_coex)/max(importance(rf_coex))))























#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Just a test:
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
    a_ij = log(subsim2/subsim1)
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
    geom_raster(aes(fill = a_ij)) + 
    scale_fill_gradient2() +
    theme_classic() +
    labs(
      fill = 'a_ij'
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
    a_ij = log(subsim2/subsim1)
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
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        species_name_sp1 = species_name,
        shade_sp1 = shade
      ) |>
      select(species_id_old, species_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        species_name_sp2 = species_name,
        shade_sp2 = shade) |>
      select(species_id_old, species_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  ggplot(aes(a_ij, fct_reorder(species_name_sp1, a_ij), fill = shade_sp2)) +
    ggridges::geom_density_ridges2(
      color = 'transparent',
      alpha = 0.7
    ) +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab('a_ij') +
    theme(
      axis.text.y = element_text(face = "italic"),
      legend.position = 'top'
    ) +
    labs(
      fill = 'Shade tolerance\nof resident species'
    )


# Coexistence metrics
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
    a_ij = log(subsim2/subsim1),
    a_ii = 1/subsim1
  ) |>
  select(sim, a_ii, a_ij) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number()
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, seed)
      )
  ) |>
  select(!sim) |>
  full_join(
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
        a_ji = log(subsim2/subsim1),
        a_jj = 1/subsim1
      ) |>
      select(sim, a_jj, a_ji) |>
      left_join(
        pars_sim |>
          mutate(
            sim = row_number()
          ) |>
          select(
            !c(n_time, deltaTime, param_method, plotSize, seed)
          )
      ) |>
      select(!sim),
    by = c(
      'replication' = 'replication',
      'sp2' = 'sp1', 'sp1' = 'sp2'
    )
  ) |>
  mutate(
    gamma_ij = a_ij/a_jj,
    gamma_ji = a_ji/a_ii,
    FD = sqrt(gamma_ij/gamma_ji),
    p = sqrt(gamma_ij * gamma_ji)
  )




# Sensitivity of resident species at equilibrium after the invasion of the focal species
sens_sp2 <- out |>
  filter(sp == 'sp2') |>
  group_by(sim) |>
  mutate(
    diff_N = (N_eq[subsim == 'subsim1'] - N_eq[subsim == 'subsim2'])/N_eq[subsim == 'subsim1'],
    diff_BA = (BA_eq[subsim == 'subsim1'] - BA_eq[subsim == 'subsim2'])/BA_eq[subsim == 'subsim1']
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

sens_sp2 |>
  ggplot(aes(diff_N, diff_BA, color = sp2)) +
    geom_point(alpha = 0.5, size = 0.5) +
    theme(legend.position = 'none')

sens_sp2 |>
  left_join(
    read_csv('data/species_id.csv') |>
      select(!species_id_new),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    species_name = fct_reorder(species_name, diff_BA)
  ) |>
  ggplot(aes(diff_BA, species_name)) +
    ggridges::geom_density_ridges2(color = rgb(0.3,0.5,0.4,0.6), fill = rgb(0.3,0.5,0.4,0.6), alpha = 0.8) +
    theme_classic() +
    ylab('Resident species') +
    xlab('Basal area sensitivity') +
    theme(
      axis.text.y = element_text(face = "italic")
    )






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate stochastic growth rate from N and BA over time
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambda_s <- function(array_id)
{
  sim_i <- readRDS(
    paste0('simulations/coexistence/output/sim_', array_id, '.RDS')
  )

  lambda_sp1_subsim1 <- sim_i[['year_summ_subsim1']] |>
    mutate(
      N_sp1_t0 = lag(N_sp1),
      BA_sp1_t0 = lag(BA_sp1),
      lambda_N = exp(mean(log(N_sp1/N_sp1_t0), na.rm = TRUE)),
      lambda_BA = exp(mean(log(BA_sp1/BA_sp1_t0), na.rm = TRUE))
    ) |>
    head(n = 1) |>
    pull(lambda_N, lambda_BA)
  
    lambda_sp2_subsim1 <- sim_i[['year_summ_sp2_subsim2']] |>
      mutate(
        N_sp2_t0 = lag(N_sp2),
        BA_sp2_t0 = lag(BA_sp2),
        lambda_N = exp(mean(log(N_sp2/N_sp2_t0), na.rm = TRUE)),
        lambda_BA = exp(mean(log(BA_sp2/BA_sp2_t0), na.rm = TRUE))
      ) |>
      head(n = 1) |>
      pull(lambda_N, lambda_BA)
    
    lambda_sp1sp2_subsim2 <- sim_i[['year_summ_sp1sp2_subsim2']] |>
      mutate(
        N_sp1_t0 = lag(N_sp1),
        N_sp2_t0 = lag(N_sp2),
        BA_sp1_t0 = lag(BA_sp1),
        BA_sp2_t0 = lag(BA_sp2),
        lambda_sp1_N = exp(mean(log(N_sp1/N_sp1_t0), na.rm = TRUE)),lambda_sp2_N = exp(mean(log(N_sp2/N_sp2_t0), na.rm = TRUE)),
        lambda_sp1_BA = exp(mean(log(BA_sp1/BA_sp1_t0), na.rm = TRUE)),
        lambda_sp2_BA = exp(mean(log(BA_sp2/BA_sp2_t0), na.rm = TRUE))
      ) |>
      head(n = 1)
      pull(lambda_sp1_N, lambda_sp2_N, lambda_sp1_BA, lambda_sp2_BA)

}


l_stoch <- parse_number(
    dir('simulations/coexistence/output/',
    pattern = 'RDS')
  ) |>
  sort() |>
  map_dfr(
    ~ lambda_s(.x) |>
    bind_cols(sim = .x)
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
  as_tibble()




l_stoch |>  
  filter(sp_pair == '18032ABIBAL_18034PICRUB') |>  
  ggplot(aes(lambda_N, lambda_BA, color = as.factor(replication))) +
    geom_path(alpha = 0.5) + 
    theme(legend.position = 'none')


# calculate variance over time
var_lambda <- l_stoch |>
  group_by(sp_pair, year) |>
  reframe(
    var_N = var(lambda_N),
    var_BA = var(lambda_BA)
  )







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check difference in size distribution of the resident species
# between alone and with competitive species
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_sizeDist <- function(array_id) 
{
  # read sim
  sim_i <- readRDS(
    paste0('simulations/coexistence/output/sim_', array_id, '.RDS')
  )

  sz_dist <- sim_i[['mat_sp1_subsim1']]

  tibble(
    bin = 1:length(sz_dist),
    size = sz_dist
  )
}


sz <- parse_number(dir('simulations/coexistence/output/', pattern = 'RDS')) |>
  sort() |>
  map_dfr(
    ~ get_sizeDist(.x) |>
    bind_cols(sim = .x)
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  ) 

sz |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp2 = shade) |>
      select(species_id_old, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  filter(shade_sp2 %in% c('intolerant', 'tolerant')) |>
  group_by(sp1, sp2, bin) |>
  reframe(
    size = mean(size),
    shade_sp2 = unique(shade_sp2)
  ) |>
  ggplot(aes(bin, size, color = shade_sp2, group = sp2)) +
    geom_line(alpha = 0.6) +
    theme_minimal() +
    facet_wrap(~sp1, scales = 'free')


get_sizeDistDiff <- function(array_id) 
{
  # read sim
  sim_i <- readRDS(
    paste0('simulations/coexistence/output/sim_', array_id, '.RDS')
  )

  sz_dist <- (sim_i[['mat_sp2_subsim2']] - sim_i[['mat_sp1sp2_subsim2']])^2

  tibble(
    bin = 1:length(sz_dist),
    size = sz_dist
  )
}


sz_diff <- parse_number(dir('simulations/coexistence/output/', pattern = 'RDS')) |>
  sort() |>
  map_dfr(
    ~ get_sizeDistDiff(.x) |>
    bind_cols(sim = .x)
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  )

sz_diff |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp1 = shade) |>
      select(species_id_old, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp2 = shade) |>
      select(species_id_old, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  group_by(sp2, bin) |>
  reframe(
    size = mean(size),
    shade_sp2 = unique(shade_sp2)
  ) |>
  ggplot(aes(bin, size, color = shade_sp2, group = sp2)) +
    geom_line(alpha = 0.6) +
    theme_minimal()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sensitivity vs competition parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


pars_sp1 <- parse_number(
  dir(
    'simulations/coexistence/output/',
    pattern = 'RDS'
  )
) |>
sort() |>
map_dfr(
  ~ readRDS(
    paste0(
      'simulations/coexistence/output/sim_', .x, '.RDS'
    )
  )[['pars_sp1']] |>
  unlist() |>
  enframe() |>
  bind_cols(sim = .x) #|>
  #filter(name %in% c('growth.Beta', 'growth.theta', 'mort.Beta', 'mort.theta', 'rec.beta_p', 'rec.optimal_BA', 'rec.sigma_BA'))
)

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
    sens = (subsim1 - subsim2)/subsim1,
    subsim1 = NULL,
    subsim2 = NULL
  ) |>
  left_join(
    pars_sp1 |>
      pivot_wider() |>
      mutate(
        growth.BetaInter = growth.Beta * growth.theta,
        mort.BetaInter = mort.Beta * mort.theta,
        growth.theta = NULL,
        mort.theta = NULL
      )
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(shade_sp1 = shade) |>
      select(species_id_old, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  pivot_longer(
    cols = c(a_ii, a_ij, sens),
    names_to = 'sens',
    values_to = 'sens_value'
  ) |>
  pivot_longer(
    cols = c(contains('.')),
    names_to = 'par',
    values_to = 'par_value'
  ) |>
  ggplot(aes(sens_value, par_value, color = shade_sp1)) +
    #geom_point(alpha = 0.5, size = 0.4) +
    geom_smooth() +
    facet_wrap(~sens + par, scales = 'free') +
    theme_minimal()



library(ranger)

rf_dat <- out |>
  filter(sp == 'sp1') |>
  select(sim, subsim, lambda, BA_eq) |>
  pivot_wider(
    names_from = subsim,
    values_from = c(lambda, BA_eq)
  ) |>
  mutate(
    # first deal with 0 BA or values too close to it
    BA_eq_subsim2 = ifelse(BA_eq_subsim2 < 0.1, 0.1, BA_eq_subsim2),
    BA_eq_subsim1 = ifelse(BA_eq_subsim1 < 0.1, 0.1, BA_eq_subsim1),
    a_ii = 1/BA_eq_subsim1,
    sens_BA = (BA_eq_subsim1 - BA_eq_subsim2)/BA_eq_subsim1,
    sens_lambda = (lambda_subsim1 - lambda_subsim2)/lambda_subsim1,
    BA_eq_subsim1 = NULL,
    BA_eq_subsim2 = NULL,
    lambda_subsim1 = NULL,
    lambda_subsim2 = NULL
  ) |>
  left_join(
    pars_sp1 |>
      pivot_wider() #|>
      # mutate(
      #   growth.BetaInter = growth.Beta * growth.theta,
      #   mort.BetaInter = mort.Beta * mort.theta,
      #   growth.theta = NULL,
      #   mort.theta = NULL
      # )
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  ) |>
  select(!c(sim, replication))

rf_aii <- ranger(a_ii ~., rf_dat |> select(c(a_ii, sp1, sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)), importance = 'impurity_corrected')
rf_sensBA <- ranger(sens_BA ~., rf_dat |> select(c(sens_BA, sp1, sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)), importance = 'impurity_corrected')
rf_sensLambda <- ranger(sens_lambda ~., rf_dat |> select(c(sens_lambda, sp1, sp2, growth.Beta, growth.theta, mort.Beta, mort.theta, rec.beta_p, rec.optimal_BA, rec.sigma_BA)), importance = 'impurity_corrected')

knitr::kable(sort(importance(rf_aii)/max(importance(rf_aii))))
knitr::kable(sort(importance(rf_sensBA)/max(importance(rf_sensBA))))
knitr::kable(sort(importance(rf_sensLambda)/max(importance(rf_sensLambda))))




# Another random test
# If we standarize the sensitivity of focal species by the abundance of the resident? 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

p1 <- out |>
  mutate(BA_eq = ifelse(BA_eq < 0.5, 0, BA_eq)) |>
  select(sim, subsim, sp, BA_eq) |>
  pivot_wider(
    names_from = c(sp, subsim),
    values_from = BA_eq,
    names_sep = '_'
  ) |>
  mutate(
    s12 = log(sp1_subsim2/sp1_subsim1),
    s12s = s12/sp2_subsim2,
  ) |>
  filter(!is.infinite(s12s)) |>
  filter(
    s12s < quantile(s12s, 0.995, na.rm = T) &
    s12s > quantile(s12s, 0.005, na.rm = T)
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp1 = shade,
        sp_name_sp1 = species_name
      ) |>
      select(species_id_old, sp_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp2 = shade,
        sp_name_sp2 = species_name
      ) |>
      select(species_id_old, sp_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant')),
  ) |>
  ggplot(aes(s12s, fct_reorder(sp_name_sp1, s12s), fill = shade_sp2)) +
    ggridges::geom_density_ridges2(
      color = 'transparent',
      alpha = 0.7
    ) +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab(expression(alpha[ij])) +
    theme(
      axis.text.y = element_text(face = "italic"),
      legend.position = 'top'
    ) +
    labs(
      fill = 'Shade tolerance\nof resident species'
    )
 
 p2 <- out |>
  mutate(BA_eq = ifelse(BA_eq < 0.5, 0, BA_eq)) |>
  select(sim, subsim, sp, BA_eq) |>
  pivot_wider(
    names_from = c(sp, subsim),
    values_from = BA_eq,
    names_sep = '_'
  ) |>
  mutate(
    s12 = log(sp1_subsim2/sp1_subsim1),
    s12s = s12/sp2_subsim2,
  ) |>
  filter(!is.infinite(s12s)) |>
  filter(
    s12s < quantile(s12s, 0.995, na.rm = T) &
    s12s > quantile(s12s, 0.005, na.rm = T)
  ) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number(),
        sp_pair = paste0(sp1, '_', sp2)
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, sp_pair, seed)
      )
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp1 = shade,
        sp_name_sp1 = species_name
      ) |>
      select(species_id_old, sp_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp2 = shade,
        sp_name_sp2 = species_name
      ) |>
      select(species_id_old, sp_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    shade_sp1 = factor(shade_sp1, levels = c('tolerant', 'intermediate', 'intolerant')),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(shade_sp1, s12s)) +
    geom_boxplot() +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_minimal() +
    xlab('Shade tolerance of invader species') +
    ylab(expression(alpha[ij])) +
    labs(fill = 'Shade tolerance of\nresident species')

ggpubr::ggarrange(p1, p2, ncol = 2)



out |>
  filter(subsim == 'subsim1' & sp == 'sp1') |>
  select(sim, BA_eq) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number()
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, seed)
      )
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp1 = shade,
        sp_name_sp1 = species_name
      ) |>
      select(species_id_old, sp_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp2 = shade,
        sp_name_sp2 = species_name
      ) |>
      select(species_id_old, sp_name_sp2, shade_sp2),
    by = c('sp2' = 'species_id_old')
  ) |>
  mutate(
    shade_sp1 = factor(shade_sp1, levels = c('tolerant', 'intermediate', 'intolerant')),
    shade_sp2 = factor(shade_sp2, levels = c('tolerant', 'intermediate', 'intolerant'))
  ) |>
  ggplot(aes(BA_eq, fct_reorder(sp_name_sp1, BA_eq), fill = shade_sp2)) +
    ggridges::geom_density_ridges2(
      color = 'transparent',
      alpha = 0.7
    ) +
    scale_fill_manual(
      values = c("#87bc45", "#edbf33", "#ea5545")
    ) +
    theme_classic() +
    ylab('Invader species') +
    xlab(expression(alpha[ij])) +
    theme(
      axis.text.y = element_text(face = "italic"),
      legend.position = 'top'
    ) +
    labs(
      fill = 'Shade tolerance\nof resident species'
    )
 


out |>
  filter(sp == 'sp1') |>
  select(sim, lambda, subsim, BA_eq) |>
  left_join(
    pars_sim |>
      mutate(
        sim = row_number()
      ) |>
      select(
        !c(n_time, deltaTime, param_method, plotSize, seed)
      )
  ) |>
  group_by(sp1, sp2, subsim) |>
  reframe(
    lambda = mean(lambda),
    BA_eq = mean(BA_eq)
  ) |>
  left_join(
    read_csv('data/species_id.csv') |>
      rename(
        shade_sp1 = shade,
        sp_name_sp1 = species_name
      ) |>
      select(species_id_old, sp_name_sp1, shade_sp1),
    by = c('sp1' = 'species_id_old')
  ) |>
  ggplot(aes(lambda, BA_eq, color = subsim)) +
    geom_point(alpha = 0.8, size = 1.3) +
    facet_wrap(~shade_sp1) +
    xlim(0.99, 1.09)
