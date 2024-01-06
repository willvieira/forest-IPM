# Output preparation and viz

library(tidyverse)
library(ggdist)
library(ggtext)
library(ggrepel)
library(ggpubr)
library(ggiraph)
library(fixest)
library(ggblend)
library(ggdensity)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load stuff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load parameters and output
sim_path <- file.path('simulations', 'model_lambdaPlot')

# parameters
sim_pars <- readRDS(file.path(sim_path, 'simulation_pars.RDS'))

# output
out_pars <- readRDS(file.path(sim_path, 'param_posterior.RDS')) |>
  left_join(sim_pars) |>
  mutate(
    sim = factor(sim),
    cond = factor(cond)
  )


# data for latitude - longitude
treeData <- readRDS(
    paste0(readLines('_data.path'), 'treeData.RDS')
  )

treeData |>
  filter(species_id %in% unique(sim_pars$species_id)) |>
  group_by(species_id) |>
  reframe(
    mid_pos = (max(bio_01_mean, na.rm=TRUE) + min(bio_01_mean, na.rm=TRUE))/2
  ) ->
  sp_pos


# load species info
spIds <- read_csv(file.path(readLines('_data.path'), 'species_id.csv')) |>
  mutate(
    shade_sylvics = factor(shade_sylvics, levels = c(
      'very-tolerant',
      'tolerant',
      'intermediate',
      'intolerant',
      'very-intolerant'
      )
    )
  ) |>
  filter(species_id_old %in% unique(sim_pars$species_id)) |>
  select(species_id_old, species_name) |>
  rename(species_id = species_id_old)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter distribution in function of rhat
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

readRDS(file.path(sim_path, 'out_diag.RDS')) |>
  imap_dfr(~.x[['rhat']] |> bind_cols(array_id = .y)) ->
rhats


out_pars |>
  left_join(
    sim_pars |>
      mutate(
        sim = factor(sim),
        cond = factor(cond)
      )
  ) |>
  left_join(rhats |> rename(par = variable)) |>
  filter(par != 'lp__') ->
out_pars



out_pars |>
  # filter(rhat < 1.03) |>
  filter(par == 'sigma_inter' & cond == 2) |>
  ggplot() +
  aes(value, fct_reorder(species_id, value)) +
  aes(color = sim) +
  facet_wrap(~border) +
  stat_pointinterval() +
  theme_classic() +
  geom_vline(xintercept = 0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Traceplot for each model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pdf(file.path(sim_path, 'traceplot.pdf'), height = 3, width = 12)
  for(Sp in unique(sim_pars$species_id)) {
    for(brd in c('cold', 'hot')) {
      for(cnd in 1:2) {
        # try({
          print(
          out_pars |>
          filter(species_id == Sp & border == brd & sim == 2 & cond == cnd) |>
          ggplot() +
          aes(iter, value) +
          facet_wrap(~par, scales = 'free', nrow = 1) +
          geom_line() +
          labs(
            title = Sp,
            subtitle = paste0(brd, '; condition: ', cnd)
          )
        )
        # })
      }
    }
  }
  dev.off()

#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Raw predicted distribution for each species and parameter and sim and cond
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  sim_lambda <- readRDS('simulations/lambda_plot/simulation_pars.RDS') |>
    ungroup() |>
    mutate(array_id = row_number())

  pdf(file.path(sim_path, 'pred_lm_cond2.pdf'), height = 6, width = 9)
  for(Sp in unique(sim_pars$species_id))
  {
    try({
    Cond = 2

    # get which db_origin is more abundant
    db_sim_cold = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 1, '_', Cond, '.RDS')) |>
      bind_rows(
        readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', Cond, '.RDS'))
      ) |>
      mutate(sim = as.factor(sim))

    db_sim_hot = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 1, '_', Cond, '.RDS')) |>
      bind_rows(
        readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', Cond, '.RDS'))
      ) |>
      mutate(sim = as.factor(sim))

    # get temperature range for X axis limits
    db_sim_cold |>
      filter(sim == 1) |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_cold

    db_sim_hot |>
      filter(sim == 1) |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_hot

    # define Y axis limits
    db_sim_hot |>
      filter(sim == 2) |>
      bind_cols(border = 'hot') |>
      bind_rows(
        db_sim_cold |>
          filter(sim == 2) |>
          bind_cols(border = 'cold')
      ) |>
      pull(lambda) |>
      log() |>
      quantile(probs = c(0.001, 0.999)) ->
    yLim
    
    out_pars |>
      filter(species_id == Sp & border == 'cold' & cond == Cond) |>
      select(sim, iter, par, value) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(sim, iter) |>
      expand_grid(
        temp = seq(temp_range_cold[1], temp_range_cold[2], length.out = 50)
      ) |>
      mutate(
        lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
      ) |>
      ggplot() +
      aes(temp, lambda) +
      facet_wrap(~sim, ncol = 1) +
      geom_point(
        data = db_sim_cold |>
          group_by(sim) |>
          slice_sample(n = 1e4),
        aes(temp, log(lambda)),
        alpha = 0.4, size = 0.3
      ) +
      stat_lineribbon(.width = 0.9, alpha = 0.6, color = 'purple') +
      geom_hline(yintercept = 0, linetype = 2) +
      theme_classic() +
      theme(
        plot.title = element_text(face = 'italic'),
        legend.position = 'none'
      ) +
      labs(
        title = spIds |> filter(species_id == Sp) |> pull(species_name),
        subtitle = 'Cold',
        x = 'Mean annual temperature (°C)',
        y = expression('ln('~lambda~')')
      ) +
      xlim(temp_range_cold) +
      ylim(yLim) ->
    p1

    out_pars |>
      filter(species_id == Sp & border == 'hot' & cond == Cond) |>
      select(sim, iter, par, value) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(sim, iter) |>
      expand_grid(
        temp = seq(temp_range_hot[1], temp_range_hot[2], length.out = 50)
      ) |>
      mutate(
        lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
      ) |>
      ggplot() +
      aes(temp, lambda) +
      facet_wrap(~sim, ncol = 1) +
      geom_point(
        data = db_sim_hot |>
          group_by(sim) |>
          slice_sample(n = 1e4),
        aes(temp, log(lambda)),
        alpha = 0.4, size = 0.3
      ) +
      stat_lineribbon(.width = 0.9, alpha = 0.6, color = 'purple') +
      geom_hline(yintercept = 0, linetype = 2) +
      theme_classic() +
      theme(legend.position = 'none') +
      labs(
        title = '',
        subtitle = 'Hot',
        x = 'Mean annual temperature (°C)',
        y = NULL
      ) +
      xlim(temp_range_hot) +
      ylim(yLim) ->
    p2

    print(ggarrange(p1, p2, ncol = 2))
  })
  }
  dev.off()

#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Raw predicted distribution for each species and parameter and sim and cond
# But here only for sim == 2 and each chain separetly (for indiv checks)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pdf(file.path(sim_path, 'pred_lm_chain.pdf'), height = 6, width = 9)
  for(Sp in unique(sim_pars$species_id))
  {
    try({

    for(Cond in 1:2) {
      # get which db_origin is more abundant
      db_sim_cold = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', Cond, '.RDS'))

      db_sim_hot = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', Cond, '.RDS'))

      # get temperature range for X axis limits
      db_sim_cold |>
        filter(sim == 2) |>
        pull(bio_01_mean) |>
        range(na.rm = TRUE) ->
      temp_range_cold

      db_sim_hot |>
        filter(sim == 2) |>
        pull(bio_01_mean) |>
        range(na.rm = TRUE) ->
      temp_range_hot

      print(
      out_pars |>
        filter(species_id == Sp & border == 'cold' & sim == 2 & cond == Cond) |>
        select(iter, par, value) |>
        pivot_wider(names_from = par, values_from = value) |>
        mutate(
          chain = as.factor(case_when(
            iter <= 1000 ~ 1,
            iter <= 2000 ~ 2,
            iter <= 3000 ~ 3,
            iter > 3000 ~ 4
          )
        )) |>
        group_by(iter) |>
        expand_grid(
          temp = seq(temp_range_cold[1], temp_range_cold[2], length.out = 50)
        ) |>
        mutate(
          lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
        ) |>
        ggplot() +
        aes(temp, lambda) +
        facet_wrap(~chain) +
        geom_point(
          data = db_sim_cold |>
            slice_sample(n = 1e4),
          aes(temp, log(lambda)),
          alpha = 0.4, size = 0.3
        ) +
        stat_lineribbon(.width = 0.9, alpha = 0.6) +
        geom_hline(yintercept = 0, linetype = 2) +
        theme_classic() +
        theme(
          plot.title = element_text(face = 'italic'),
          legend.position = 'none'
        ) +
        labs(
          title = Sp,
          subtitle = paste0('Cold; condition: ', Cond),
          x = 'Mean annual temperature (°C)',
          y = expression('ln('~lambda~')')
        ) +
        xlim(temp_range_cold)
      )

      print(
      out_pars |>
        filter(species_id == Sp & border == 'hot' & sim == 2 & cond == Cond) |>
        select(iter, par, value) |>
        pivot_wider(names_from = par, values_from = value) |>
        mutate(
          chain = as.factor(case_when(
            iter <= 1000 ~ 1,
            iter <= 2000 ~ 2,
            iter <= 3000 ~ 3,
            iter > 3000 ~ 4
          )
        )) |>
        group_by(iter) |>
        expand_grid(
          temp = seq(temp_range_hot[1], temp_range_hot[2], length.out = 50)
        ) |>
        mutate(
          lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
        ) |>
        ggplot() +
        aes(temp, lambda) +
        facet_wrap(~chain) +
        geom_point(
          data = db_sim_hot |>
            slice_sample(n = 1e4),
          aes(temp, log(lambda)),
          alpha = 0.4, size = 0.3
        ) +
        stat_lineribbon(.width = 0.9, alpha = 0.6, color = 'purple') +
        geom_hline(yintercept = 0, linetype = 2) +
        theme_classic() +
        theme(legend.position = 'none') +
        labs(
          title = Sp,
          subtitle = paste0('Hot; condition: ', Cond),
          x = 'Mean annual temperature (°C)',
          y = NULL
        ) +
        xlim(temp_range_hot)
      )
    }})
  }
  dev.off()

#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Raw predicted distribution for each species and cond (sim == 2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  pdf(file.path(sim_path, 'pred.pdf'), height = 6, width = 9)
  for(Sp in unique(sim_pars$species_id))
  {
    try({

    # get which db_origin is more abundant
    db_sim_cold = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', 1, '.RDS')) |>
      bind_rows(
        readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', 2, '.RDS'))
      )

    db_sim_hot = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', 1, '.RDS')) |>
      bind_rows(
        readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', 2, '.RDS'))
      )

    # get temperature range for X axis limits
    db_sim_cold |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_cold

    db_sim_hot |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_hot

    # define Y axis limits
    db_sim_hot |>
      bind_cols(border = 'hot') |>
      bind_rows(
        db_sim_cold |>
          bind_cols(border = 'cold')
      ) |>
      pull(lambda) |>
      log() |>
      quantile(probs = c(0.001, 0.999)) ->
    yLim
    
    out_pars |>
      filter(species_id == Sp & border == 'cold' & sim == 2 & cond != 3) |>
      select(cond, iter, par, value) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond, iter) |>
      expand_grid(
        temp = seq(temp_range_cold[1], temp_range_cold[2], length.out = 50)
      ) |>
      mutate(
        lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
      ) |>
      ggplot() +
      aes(temp, lambda) +
      facet_wrap(
        ~cond,
        ncol = 1,
        labeller = as_labeller(setNames(c('No competition', 'Heterospecific competition'), 1:2))
      ) +
      geom_point(
        data = db_sim_cold |>
          group_by(cond) |>
          slice_sample(n = 1e4),
        aes(temp, log(lambda)),
        alpha = 0.4, size = 0.3
      ) +
      stat_lineribbon(.width = 0.9, alpha = 0.6, color = 'purple') +
      geom_hline(yintercept = 0, linetype = 2) +
      theme_classic() +
      theme(
        plot.title = element_text(face = 'italic'),
        legend.position = 'none'
      ) +
      labs(
        title = spIds |> filter(species_id == Sp) |> pull(species_name),
        subtitle = 'Cold',
        x = 'Mean annual temperature (°C)',
        y = expression('ln('~lambda~')')
      ) +
      xlim(temp_range_cold) +
      ylim(yLim) ->
    p1

    out_pars |>
      filter(species_id == Sp & border == 'hot' & sim == 2 & cond != 3) |>
      select(cond, iter, par, value) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond, iter) |>
      expand_grid(
        temp = seq(temp_range_hot[1], temp_range_hot[2], length.out = 50)
      ) |>
      mutate(
        lambda = sn::rsn(n(), beta_mean * temp + inter, exp(sigma_inter + beta_sigma * temp), alpha)
      ) |>
      ggplot() +
      aes(temp, lambda) +
      facet_wrap(
        ~cond,
        ncol = 1,
        labeller = as_labeller(setNames(c('No competition', 'Heterospecific competition'), 1:2))
      ) +
      geom_point(
        data = db_sim_hot |>
          group_by(cond) |>
          slice_sample(n = 1e4),
        aes(temp, log(lambda)),
        alpha = 0.4, size = 0.3
      ) +
      stat_lineribbon(.width = 0.9, alpha = 0.6, color = 'purple') +
      geom_hline(yintercept = 0, linetype = 2) +
      theme_classic() +
      theme(legend.position = 'none') +
      labs(
        title = '',
        subtitle = 'Hot',
        x = 'Mean annual temperature (°C)',
        y = NULL
      ) +
      xlim(temp_range_hot) +
      ylim(yLim) ->
    p2

    print(ggarrange(p1, p2, ncol = 2))
  })
  }
  dev.off()

#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Raw predicted distribution for each species and cond (sim == 2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  get_ext <- function(mean_y, var_y, alpha) {
    P = ecdf(sn::rsn(1000, mean_y, var_y, alpha))
    return(1 - P(0))
  }


  pdf(file.path(sim_path, 'extinction_risk.pdf'), height = 3.5, width = 9)
  for(Sp in unique(sim_pars$species_id))
  {
    # get which db_origin is more abundant
    db_sim_cold = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', 1, '.RDS'))
    db_sim_hot = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', 1, '.RDS'))

    # get range max 
    db_sim_cold |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_cold
    db_sim_hot |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_hot

    # get species range from db
    treeData |>
      filter(species_id == Sp) |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    species_range
    
    out_pars |>
      filter(species_id == Sp & border == 'cold' & sim == 2 & cond != 3) |>
      group_by(cond, par) |>
      reframe(value = mean(value)) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond) |>
      expand_grid(
        temp = seq(species_range[1] - abs(species_range[1] * 0.2), temp_range_cold[2], length.out = 50)
      ) |>
      mutate(
        mean = beta_mean * temp + inter,
        sd = exp(sigma_inter + beta_sigma * temp)
      ) |>
      rowwise() |>
      mutate(
        ext = get_ext(mean, sd, alpha)
      ) |>
      ggplot() +
      aes(temp, ext) +
      aes(color = cond) +
      geom_line(linewidth = 1.3, alpha = 0.8) +
      geom_vline(xintercept = species_range[1], linetype = 2, alpha = 0.6) +
      theme_classic() +
      scale_color_manual(values = c('#5ab4ac', '#d8b365')) +
      theme(
        plot.title = element_text(face = 'italic'),
        legend.position = 'none'
      ) +
      labs(
        title = spIds |> filter(species_id == Sp) |> pull(species_name),
        subtitle = 'Cold',
        x = 'Mean annual temperature (°C)',
        y = 'Suitable probability'
      ) +
      ylim(0, 1) ->
    p1

    out_pars |>
      filter(species_id == Sp & border == 'hot' & sim == 2 & cond != 3) |>
      group_by(cond, par) |>
      reframe(value = mean(value)) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond) |>
      expand_grid(
        temp = seq(temp_range_hot[1], species_range[2] + species_range[2] * 0.2, length.out = 50)
      ) |>
      mutate(
        mean_y = beta_mean * temp + inter,
        sd_y = exp(sigma_inter + beta_sigma * temp)
      ) |>
      rowwise() |>
      mutate(
        ext = get_ext(mean_y, sd_y, alpha)
      ) |>
      ggplot() +
      aes(temp, ext) +
      aes(color = cond) +
      geom_line(linewidth = 1.3, alpha = 0.8) +
      geom_vline(xintercept = species_range[2], linetype = 2, alpha = 0.6) +
      theme_classic() +
      scale_color_manual(values = c('#5ab4ac', '#d8b365')) +
      theme(
        axis.text.y = element_blank(),
        legend.position = 'none'
      ) +
      labs(
        title = '',
        subtitle = 'Hot',
        x = 'Mean annual temperature (°C)',
        y = NULL
      ) +
      ylim(0, 1) ->
    p2

    print(ggarrange(p1, p2, ncol = 2))
  }
  dev.off()

#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate suitable probability at both center and range limit (cold + hot)
# So I can compute the difference (reduce, 0, or increase?) in suitable prob
# between center vs border
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  sp_dt <- tibble()
  for(Sp in unique(sim_pars$species_id))
  {

    # get which db_origin is more abundant
    db_sim_cold = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_cold_', 2, '_', 1, '.RDS'))
    db_sim_hot = readRDS(paste0('simulations/model_lambdaPlot/out_species/', Sp, '_hot_', 2, '_', 1, '.RDS'))

    # get range max 
    db_sim_cold |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_cold
    db_sim_hot |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    temp_range_hot

    # get species range from db
    treeData |>
      filter(species_id == Sp) |>
      pull(bio_01_mean) |>
      range(na.rm = TRUE) ->
    species_range
    
    out_pars |>
      filter(species_id == Sp & border == 'cold' & sim == 2 & cond != 3) |>
      select(iter, par, value, cond) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond, iter) |>
      expand_grid(
        temp = c(species_range[1], temp_range_cold[2])
      ) |>
      mutate(
        mean_y = beta_mean * temp + inter,
        sd_y = exp(sigma_inter + beta_sigma * temp)
      ) |>
      rowwise() |>
      mutate(
        ext = get_ext(mean_y, sd_y, alpha)
      ) |>
      mutate(
        range_pos = if_else(temp == species_range[1], 'border', 'center')
      ) |>
      bind_cols(
        border = 'Cold'
      ) ->
    db_cold

    out_pars |>
      filter(species_id == Sp & border == 'hot' & sim == 2 & cond != 3) |>
      select(iter, par, value, cond) |>
      pivot_wider(names_from = par, values_from = value) |>
      group_by(cond, iter) |>
      expand_grid(
        temp = c(temp_range_hot[1], species_range[2])
      ) |>
      mutate(
        mean_y = beta_mean * temp + inter,
        sd_y = exp(sigma_inter + beta_sigma * temp)
      ) |>
      rowwise() |>
      mutate(
        ext = get_ext(mean_y, sd_y, alpha)
      ) |>
      mutate(
        range_pos = if_else(temp == species_range[2], 'border', 'center')
      ) |>
      bind_cols(
        border = 'Hot'
      ) ->
    db_hot

    sp_dt <- rbind(
      sp_dt,
      db_cold |>
        bind_rows(db_hot) |>
        bind_cols(species_id = Sp)
    )
  }

  saveRDS(sp_dt, file.path(sim_path, 'suitable_prob.RDS'))
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figures for suitable probability at the border and center of
# hot and cold range positions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sp_dt <- readRDS(file.path(sim_path, 'suitable_prob.RDS'))

sp_dt |>
  mutate(
    cond = case_match(
      as.character(cond),
      '1' ~ 'Fundamental niche',
      '2' ~ 'Realized niche'
    )
  ) |>
  ggplot() +
  aes(temp, ext) +
  aes(color = range_pos, fill = range_pos) +
  facet_grid(cond~border, scales = 'free_x') +
  geom_hdr(
    aes(color = NULL),
    probs = .9, alpha = .4,
    method = 'mvnorm'
  ) +
  ylim(0, 1.1) +
  stat_pointinterval() +
  theme_classic() +
  theme(legend.position = 'top') +
  scale_color_manual(values = c('#fc8d59', '#91bfdb')) +
  scale_fill_manual(values = c('#fc8d59', '#91bfdb')) +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'Suitable probability',
    color = NULL,
    fill = NULL
  )

sp_dt |>
  filter(range_pos == 'border') |>
  mutate(
    cond = case_match(
      as.character(cond),
      '1' ~ 'Fundamental niche',
      '2' ~ 'Realized niche'
    )
  ) |>
  ggplot() +
  aes(temp, ext) +
  aes(color = cond, fill = cond) +
  facet_wrap(~border, scales = 'free_x') +
  geom_hdr(
    aes(color = NULL),
    probs = .9, alpha = .4,
    method = 'mvnorm'
  ) +
  ylim(0, 1.1) +
  stat_pointinterval() +
  theme_classic() +
  theme(legend.position = 'top') +
  scale_color_manual(values = c('#5ab4ac', '#d8b365')) +
  scale_fill_manual(values = c('#5ab4ac', '#d8b365')) +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'Suitable probability',
    color = NULL,
    fill = NULL
  )

sp_dt |>
  filter(range_pos == 'border') |>
  select(species_id, temp, iter, cond, range_pos, border, ext) |>
  pivot_wider(
    names_from = cond,
    values_from = ext,
    names_prefix = 'cond_'
  ) |>
  mutate(
    sp_diff = cond_2-cond_1
  ) |>
  # filter(species_id != '19462FAGGRA') |>
  ggplot() +
  aes(temp, sp_diff) +
  aes(color = border, fill = border) +
  facet_wrap(~border, scales = 'free_x') +
  geom_hdr(
    aes(color = NULL),
    probs = .75, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.7) +
  stat_pointinterval() +
  theme_classic() +
  theme(legend.position = 'top') +
  scale_color_manual(values = c('#91bfdb', '#fc8d59')) +
  scale_fill_manual(values = c('#91bfdb', '#fc8d59')) +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'Suitable probability difference (realized - fundamental)',
    color = NULL,
    fill = NULL
  )


sp_dt |>
  filter(range_pos == 'center') |>
  select(species_id, temp, iter, cond, range_pos, border, ext) |>
  pivot_wider(
    names_from = cond,
    values_from = ext,
    names_prefix = 'cond_'
  ) |>
  mutate(
    sp_diff = cond_2-cond_1
  ) |>
  # filter(species_id != '19462FAGGRA') |>
  ggplot() +
  aes(temp, sp_diff) +
  aes(color = border, fill = border) +
  facet_wrap(~border, scales = 'free_x') +
  geom_hdr(
    aes(color = NULL),
    probs = .75, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.7) +
  stat_pointinterval() +
  theme_classic() +
  theme(legend.position = 'top') +
  scale_color_manual(values = c('#91bfdb', '#fc8d59')) +
  scale_fill_manual(values = c('#91bfdb', '#fc8d59')) +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'Suitable probability difference (realized - fundamental)',
    color = NULL,
    fill = NULL
  )


sp_dt |>
  select(species_id, iter, cond, range_pos, border, ext) |>
  pivot_wider(
    names_from = cond,
    values_from = ext,
    names_prefix = 'cond_'
  ) |>
  mutate(
    sp_diff = cond_2-cond_1
  ) |>
  group_by(species_id, range_pos, border) |>
  reframe(
    cond_1 = mean(cond_1),
    cond_2 = mean(cond_2)
  ) |>
  ggplot() +
  aes(cond_1, cond_2) +
  aes(color = range_pos, fill = range_pos) +
  facet_wrap(~border) +
  geom_hdr(
    aes(color = NULL),
    probs = .75, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, alpha = 0.7) +
  scale_color_manual(values = c('#fc8d59', '#91bfdb')) +
  scale_fill_manual(values = c('#fc8d59', '#91bfdb')) +
  theme_classic() +
  labs(
    x = 'Suitable probability - Fundamental niche',
    y = 'Suitable probability - Realized niche',
    color = NULL,
    fill = NULL
  ) +
  theme(legend.position = 'top')

sp_dt |>
  select(species_id, iter, cond, range_pos, border, ext) |>
  mutate(
    cond = case_match(
      as.character(cond),
      '1' ~ 'Fundamental niche',
      '2' ~ 'Realized niche'
    )
  ) |>
  pivot_wider(
    names_from = range_pos,
    values_from = ext,
    names_prefix = 'cond_'
  ) |>
  group_by(species_id, border, cond) |>
  reframe(
    cond_center = mean(cond_center),
    cond_border = mean(cond_border)
  ) |>
  ggplot() +
  aes(cond_center, cond_border) +
  aes(color = cond) +
  # geom_point(size = 0.2, alpha = 0.5) +
  aes(fill = cond) +
  geom_hdr(
    aes(color = NULL),
    probs = .75, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_point() +
  facet_wrap(~border) +
  scale_color_manual(values = c('#5ab4ac', '#d8b365')) +
  scale_fill_manual(values = c('#5ab4ac', '#d8b365')) +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    y = 'Suitable probability at border',
    x = 'Suitable probability at center',
    color = NULL,
    fill = NULL
  ) +
  xlim(0.05, 1) +
  ylim(0.05, 1) +
  theme(legend.position = 'top')
