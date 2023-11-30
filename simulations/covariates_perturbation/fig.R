# Output preparation and viz

library(tidyverse)
library(ggdist)
library(ggtext)
library(ggrepel)
library(ggpubr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load stuff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load parameters and output
sim_path <- file.path('simulations', 'covariates_perturbation')

# parameters
sim_pars <- readRDS(file.path(sim_path, 'simulation_pars.RDS'))

# output
sim_out <- readRDS(file.path(sim_path, 'final_output.RDS'))

# data for latitude - longitude
treeData <- readRDS(paste0(readLines('_data.path'), 'treeData.RDS')) |>
  filter(species_id %in% unique(sim_pars$species_id))

treeData |> 
  group_by(species_id, plot_id) |>
  reframe(
    latitude = unique(latitude),
    longitude = unique(longitude)
  ) ->
plot_position

# merge pars with output
sim_out |>
  bind_rows() |>
  left_join(
    sim_pars |>
      ungroup() |>
      mutate(array_id = row_number()) |>
      select(species_id, plot_id, year_measured, plot_size, array_id)
  ) |>
  select(!array_id) |>
  left_join(plot_position) ->
out

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
  filter(species_id_old %in% unique(out$species_id))



# /!\ important /!\
# I made a mistake in the ipm_i.R code where I wrongly diveded the lambda
# difference by the difference of the covariate using the scaled version
# This happened for temp and prec covariates.
# I should instead divide by the difference in the covariate in the natural scal

out |>
  # FIX TEMP and PREC partial deriv
  mutate(
    par.temp = (par.temp * 0.01)/(temp_pertb - temp),
    par.prec = (par.prec * 0.01)/(prec_pertb - prec)
  ) |>
  # For some reason (I think it's because of inviction after the small pertubation in size), some plots had a reduction in the plot basal area instead of an increase
  filter(BA_het_pertb > BA_het) ->
out


# Average marginal effect for each species
out |>
  group_by(species_id, rep) |>
  reframe(across(contains('par.'), mean, na.rm = TRUE)) ->
AME_sp


out |>
  group_by(species_id, plot_id, year_measured) |>
  mutate(across(contains('par.'), \(x) mean(log(x)))) |>
  slice_head(n = 1) ->
out_summ

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Species level evalution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf('simulations/covariates_perturbation/out.pdf', width = 10, height = 7)

# AME across species
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AME_sp |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  pivot_longer(cols = contains('par.')) |>
  ggplot() +
  aes(log(value), fct_reorder(species_name, value)) +
  aes(color = name) +
  stat_pointinterval() +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = 'italic'),
    legend.position = 'none'
  ) +
  scale_color_manual(
      values = c('#7fc97f', '#ffff99', '#386cb0', '#fdc086'),
      labels = c(
        expression(plain('BA')[conspecific]),
        expression(plain('BA')[heterospecific]),
        'Precipitaiton',
        'Temperature'
      )
    ) +
  labs(
    x = 'ln(Sensitivity)',
    y = NULL,
    color = NULL
  ) ->
p1

AME_sp |>
  group_by(species_id) |>
  reframe(across(contains('par.'), mean)) |>
  pivot_longer(cols = contains('par.')) |>
  mutate(
    name = str_replace(name, 'par.', ''),
    name = case_match(
      name,
      'BA_con' ~ 'BA cons',
      'BA_het' ~ 'BA het',
      'temp' ~ 'Temperature',
      'prec' ~ 'Precipitation'
    )
  ) |>
  ggplot() +
  aes(fct_reorder(name, log(value), .desc = TRUE), log(value)) +
  aes(fill = name) +
  geom_boxplot() +
  theme_classic() +
  labs(
    x = '',
    y = 'ln(Sensitivity)'
  ) +
  scale_fill_manual(values = c('#7fc97f', '#ffff99', '#386cb0', '#fdc086')) +
  theme(legend.position = 'none') ->
p2

ggarrange(p1, p2, ncol = 2)


# AME by the position in range
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AME_sp |>
  select(!contains('BA')) |>
  pivot_longer(cols = contains('par.')) |>
  left_join(
    treeData |>
      group_by(species_id, plot_id) |>
      reframe(
        bio_01_mean = mean(bio_01_mean, na.rm = TRUE)
      ) |>
      group_by(species_id) |>
      reframe(
        range_pos = median(bio_01_mean, na.rm = TRUE)
      )
  ) |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  group_by(species_name, name, range_pos) |>
  mutate(
    name = case_match(
      name,
      'par.temp' ~ 'Temperature',
      'par.prec' ~ 'Precipitation'
    )
  ) |>
  reframe(
    AME = mean(value),
    range_pos = mean(range_pos)
  ) ->
mean_rangeValue

AME_sp |>
  select(!contains('BA')) |>
  pivot_longer(cols = contains('par.')) |>
  left_join(
    treeData |>
      group_by(species_id, plot_id) |>
      reframe(
        bio_01_mean = mean(bio_01_mean, na.rm = TRUE)
      ) |>
      group_by(species_id) |>
      reframe(
        range_pos = median(bio_01_mean, na.rm = TRUE)
      )
  ) |>
  mutate(
    name = case_match(
      name,
      'par.temp' ~ 'Temperature',
      'par.prec' ~ 'Precipitation'
    )
  ) |>
  ggplot() +
  aes(range_pos, log(value)) +
  facet_wrap(~name) +
  stat_pointinterval() +
  geom_smooth(method = 'lm') +
  geom_text_repel(
    data = mean_rangeValue,
    aes(x = range_pos, y = log(AME), label = species_name),
    alpha = 0.8,
    size = 2,
    fontface = 'italic'
  ) +
  theme_classic() +
  labs(
    x = 'Range position in Mean Annual Temperature (°C)',
    y = 'ln(Sensitivity)'
  )


# Comp and clim effect across the temperature range 
out_summ |>
  ungroup() |>
  mutate(
    clim = log(exp(par.temp) + exp(par.prec)),
    comp = log(exp(par.BA_con) + exp(par.BA_het)),
    ccr = comp - clim
  ) |>
  select(species_id, temp, clim, comp, ccr) |>
  pivot_longer(cols = c(clim, comp, ccr)) |>
  mutate(
    name = factor(case_match(
      name,
      'comp' ~ 'Competition',
      'clim' ~ 'Climate',
      'ccr' ~ 'Competition/climate ratio'
    )),
    name = factor(name, levels = levels(name)[c(2, 1, 3)])
  ) |>
  left_join(spIds, by = c('species_id' = 'species_id_old')) |>
  ggplot() +
  aes(temp, value) +
  aes(color = name) +
  geom_smooth() +
  facet_wrap(~species_name) +
  theme_classic() +
  scale_color_manual(values = c('#d8b365', '#5ab4ac', 'black')) +
  geom_hline(yintercept = 0, alpha = 0.6) +
  theme(
    strip.text = element_text(face = "italic"),
    strip.background = element_blank(),
    legend.position = 'top'
  ) +
  labs(
    x = 'Mean annual temperature (°C)',
    y = NULL,
    color = NULL
  )



# function to determine if a plot is at lower, center, or upper border
# given how much of the tail we consider border (prob arg) [0-0.5]
which_border <- function(temp, prob = 0.1, naRm = TRUE) {
  temp_range = quantile(temp, probs = c(prob, 1 - prob), na.rm = naRm)
  # output vector with 'center' class
  out_pos = rep('Center', length(temp))
  # lower border
  out_pos[temp < temp_range[1]] = 'Cold'
  out_pos[temp > temp_range[2]] = 'Hot'
  
  return(out_pos)
}

# class of range position (cold, center, hot)
treeData |>
  group_by(species_id, plot_id) |>
  # get a single obs per plot to remove abundance effect
  reframe(bio_01_mean = mean(bio_01_mean)) |>
  group_by(species_id) |>
  mutate(
    border_cl = which_border(bio_01_mean, prob = 0.1)
  ) |>
  select(species_id, plot_id, border_cl) |>
  mutate(border_cl = factor(border_cl, levels = c('Cold', 'Center', 'Hot'))) ->
plotBorder_class


out_summ |>
  left_join(plotBorder_class) |>
  pivot_longer(cols = contains('par.')) |>
  group_by(species_id, name, border_cl) |>
  reframe(ame_ln = mean(value)) |>
  mutate(
    name = str_replace(name, 'par.', ''),
    name = factor(case_match(
      name,
      'BA_con' ~ 'BA conspecific',
      'BA_het' ~ 'BA heterospecific',
      'temp' ~ 'Temperature',
      'prec' ~ 'Precipitation'
    )),
    name = factor(name, levels = levels(name)[c(4, 1, 2, 3)])
  ) |>
  ggplot() +
  aes(border_cl, ame_ln) +
  aes(fill = border_cl) +
  facet_grid(~name) +
  geom_boxplot() +
  scale_fill_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
  theme_classic() +
  labs(
    x = NULL,
    y = 'ln(Sensitivity)',
    fill = NULL
  ) +
  theme(legend.position = 'top')

out_summ |>
  mutate(
    clim = log(exp(par.temp) + exp(par.prec)),
    comp = log(exp(par.BA_con) + exp(par.BA_het)),
    ccr = comp - clim
  ) |>
  left_join(plotBorder_class) |>
  group_by(species_id, border_cl) |>
  reframe(ccr = mean(ccr)) |>
  ggplot() +
  aes(border_cl, ccr) +
  aes(fill = border_cl) +
  geom_boxplot() +
  scale_fill_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
  theme_classic() +
  labs(
    x = NULL,
    y = 'Competition/climate ratio',
    fill = NULL
  )


treeData |>
      left_join(plotBorder_class) |>
      group_by(species_id, border_cl) |>
      reframe(
        range_pos = median(bio_01_mean, na.rm = TRUE)
      ) |>
      # add small noise to range position because if two species have the exact
      # same mean, ggdist will interpret as the a single distribution (bug?!)
      mutate(range_pos = range_pos + rnorm(n(), 0, 0.001)) ->
sp_range_pos


out_summ |>
  left_join(plotBorder_class) |>
  mutate(par.BA_het = ifelse(is.na(par.BA_het), 0, par.BA_het)) |>
  mutate(
    clim = log(exp(par.temp) + exp(par.prec)),
    comp = log(exp(par.BA_con) + exp(par.BA_het)),
    ccr = comp - clim
  ) |>
  left_join(spIds, by = c('species_id' = 'species_id_old')) |>
  ggplot() +
  aes(border_cl, ccr) +
  aes(fill = border_cl) +
  facet_wrap(~species_name) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  geom_boxplot() +
  scale_fill_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
  theme_classic() +
  theme(strip.text = element_text(face = "italic")) +
  labs(
    x = NULL,
    y = 'Competion/climate ratio',
    fill = NULL
  ) +
  ylim(-4, 2)


# relative effect of comp and clim across the range
out |>
  # three species that had the lowerst AME (low sensitivity to covariates)
  # filter(!species_id %in% c('NAQUEPRI', '19290QUEALB', '27821NYSSYL')) |>
  mutate(
    clim = log(par.temp + par.prec),
    comp = log(par.BA_con + par.BA_het),
    ccr = comp - clim
  ) |>
  left_join(plotBorder_class) |>
  group_by(species_id, border_cl, rep) |>
  mutate(ccr = mean(ccr)) |>
  slice_head(n = 1) |>
  left_join(sp_range_pos) |>
  filter(border_cl != 'Center') |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  ggplot() +
  aes(range_pos, ccr) +
  aes(color = border_cl) +
  stat_pointinterval(alpha = 0.7) +
  geom_quantile(quantiles = 0.5, linewidth = 1.2) +
  geom_quantile(quantiles = c(0.25, 0.75), alpha = 0.8, linewidth = .8) +
  geom_hline(yintercept = 0, alpha = 0.5, linetype = 2) +
  scale_color_manual(values = c('#91bfdb', '#fc8d59')) +
  theme_classic() +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'Competition/climate ratio',
    color = ''
  )


out |>
  # three species that had the lowerst AME (low sensitivity to covariates)
  # filter(!species_id %in% c('NAQUEPRI', '19290QUEALB', '27821NYSSYL')) |>
  mutate(
    clim = log(par.temp + par.prec),
    comp = log(par.BA_con + par.BA_het)
  ) |>
  left_join(plotBorder_class) |>
  group_by(species_id, border_cl, rep) |>
  mutate(
    Climate = mean(clim),
    Competition = mean(comp)
  ) |>
  slice_head(n = 1) |>
  pivot_longer(cols = c(Climate, Competition)) |>
  left_join(    
    treeData |>
      left_join(plotBorder_class) |>
      group_by(species_id, border_cl) |>
      reframe(
        range_pos = median(bio_01_mean, na.rm = TRUE)
      ) |>
      mutate(range_pos = range_pos + rnorm(n(), 0, 0.001))
  ) |>
  # filter(border_cl != 'Center') |>
  ggplot() +
  aes(range_pos, value) +
  aes(fill = name, color = name) +
  facet_wrap(~border_cl) +
  stat_pointinterval(alpha = 0.7) +
  ggdensity::geom_hdr(
    aes(color = NULL),
    probs = .75, alpha = .4,
    method = 'mvnorm'
  ) +
  # geom_smooth(method = 'lm', alpha = .8) +
  scale_color_manual(values = c('#d8b365', '#5ab4ac')) +
  scale_fill_manual(values = c('#d8b365', '#5ab4ac')) +
  theme_classic() +
  labs(
    x = 'Mean annual temperature (°C)',
    y = 'ln(Sensitivity)',
    color = '',
    fill = ''
  ) +
  theme(legend.position = 'top')

library(sf)

plot_spatial_sens <- function(sp)
{
  out_summ |>
    filter(species_id == sp) |>
    mutate(
      clim_ln = log(exp(par.temp) + exp(par.prec)),
      comp_ln = log(exp(par.BA_con) + exp(par.BA_het)),
      ccr_ln = comp_ln - clim_ln
    ) |>
    select(contains('tude'), contains('ln')) |>
    ungroup() ->
  sp_dt

  theme_plot <- list(
    theme_classic() +
    theme(
      legend.position = 'bottom',
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
			legend.key.height = unit(3, 'mm')
    )
  )

  sp_dt |>
    filter(
      comp_ln > quantile(comp_ln, 0.05) & comp_ln < quantile(comp_ln, 0.95)
    ) |>
    st_as_sf(
      coords = c('longitude', 'latitude'),
      crs = 4326
    ) |>
    ggplot() +
    geom_sf(aes(color = comp_ln), alpha = 0.65, size = 0.1) +
    scale_color_viridis_c() +
    labs(color = expression('ln('~S[competition]~')')) +
    theme_plot ->
  p1

  sp_dt |>
    filter(
      clim_ln > quantile(clim_ln, 0.05) & clim_ln < quantile(clim_ln, 0.95)
    ) |>
    st_as_sf(
      coords = c('longitude', 'latitude'),
      crs = 4326
    ) |>
    ggplot() +
    geom_sf(aes(color = clim_ln), alpha = 0.65, size = 0.1) +
    scale_color_viridis_c() +
    labs(color = expression('ln('~S[climate]~')')) +
    theme_plot ->
  p2


  sp_dt |>
    filter(
      ccr_ln > quantile(ccr_ln, 0.05) & ccr_ln < quantile(ccr_ln, 0.95)
    ) |>
    st_as_sf(
      coords = c('longitude', 'latitude'),
      crs = 4326
    ) |>
    ggplot() +
    geom_sf(aes(color = ccr_ln), alpha = 0.65, size = 0.1) +
    scale_color_gradient2() +
    labs(color = expression('ln(CCR)')) +
    theme_plot ->
  p3

  print(
      annotate_figure(
        ggarrange(p1, p2, p3, nrow = 1),
        top = text_grob(
          spIds |> filter(species_id_old == sp) |> pull(species_name),
          face = 'italic', size = 10
        )
      )
  )
}

plot_spatial_sens('28728ACERUB')
