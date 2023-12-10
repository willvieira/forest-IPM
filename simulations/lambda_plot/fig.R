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
sim_path <- file.path('simulations', 'lambda_plot')

# parameters
sim_pars <- readRDS(file.path(sim_path, 'simulation_pars.RDS'))

# output
sim_out <- readRDS(file.path(sim_path, 'final_output.RDS'))

# data for latitude - longitude
treeData <- readRDS(
    paste0(readLines('_data.path'), 'treeData.RDS')
  ) |>
  filter(species_id %in% unique(sim_pars$species_id)) |>
  group_by(species_id, plot_id) |>
  reframe(
    latitude = unique(latitude),
    longitude = unique(longitude)
  ) ->
plot_position

# merge pars with output
sim_out |>
  left_join(
    sim_pars |>
      ungroup() |>
      mutate(array_id = row_number()) |>
      select(species_id, plot_id, year_measured, plot_size, array_id)
  ) |>
  mutate(sim = factor(sim)) |>
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
  filter(species_id_old %in% unique(out$species_id)) |>
  select(species_id_old, species_name) |>
  rename(species_id = species_id_old)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compute summires (mean, sd, extinction risk) across replications
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to get area under the curve from threshold
area_thr <- function(x, threshold = 1) {
  dd = density(x)
  den_weight =  dd$y
  sum( (den_weight/sum(den_weight))[dd$x < threshold] )
}

out |>
  mutate(lambda = log(lambda)) |>
  group_by(species_id, plot_id, year_measured, sim) |>
  mutate(
    mean_lambda = mean(lambda),
    sd_lambda = sd(lambda),
    sd_temp = mean(temp_sd),
    ext = area_thr(lambda, threshold = 0)
  ) |>
  slice_head(n = 1) ->
out_summ

# group plots as bellow or above the center of the distribution (50% quantile)
out_summ |>
  group_by(species_id) |>
  mutate(
    plot_pos = ifelse(temp < quantile(temp, 0.5), 'Cold range', 'Hot range')
  ) |>
  group_by(species_id, plot_id, year_measured, sim) ->
out_summ



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# How the average lambda changes from demographic only stochasticity
# to demo + environment (sim2) and demo + envir + BA (sim3)
# Delta average lambda from sim1 to sim2 and sim3

# get delta diff for average lambda, sd lambda, and ext
out_summ |>
  group_by(species_id, plot_id, year_measured) |>
  mutate(
    d_meanLambda_sim2 = mean_lambda[sim == 2] - mean_lambda[sim == 1],
    d_meanLambda_sim3 = mean_lambda[sim == 3] - mean_lambda[sim == 1],
    d_sdLambda_sim2 = sd_lambda[sim == 2] - sd_lambda[sim == 1],
    d_sdLambda_sim3 = sd_lambda[sim == 3] - sd_lambda[sim == 1],
    d_ext_sim2 = ext[sim == 2] - ext[sim == 1],
    d_ext_sim3 = ext[sim == 3] - ext[sim == 1]
  ) |>
  slice_head(n = 1) ->
out_delta

# mean lambda
out_delta |>
  select(species_id, latitude, contains('d_meanLambda')) |>
  pivot_longer(contains('d_')) |>
  ggplot() +
  aes(value, fct_reorder(species_id, value)) +
  facet_wrap(~name, scales = 'free_x') +
  stat_pointinterval()

# sd lambda
out_delta |>
  select(species_id, latitude, contains('d_sdLambda')) |>
  pivot_longer(contains('d_')) |>
  ggplot() +
  aes(value, fct_reorder(species_id, value)) +
  facet_wrap(~name, scales = 'free_x') +
  stat_pointinterval()

# extinction risk
out_delta |>
  select(species_id, latitude, contains('d_ext')) |>
  pivot_longer(contains('d_')) |>
  ggplot() +
  aes(value, fct_reorder(species_id, value)) +
  facet_wrap(~name, scales = 'free_x') +
  stat_pointinterval()


# extinction in function of latitude
out_summ |>
  ggplot() +
  aes(latitude, ext, color = factor(sim), group = interaction(sim, plot_pos)) +
  geom_smooth() +
  facet_wrap(~species_id, scales = 'free')


# Temperature variability vs lambda variability
out_summ |>
  # clim stochasticity only
  filter(sim == 2) |>
  ggplot() +
  aes(temp_sd, sd_lambda) +
  facet_wrap(~species_id) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(level = NA) +
  theme(legend.position = 'none')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GLM with quasibinomial logit link between latitude and extinction risk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# slope for each species-model combination
out_summ |>
  group_by(species_id, plot_pos, sim) |>
  select(latitude, ext) |>
  mutate(
    ext = case_when(
      ext > 1 ~ 1,
      .default = ext
    )
  ) |>
  nest() |>
  mutate(
    model_out = map(data, ~glm(data = ., ext~latitude, family = 'quasibinomial')),
    coeffs = map(model_out, coefficients)
  ) |>
  select(!c(data, model_out)) |>
  unnest(coeffs) |>
  group_by(species_id, plot_pos, sim) |>
  mutate(par = c('int', 'slope')) |>
  filter(par == 'slope') |>
  left_join(spIds) |>
  ggplot() +
  aes(coeffs, fct_reorder(species_name, coeffs)) +
  aes(color = factor(sim)) +
  facet_wrap(~plot_pos) +
  stat_pointinterval() +
  theme_classic() +
  geom_vline(xintercept = 0, alpha = 0.6) +
  labs(
    x = 'Slope',
    y = NULL,
    color = 'Model'
  ) +
  theme(axis.text.y = element_text(face = 'italic'))


# slope for each species-model combination
out_summ |>
  group_by(species_id, plot_pos, sim) |>
  select(latitude, ext) |>
  mutate(
    ext = case_when(
      ext > 1 ~ 1,
      .default = ext
    )
  ) |>
  nest() |>
  mutate(
    model_out = map(data, ~glm(data = ., ext~latitude, family = 'quasibinomial')),
    coeffs = map(model_out, coefficients)
  ) |>
  select(!c(data, model_out)) |>
  unnest(coeffs) |>
  group_by(species_id, plot_pos, sim) |>
  mutate(par = c('int', 'slope')) |>
  filter(par == 'slope') |>
  left_join(spIds) |>
  mutate(
    sim = case_match(
      as.numeric(sim),
      1 ~ 'Demog',
      2 ~ 'Demog + envir',
      3 ~ 'Demog + envir + comp'
    )
  ) |>
  ggplot() +
  aes(coeffs, fct_reorder(species_name, coeffs)) +
  aes(color = plot_pos) +
  aes(shape = sim) +
  geom_jitter(size = 1.8) +
  theme_classic() +
  geom_vline(xintercept = 0, alpha = 0.6) +
  scale_color_manual(values = c('#91bfdb', '#fc8d59')) +
  labs(
    x = 'Slope',
    y = NULL,
    color = 'Plots range position',
    shape = 'Model'
  ) +
  theme(
    axis.text.y = element_text(face = 'italic'),
    panel.grid.major.y = element_line(colour = rgb(0,0,0,.1))
  )






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lambda, lambda sd, and extinction risk in function of mean annual temperature
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(
  file = file.path(sim_path, 'lambda_mean.pdf'),
  width = 14, height = 12
)

vars_to_plot <- c('mean_lambda', 'sd_lambda', 'ext')
vars_name <- c(expression(bar(lambda)), expression(sigma[lambda]), 'Extinction risk')

for(var in 1:length(vars_to_plot))
{
  print(
  out_summ |>
  # filter(species_id %in% '18032ABIBAL') |>
    left_join(spIds) |>
    ggplot() +
    aes(temp, !! rlang::sym(vars_to_plot[var])) +
    aes(color = sim) +
    facet_wrap(~species_name, scales = 'free') +
    geom_smooth() +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.8) +
    theme_classic() +
    labs(
      x = 'Mean annual temperature (°C)',
      y = vars_name[var],
      color = 'Simulation'
    ) +
    theme(
      legend.position = 'top',
      axis.text.y = element_text(face = 'italic')
    )
  )
}

for(var in 1:length(vars_to_plot))
{
  print(
  out_summ |>
  # filter(species_id %in% '18032ABIBAL') |>
    left_join(spIds) |>
    ggplot() +
    aes(temp, !! rlang::sym(vars_to_plot[var])) +
    aes(color = sim, group=interaction(sim, plot_pos)) +
    facet_wrap(~species_name, scales = 'free_y') +
    geom_smooth(method = 'lm') +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.8) +
    theme_classic() +
    labs(
      x = 'Mean annual temperature (°C)',
      y = vars_name[var],
      color = 'Simulation'
    ) +
    theme(
      legend.position = 'top',
      axis.text.y = element_text(face = 'italic')
    )
  )
}

# temperature variance
out_summ |>
  filter(sim == 2) |>
  left_join(spIds) |>
  ggplot() +
  aes(temp, sd_temp) +
  aes(group = plot_pos) +
  facet_wrap(~species_name) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  labs(
    x = 'Mean annual temperature (°C)',
    y = expression(sigma[temp])
  ) +
  theme(
    legend.position = 'top',
    axis.text.y = element_text(face = 'italic')
  )

dev.off()





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lambda, lambda sd, and extinction risk in function of mean annual temperature
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# slope for each species-model combination
out_summ |>
  group_by(species_id, plot_pos, sim) |>
  select(temp, mean_lambda, sd_lambda, ext, sd_temp) |>
  mutate(
    ext = case_when(
      ext > 1 ~ 1,
      .default = ext
    )
  ) |>
  nest() |>
  mutate(
    model_out = map(data, ~glm(data = ., mean_lambda~temp, family = 'gaussian')),
    coeffs_lambda = map(model_out, coefficients),
    model_out = map(data, ~glm(data = ., sd_lambda~temp, family = 'gaussian')),
    coeffs_sd = map(model_out, coefficients),
    model_out = map(data, ~glm(data = ., ext~temp, family = 'quasibinomial')),
    coeffs_ext = map(model_out, coefficients),
    coeffs_sd = map(model_out, coefficients),
    model_out = map(data, ~glm(data = ., sd_temp~temp, family = 'gaussian')),
    coeffs_tempSD = map(model_out, coefficients) 
  ) |>
  select(!c(data, contains('model_'))) |>
  # get slope parameter
  mutate(across(contains('coeffs'), \(x) unlist(x)[2])) |>
  pivot_longer(cols = contains('coeffs')) |>
  mutate(
    value = if_else(plot_pos == 'Cold range', value * -1, value)
  ) |>
  pivot_wider() |>
  left_join(spIds) |>
  mutate(
    sim = factor(case_match(
      as.numeric(sim),
      1 ~ 'Parameter',
      2 ~ 'Param + environment',
      3 ~ 'Param + envir + competition'
    ))
  ) |>
  filter(sim %in% c('Parameter', 'Param + environment')) ->
plot_data


plot_data |>
  filter(sim == 'Parameter') |>
  ggplot() +
  aes(coeffs_lambda, coeffs_ext) +
  aes(color = sim, fill = sim) +
  facet_wrap(~plot_pos) +
  geom_hdr(
    aes(color = NULL),
    probs = .9, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, alpha = 0.3, linetype = 2) +
  theme_classic() +
  scale_color_manual(values = c('#91bfdb', '#fc8d59')) +
  scale_fill_manual(values = c('#91bfdb', '#fc8d59')) +
  labs(
    x = expression('Slope ('~bar(lambda)%~%~textstyle(Mean~annual~temperature)~')'),
    y = expression('Slope ('~textstyle(Extinction~risk)%~%~textstyle(Mean~annual~temperature)~')'),
    color = 'Uncertainty simulation',
    fill = 'Uncertainty simulation'
  ) +
  xlim(-0.075, 0.075) +
  ylim(-1, 1) + 
  annotate(geom = 'text', x = -0.06, y = 0.95, label = expression(lambda %down% E %up% '')) +
  annotate(geom = 'text', x = 0.06, y = 0.95, label = expression(lambda %up% E %up% '')) +
  annotate(geom = 'text', x = -0.06, y = -0.95, label = expression(lambda %down% E %down% '')) +
  annotate(geom = 'text', x = 0.06, y = -0.95, label = expression(lambda %up% E %down% '')) ->
p1

plot_data |>
  ggplot() +
  aes(coeffs_lambda, coeffs_ext) +
  aes(color = sim, fill = sim) +
  facet_wrap(~plot_pos) +
  geom_hdr(
    aes(color = NULL),
    probs = .9, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, alpha = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, alpha = 0.3, linetype = 2) +
  theme_classic() +
  scale_color_manual(values = c('#91bfdb', '#fc8d59', 'yellow')) +
  scale_fill_manual(values = c('#91bfdb', '#fc8d59', 'yellow')) +
  labs(
    x = expression('Slope ('~bar(lambda)%~%~textstyle(Mean~annual~temperature)~')'),
    y = expression('Slope ('~textstyle(Extinction~risk)%~%~textstyle(Mean~annual~temperature)~')'),
    color = 'Uncertainty simulation',
    fill = 'Uncertainty simulation'
  ) +
  xlim(-0.075, 0.075) +
  ylim(-1, 1) + 
  annotate(geom = 'text', x = -0.06, y = 0.95, label = expression(lambda %down% E %up% '')) +
  annotate(geom = 'text', x = 0.06, y = 0.95, label = expression(lambda %up% E %up% '')) +
  annotate(geom = 'text', x = -0.06, y = -0.95, label = expression(lambda %down% E %down% '')) +
  annotate(geom = 'text', x = 0.06, y = -0.95, label = expression(lambda %up% E %down% '')) ->
p2

p2 +
  geom_text_repel(aes(color = NULL, label = species_name), size = rel(1.8), fontface = 'italic') ->
p3

p2 +
  aes(color = coeffs_sd, shape = sim) +
  scale_colour_gradient2() +
  labs(
    fill = 'Uncertainty simulation',
    shape = 'Uncertainty simulation',
    color = expression(sigma[lambda]) 
  ) ->
p4


pdf(
  file = file.path(sim_path, 'lambda_slope2.pdf'),
  width = 10, height = 5
)
for(i in 1:4) print(get(paste0('p', i)))
dev.off()
