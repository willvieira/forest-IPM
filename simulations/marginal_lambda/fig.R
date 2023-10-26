# Marginal lambda

library(tidyverse)
library(ggdist)
library(ggtext)
library(ggrepel)
library(ggiraph)


# load parameters and output
sim_path <- file.path('simulations', 'marginal_lambda')

lambdas <- readRDS(file.path(sim_path, 'output_complete.RDS'))


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
  )

# remove parameters value to remove file size so
# we can push to github
dir.create(file.path(sim_path, 'output_processed'))
lambdas |>
  select(species_id, clim, var, comp, replication, sim, lambda) |>
  saveRDS(file.path(sim_path, 'output_processed', 'lambdas.RDS'))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Competition effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambdas |>
  filter(sim == 'competition') |>
  # return comp value back to numeric as climate sim got out
  mutate(comp = as.numeric(comp)) |>
  mutate(
    lambda = log(lambda),
    var = case_match(
      var,
      'cons' ~ 'Conspecific',
      'het' ~ 'Heterospecific'
    )
  ) |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  group_by(species_name, var, comp) |>
  reframe(
    mean_lambda = mean(lambda),
    sd_lambda = sd(lambda),
    ci_9 = mean_lambda + qt( c(0.9), n() - 1) * sd_lambda,
    ci_1 = mean_lambda + qt( c(0.1), n() - 1) * sd_lambda
  ) |> 
  ggplot() +
  aes(comp, mean_lambda) +
  geom_ribbon_interactive(aes(ymin = ci_1, ymax = ci_9), alpha = 0.3,color=NA) +
  geom_line_interactive() +
  aes(tooltip = species_name, data_id = species_name) +
  facet_wrap(~var) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  theme_classic() +
  labs(
    x = 'Plot basal area (m2/ha)',
    y = expression('ln('~lambda~')')
  ) +
  theme(legend.position = 'none') ->
p1

tooltip_css <- "background-color:#fc8d59;padding:5px;border-radius:3px;font-style:italic;"

girafe(
  ggobj = p1,
  options = list(
    opts_tooltip(css = tooltip_css, opacity = 1),
    opts_sizing(width = .7),
    opts_hover_inv(css = "opacity:0.1;"),
    opts_hover(css = "stroke-width:2.5px;")
  )
)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Climate effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lambdas |>
  filter(species_id %in% spIds$species_id_old) |>
  filter(sim == 'climate') |>
  mutate(
    var = case_match(
      var,
      'temp' ~ 'Temperature',
      'prec' ~ 'Precipitation'
    ),
    comp = case_match(
      comp,
      'high' ~ 'High competition',
      'low' ~ 'Low competition'
    ),
    comp = factor(comp, levels = c('Low competition', 'High competition'))
  ) |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  group_by(species_name, var, comp, clim) |>
  reframe(
    mean_lambda = mean(lambda),
    sd_lambda = sd(lambda),
    ci_9 = mean_lambda + qt( c(0.9), n() - 1) * sd_lambda,
    ci_1 = mean_lambda + qt( c(0.1), n() - 1) * sd_lambda
  ) ->
clim_dt

clim_dt |>
  filter(var == 'Temperature') |>
  ggplot() +
  aes(clim, mean_lambda) +
  geom_ribbon_interactive(aes(ymin = ci_1, ymax = ci_9), alpha = 0.3,color=NA) +
  geom_line_interactive() +
  aes(tooltip = species_name, data_id = species_name) +
  facet_grid(~comp) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  theme_classic() +
  labs(
    x = 'Mean annual temperature (scaled)',
    y = expression('ln('~lambda~')')
  ) +
  theme(legend.position = 'none') ->
p_temp

clim_dt |>
  filter(var == 'Precipitation') |>
  ggplot() +
  aes(clim, mean_lambda) +
  geom_ribbon_interactive(aes(ymin = ci_1, ymax = ci_9), alpha = 0.3,color=NA) +
  geom_line_interactive() +
  aes(tooltip = species_name, data_id = species_name) +
  facet_grid(~comp) +
  geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  theme_classic() +
  labs(
    x = 'Mean annual precipitation (scaled)',
    y = expression('ln('~lambda~')')
  ) +
  theme(legend.position = 'none') ->
p_prec
