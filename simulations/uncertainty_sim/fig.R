# Parameters vs output (sensitivity analysis)

library(tidyverse)
library(ggdist)
library(ggtext)
library(ggrepel)
library(ggpubr)


# load parameters and output
sim_path <- file.path('simulations', 'uncertainty_sim')

# parameters
readRDS(file.path(sim_path, 'output_complete.RDS')) |>
  # remove parameters
  map(~ .x |> select(!contains('.'))) ->
sim_out

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
  filter(species_id_old %in% unique(sim_out[[1]]$species_id))

# make the only parameters simulations be the "non variance" reference
sim_out[['temp']] |>
  bind_rows(
    sim_out[['pars']] |>
      bind_cols(
        temp_mean = NA,
        temp_sd = 0
      )
  ) |>
  left_join(
    spIds |> select(species_id_old, species_name),
    by = c('species_id' = 'species_id_old')
  ) ->
temp_sim

sim_out[['BA']] |>
  bind_rows(
    sim_out[['pars']] |>
      bind_cols(
        BA_mean = NA,
        BA_sd = 0,
        BA = NA,
        temp_mean = NA,
        temp_sd = 0
      )
  ) |>
  left_join(
    spIds |> select(species_id_old, species_name),
    by = c('species_id' = 'species_id_old')
  ) ->
BA_sim


sim_out[['BAtemp']] |>
  bind_rows(
    sim_out[['pars']] |>
      bind_cols(
        BA_mean = NA,
        BA_sd = 0,
        BA = NA,
        temp_mean = NA,
        temp_sd = 0
      )
  ) |>
  left_join(
    spIds |> select(species_id_old, species_name),
    by = c('species_id' = 'species_id_old')
  ) ->
BAtemp_sim


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some exploraratory figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

temp_sim |>
  mutate(noise_size = factor(temp_sd, labels = 'noise')) |>
  select(species_name, noise_size, lambda) |>
  bind_cols(sim = 'climate') |>
  bind_rows(
    BA_sim |>
      mutate(noise_size = factor(BA_sd, labels = 'noise')) |>
      select(species_name, noise_size, lambda) |>
      bind_cols(sim = 'competition')
  ) |>
  bind_rows(
    BAtemp_sim |>
      mutate(noise_size = factor(BA_sd, labels = 'noise')) |>
      select(species_name, noise_size, lambda) |>
      bind_cols(sim = 'climate + competition')
  ) |>
  mutate(
    sim = paste0('paramater + ', sim),
    sim = factor(sim, levels = c('paramater + climate', 'paramater + competition', 'paramater + climate + competition'))
  ) ->
lambda_dt

lambda_dt |>
  filter(lambda > quantile(lambda, 0.001)) |>
  ggplot() +
  aes(noise_size, log(lambda)) +
  aes(fill = species_name) +
  facet_wrap(~sim) +
  stat_slab(alpha = 0.6) +
  geom_hline(yintercept = 0, alpha = 0.4, linetype = 2) +
  scale_fill_manual(
    values = c('#7fc97f', '#beaed4', '#fdc086', '#e78ac3', '#386cb0')
  ) +
  theme_classic() +
  labs(
    y = expression('ln('~lambda~')'),
    x = expression('Size of noise ('~sigma~')'),
    fill = ''
  ) +
  theme(
    legend.text = element_text(face = 'italic'),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top'
  )

lambda_dt |>
  group_by(species_name, sim, noise_size) |>
  reframe(lambda_sd = sd(lambda)) |>
  ggplot() +
  aes(noise_size, lambda_sd) +
  aes(group = sim, color = sim) +
  # geom_path(linewidth = .8, alpha = 0.8) +
  geom_smooth(level = NA) +
  facet_wrap(~species_name, scales = 'free') +
  scale_color_manual(
    values = c('#7fc97f', '#fb8072', '#fdc086')
  ) +
  theme_classic() +
  labs(
    y = expression(sigma[lambda]),
    x = expression('Size of noise ('~sigma~')'),
    color = NULL
  ) +
  theme(
    strip.text= element_text(face = 'italic'),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top'
  )


# function to get area under the curve from threshold
area_thr <- function(x, threshold = 1, normalized = TRUE) {
  dd = density(x)
  den_weight =  dd$y
  if(normalized) {
    return( sum( (den_weight/sum(den_weight))[dd$x < threshold] ) )
  }else{
    return( sum( den_weight[dd$x < threshold] ) )
  }
}

lambda_dt |>
  filter(species_name != 'Picea mariana') |>
  group_by(species_name, sim, noise_size) |>
  reframe(A = area_thr(lambda)) |>
  ggplot() +
  aes(noise_size, A) +
  aes(group = sim, color = sim) +
  # geom_path(linewidth = .8, alpha = 0.8) +
  geom_smooth(level = NA) +
  facet_wrap(~species_name, scales = 'free') +
  scale_color_manual(
    values = c('#7fc97f', '#fb8072', '#fdc086')
  ) +
  theme_classic() +
  labs(
    y = 'Extinction risk',
    x = expression('Size of noise ('~sigma~')'),
    color = NULL
  ) +
  theme(
    strip.text= element_text(face = 'italic'),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'top'
  )

  

# save processed data
dir.create(file.path(sim_path, 'output_processed'))
saveRDS(lambdas, file.path(sim_path, 'output_processed', 'lambdas.RDS'))
saveRDS(outRF_pars, file.path(sim_path, 'output_processed', 'outRF_pars.RDS'))
saveRDS(outRF_cov, file.path(sim_path, 'output_processed', 'outRF_cov.RDS'))
