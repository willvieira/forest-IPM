# Parameters vs output (sensitivity analysis)

library(tidyverse)
library(ggdist)
library(ggtext)
library(ggrepel)
library(ggtern)


# load parameters and output
sim_path <- file.path('simulations', 'sensAnalysis_V3')

readRDS(file.path(sim_path, 'output_complete.RDS')) |>
  mutate(sim = 1:n()) |>
  left_join(
    read_csv(file.path(sim_path, 'simulation_pars.csv')) |>
      mutate(sim = 1:n()) |>
      select(species_id, comp, clim, replication, sim)
  ) |>
  select(!sim) ->
lambdas


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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some exploraratory figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(
  file = file.path(sim_path, 'dist_lambda.pdf'),
  width = 16,
  height = 12
)
# Initial lambda in function of climate and competition
lambdas |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  mutate(
    `Temperature\ncondition` = factor(clim, levels = c('cold', 'center', 'hot')),
    `Competition` = comp
  ) |>
  ggplot(aes(log(lambda), color = `Temperature\ncondition`, linetype = Competition)) +
    geom_density(linewidth = .8) +
    facet_wrap(~species_name, scales = 'free') +
    geom_vline(xintercept = 0) +
    theme_classic() +
    scale_color_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
    labs(
      x = expression('ln('~lambda~')'),
      y = ''
    ) +
    theme(
      strip.text = element_text(face = 'italic'),
      strip.background = element_blank()
    )

# Summary of lambda mean and variance across species, competition, and climate
lambdas |>
  mutate(
    `Temperature\ncondition` = factor(clim, levels = c('cold', 'center', 'hot')),
    Competition = comp
  ) |>
  group_by(species_id, Competition, `Temperature\ncondition`) |>
  reframe(
    `mean lambda` = mean(lambda),
    `sd lambda` = sd(lambda)
  ) |>
  pivot_longer(
    cols = c(`mean lambda`, `sd lambda`)
  ) |>
  ggplot(aes(Competition, value, fill = `Temperature\ncondition`)) +
    geom_boxplot() +
    facet_wrap(~name, scales = 'free') +
    scale_fill_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
    theme_minimal() +
    xlab('Competition intensity') +
    ylab('') +
    labs(fill = 'Climate')

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random forest
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ranger)

rfdat <- lambdas |>
  select(!replication)

outRF_pars <- tibble()
count = 1;totalCount <- length(unique(rfdat$species_id)) * 20 * 6
for(Sp in unique(rfdat$species_id))
{
  rfdat_sp <- rfdat |>
      filter(species_id == Sp) |>
      select(!species_id)

  for(Rep in 1:20)
  {
    for(Comp in c('low', 'high'))
    {
      for(Clim in c('cold', 'center', 'hot'))
      {
        RF = ranger(
          lambda ~ .,
          data = rfdat_sp |>
            filter(comp == Comp & clim == Clim) |>
            select(!c(comp, clim)),
          importance = 'permutation'
        )
        import <- importance(RF)/sum(importance(RF))
        
        outRF_pars <- rbind(
          outRF_pars,
          tibble(
            species_id = Sp,
            rep = Rep,
            comp = Comp,
            clim = Clim,
            R2 = RF$r.squared,
            par = names(import),
            imp = import
          )
        )

        cat('Progress', round(count/totalCount * 100, 1), '%\r')
        count = count + 1
      }
    }
  }
}

# calculate vital rate importance
outRF_pars <- outRF_pars |>
  mutate(parvr = gsub('\\..*', '', par)) |>
  group_by(species_id, comp, clim, rep, parvr) |>
  mutate(
    impvr = sum(imp)
  ) |>
  ungroup()


# Random forest for predictive variables
outRF_cov <- tibble()
count = 1;totalCount <- length(unique(rfdat$species_id)) * 20
for(Sp in unique(rfdat$species_id))
{
  rfdat_sp <- rfdat |>
      filter(species_id == Sp) |>
      select(lambda, comp, clim)

  for(Rep in 1:20)
  {
    RF = ranger(
      lambda ~ .,
      data = rfdat_sp,
      importance = 'permutation'
    )
    import <- importance(RF)/sum(importance(RF))
    
    outRF_cov <- rbind(
      outRF_cov,
      tibble(
        species_id = Sp,
        rep = Rep,
        R2 = RF$r.squared,
        par = names(import),
        imp = import
      )
    )

    cat('Progress', round(count/totalCount * 100, 1), '%\r')
    count = count + 1
  }
}



# Simplex plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# I need is to get the proportion of Growth, Mort, and Recruit
# First issue is that I have a forth variable (sizeIngrowth) which
# is insignifcant close to the others. When removing it, should I normalize
# the remaining 3 vars to 1 or keep their orinial values?
#(I start by the former)
outRF_pars |>
  group_by(comp, clim, rep, parvr) |>
  slice_head(n = 1) |>
  ungroup() |>
  filter(parvr == 'sizeIngrowth') |>
  pull(impvr) |>
  quantile()
# Importance of size ingrowth is very low even in the best case scenario
# 100% upper quantile = 1.19%, so we can ignore it and just use the
# growth, mortality, and recruitment vital rates

pdf(
  file = file.path(sim_path, 'importance_summ.pdf'),
  width = 14,
  height = 10
)

outRF_pars |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  group_by(species_name, rep, comp, clim) |>
  slice_head(n = 1) |>
  ggplot() +
  aes(R2, fct_reorder(species_name, R2)) +
  stat_pointinterval() +
  theme_classic() +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = 'none'
  ) +
    theme(
    axis.text.y = element_text(face = "italic"),
    panel.grid.major.y = element_line(colour = rgb(0,0,0,.1))
  ) +
  labs(
    x = expression(R^2),
    y = ''
  )


outRF_pars |>
  group_by(species_id, comp, clim, rep, parvr) |>
  slice_head(n = 1) |>
  select(!c(par, imp, rep)) |>
  group_by(species_id, comp, clim, parvr) |>
  reframe(
    R2 = mean(R2),
    impvr = mean(impvr)
  ) |>
  pivot_wider(
    names_from = parvr,
    values_from = impvr
  ) |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  mutate(
    clim = factor(clim, levels = c('cold', 'center', 'hot')),
    comp = paste0(comp, ' competition')
  ) |>
  ggtern(
    aes(x = growth, y = mort, z = rec, color = shade_sylvics)
  ) +
  geom_point() +
  facet_grid(
    comp ~ clim
  ) +
  scale_color_manual(
    values = c("#20bc45", "#87bc45", "#edbf33", "#ea5545", "#ba0000")
  ) +
  theme(legend.position = 'top') +
  labs(color = NULL)


# Mortality vs recruitment  
outRF_pars |>
  group_by(species_id, comp, clim, rep, parvr) |>
  slice_head(n = 1) |>
  select(!c(par, imp, rep)) |>
  group_by(species_id, comp, clim, parvr) |>
  reframe(
    R2 = mean(R2),
    impvr = mean(impvr)
  ) |>
  pivot_wider(
    names_from = parvr,
    values_from = impvr
  ) |>
  left_join(
    spIds,
    by = c('species_id' = 'species_id_old')
  ) |>
  mutate(
    clim = factor(clim, levels = c('cold', 'center', 'hot')),
    comp = paste0(comp, ' competition'),
    comp = factor(comp, levels = c('low competition', 'high competition')),
    
  ) |>
  ggplot() +
  aes(mort, rec, fill = clim, color = R2) +
  facet_wrap(~comp) +
  geom_hdr(
    aes(color = NULL),
    probs = .9, alpha = .4,
    method = 'mvnorm'
  ) +
  geom_point() +
  theme_classic() +
  scale_fill_manual(values = c('#91bfdb', '#99d594', '#fc8d59')) +
  scale_color_binned(type = "viridis") +
  labs(
    x = 'Mortality importance',
    y = 'Recruitment importance',
    color = 'Growth\nimportance',
    fill = 'Temperature\nrange'
  )

  # Predictive variables
  outRF_cov |>
    filter(par == 'clim') |>
    left_join(
      spIds,
      by = c('species_id' = 'species_id_old')
    ) |>
    group_by(species_name) |>
    mutate(
      R2_mean = mean(R2),
      imp_mean = mean(imp)
    ) |>
    slice_head(n = 1) |>
    ggplot() +
    aes(imp, R2) +
    geom_point() +
    geom_text_repel(
      aes(label = species_name),
      alpha = 0.8,
      size = 2.3,
      fontface = 'italic'
    ) +
    theme_classic() +
    geom_hline(yintercept = 0.6, alpha = 0.4, linetype = 2) +
    labs(
      y = expression(R^2),
      x = 'Climate vs competition relative importance'
    )

dev.off()


# save processed data
dir.create(file.path(sim_path, 'output_processed'))
saveRDS(lambdas, file.path(sim_path, 'output_processed', 'lambdas.RDS'))
saveRDS(outRF_pars, file.path(sim_path, 'output_processed', 'outRF_pars.RDS'))
saveRDS(outRF_cov, file.path(sim_path, 'output_processed', 'outRF_cov.RDS'))

  list(
    lambdas = lambdas,
    outRF_pars = out_dt
    outRF_cov = out_pred
  ),
  file.path(sim_path, 'output_processed', 'output.RDS')
)