# Parameters vs output (sensitivity analysis)

library(tidyverse)
library(ggdist)

# load parameter simulations
sim_pars <- read_csv('simulations/sensAnalysis/output/simulation_pars.csv') |>
  mutate(
    clim = fct_recode(clim, border = 'bad', optimal = 'good')
  )

out <- parse_number(
    dir(
      'simulations/sensAnalysis/output/', pattern = '.RDS'
    )
  ) |>
  sort() |>
  map_dfr(
    ~readRDS(
      paste0('simulations/sensAnalysis/output/sim_', .x, '.RDS')
    )[['year_summ']] |>
    bind_cols(sim = .x)
  ) |>
  # add parameter simulations
  left_join(
    sim_pars |>
      select(!c(param_method, seed, contains('random'))) |>
      mutate(sim = row_number())
  ) |>
  as_tibble()

out_pars <- parse_number(
    dir(
      'simulations/sensAnalysis/output/', pattern = '.RDS'
    )
  ) |>
  sort() |>
  map_dfr(
    ~readRDS(
      paste0('simulations/sensAnalysis/output/sim_', .x, '.RDS')
    )[[1]] |>
    unlist() |>
    enframe() |>
    bind_cols(sim = .x)
  )

lambdas <- parse_number(
    dir(
      'simulations/sensAnalysis/output/', pattern = '.RDS'
    )
  ) |>
  sort() |>
  map_dfr(
    ~readRDS(
      paste0('simulations/sensAnalysis/output/sim_', .x, '.RDS')
    )[['lambda']] |>
    enframe(name = 'lambda') |>
    bind_cols(sim = .x)
  ) |>  
  left_join( 
    sim_pars |> 
      select(!c(param_method, seed, contains('random'))) |> 
      mutate(sim = row_number()) 
  )





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some exploraratory figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(
  file = 'simulations/sensAnalysis/dist_lambda_N_BA.pdf',
  width = 16,
  height = 12
)
# Initial lambda in function of climate and competition
lambdas |> 
  filter(lambda == 1) |> 
  mutate(
    species_id = gsub('[0-9]', '', species_id)
  ) |>
  ggplot(aes(value, color = comp, linetype = clim)) +
    geom_density() +
    facet_wrap(~species_id, scales = 'free') +
    geom_vline(xintercept = 1) +
    theme_minimal() +
    xlab('Lambda')

# Final Total BA in func of clim and competition
out |> 
  group_by(sim) |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(
    species_id = gsub('[0-9]', '', species_id)
  ) |>
  ggplot(aes(BA_sp1, color = comp, linetype = clim)) +
    geom_density() +
    facet_wrap(~species_id, scales = 'free') +
    theme_minimal() +
    xlab('Basal area at equilibrium')

# Final Total N in func of clim and competition
out |> 
  group_by(sim) |>
  filter(year == max(year)) |>
  ungroup() |>
  mutate(
    species_id = gsub('[0-9]', '', species_id)
  ) |>
  ggplot(aes(N_sp1, color = comp, linetype = clim)) +
    geom_density() +
    facet_wrap(~species_id, scales = 'free') +
    theme_minimal() +
    xlab('Population size at equilibrium')

# Total BA in func of clim and competition
out |> 
  group_by(sim) |>
  filter(year == max(year)) |>
  ungroup() |>
  left_join(
    lambdas |>
      filter(lambda == 1)
  ) |>
  filter(BA_sp1 > 1) |>
  ggplot(aes(BA_sp1, value, color = comp, shape = clim)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~species_id, scales = 'free') +
    theme_minimal() +
    xlab('Total basal area') +
    ylab('Lambda')

# Summary of lambda mean and variance across species, competition, and climate
lambdas |> 
  filter(lambda == 1) |> 
  group_by(species_id, comp, clim) |>
  reframe(
    `mean lambda` = mean(value),
    `sd lambda` = sd(value)
  ) |>
  pivot_longer(
    cols = c(`mean lambda`, `sd lambda`)
  ) |>
  ggplot(aes(comp, value, fill = clim)) +
    geom_boxplot() +
    facet_wrap(~name, scales = 'free') +
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
  filter(lambda == 1) |> # initial lambda == 1, final lambda == 2
  select(!c(lambda, n_time, deltaTime, plotSize)) |>
  rename(lambda = value) |>
  left_join(
    out_pars |>
      pivot_wider()
  )

out_dt <- tibble()
count = 1;totalCount <- length(unique(rfdat$species_id)) * 20 * 4
for(Sp in unique(rfdat$species_id))
{
  rfdat_sp <- rfdat |>
      filter(species_id == Sp) |>
      select(!c(sim, species_id, replication))

  for(Rep in 1:20)
  {
    for(Comp in c('low', 'high'))
    {
      for(Clim in c('optimal', 'border'))
      {
        RF = ranger(
          lambda ~ .,
          data = rfdat_sp |>
            filter(comp == Comp & clim == Clim) |>
            select(!c(comp, clim)),
          importance = 'permutation'
        )
        import <- importance(RF)/sum(importance(RF))
        
        out_dt <- rbind(
          out_dt,
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
out_dt <- out_dt |>
  mutate(parvr = gsub('\\..*', '', par)) |>
  group_by(species_id, comp, clim, rep, parvr) |>
  mutate(
    impvr = sum(imp)
  ) |>
  ungroup()


# Figures
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Importance partition between vital rares for the RandomForest 20 replications
pdf(
  file = 'simulations/sensAnalysis/importance.pdf',
  width = 12,
  height = 10
)
for(Sp in unique(out_dt$species_id))
{
  dt_sp <- out_dt |>
      filter(species_id == Sp) |>
      filter(par != 'Rep.Loc') |>
      #filter(clim == 'good' & comp == 'low') |>
      group_by(comp, clim, rep, parvr) |>
      slice_head(n = 1)

  R2 <- dt_sp |>
    group_by(comp, clim) |>
    reframe(
      R2 = mean(R2)
    ) |>
    mutate(
      name = paste0(comp, ', ', clim, ': ', round(R2, 2))
    ) |>
    pull(name) |>
    paste0(collapse = '; ')
   
   print(
    dt_sp |>
      ggplot(aes(y = as.factor(rep), x = impvr, fill = parvr)) +
        geom_bar(stat = 'identity') +
        facet_wrap(
          ~comp + clim,
          labeller = labeller(
            comp = setNames(
              paste0('Competition: ', c('low', 'high')),
              c('low', 'high')
            ),
            clim = setNames(
              paste0('Climate: ', c('border', 'optimal')),
              c('border', 'optimal')
            )
          )
        ) +
        xlab('Relative importance') +
        ylab('Replication') +
        labs(
          subtitle = paste0(Sp, '\nAverage R2: ', R2)
        )
  )
}
dev.off()




# Simplex plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# I need is to get the proportion of Growth, Mort, and Recruit
# First issue is that I have a forth variable (sizeIngrowth) which
# is insignifcant close to the others. When removing it, should I normalize
# the remaining 3 vars to 1 or keep their orinial values?
#(I start by the former)
out_dt |>
  filter(par != 'Rep.Loc') |>
  group_by(comp, clim, rep, parvr) |>
  slice_head(n = 1) |>
  ungroup() |>
  filter(parvr == 'sizeIngrowth') |>
  pull(impvr) |>
  quantile()
# Importance of size ingrowth is very low even in the best case scenario
# 100% upper quantile = 3.1%, so we can ignore it and just use the
# growth, mortality, and recruitment vital rates

library(ggtern)

pdf(
  file = 'simulations/sensAnalysis/importance_summ.pdf',
  width = 14,
  height = 10
)
out_dt |>
  filter(par != 'Rep.Loc') |>
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
  ggtern(
    aes(x = growth, y = mort, z = rec, color = species_id)
  ) +
  geom_point() +
  facet_wrap(
    ~comp + clim,
    labeller = labeller(
      comp = setNames(
        paste0('Competition: ', c('low', 'high')),
        c('low', 'high')
      ),
      clim = setNames(
        paste0('Climate: ', c('border', 'optimal')),
        c('border', 'optimal')
      )
    )
  )


# Mortality vs recruitment
out_dt |>
  filter(par != 'Rep.Loc') |>
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
    read_csv('data//species_id.csv') |> 
      select(species_id_old, group, shade) |> 
      rename(species_id = species_id_old)
  ) |>
  ggplot(aes(mort, rec, color = shade)) +
    geom_point() +
    facet_wrap(~clim+comp) +
    theme_minimal()
  
dev.off()


