
library(tidyverse)

pars_sim <- read_csv('simulations/species_pair/output/simulation_pars.csv')

out_pair <- parse_number(dir('simulations/species_pair/output/', pattern = 'RDS')) |>
  sort() |>
  map_dfr(
    ~ readRDS(
      paste0('simulations/species_pair/output/sim_', .x, '.RDS')
    )[['year_summ']] |>
    bind_cols(sim = .x)
  ) |>
  left_join(
    pars_sim |>
      mutate(sim = row_number()) |>
      select(
        !c(contains('Ninit_'), n_time, deltaTime, param_method, plot_random, seed)
      )
  ) |>
  mutate(
    sp_pair = paste0(sp1, '_', sp2)
  ) |>
  as_tibble()

pars_sim_re <- read_csv('simulations/species_pair_randomEffect/output/simulation_pars.csv')

out_pair_re <- parse_number(dir('simulations/species_pair_randomEffect/output/', pattern = 'RDS')) |>
  sort() |>
  map_dfr(
    ~ readRDS(
      paste0('simulations/species_pair_randomEffect/output/sim_', .x, '.RDS')
    )[['year_summ']] |>
    bind_cols(sim = .x)
  ) |>
  left_join(
    pars_sim_re |>
      mutate(sim = row_number()) |>
      select(
        !c(contains('Ninit_'), n_time, deltaTime, param_method, seed)
      )
  ) |>
  mutate(
    sp_pair = paste0(sp1, '_', sp2)
  ) |>
  as_tibble()


# evolution of pop N over time
pdf(
  file = 'simulations/species_pair_randomEffect/N_over_time.pdf',
  width = 20,
  height = 16
)
for(Sp in unique(out_pair_re$sp1)) {
  print(
    out_pair_re |>
      filter(sp1 == Sp) |>
      select(!contains('BA_')) |>
      pivot_longer(
        cols = contains('N_')
      ) |>
      ggplot(
        aes(
          year,
          value,
          color = name,
          group = interaction(replication, name)
        )
      ) +
      geom_line(alpha = 0.2) +
      facet_wrap(~sp2, scales = 'free') +
      theme_minimal() +
      labs(title = paste0('Focal species (sp1): ', Sp)) +
      theme(legend.position = 'top')
  )
}
dev.off()


out_Ndiff_pair <- out_pair |>
  group_by(sp_pair, replication) |>
  mutate(
    diffN = N_sp1[year == max(year)] - N_sp2[year == max(year)]
  ) |>
  slice_head(n = 1)

out_Ndiff_pair_re <- out_pair_re |>
  group_by(sp_pair, replication) |>
  mutate(
    diffN = N_sp1[year == max(year)] - N_sp2[year == max(year)]
  ) |>
  slice_head(n = 1)


pdf(
  file = 'simulations/species_pair_randomEffect/N_diff.pdf',
  width = 6,
  height = 10
)
for(Sp in unique(out$sp1)) {
  print(
    out_Ndiff_pair_re |>
      filter(sp1 == Sp) |>
      ggplot(aes(diffN, y = fct_reorder(sp2, diffN), fill = after_stat(x))) +
        ggridges::geom_density_ridges_gradient() + 
        scale_fill_gradient2() +
        theme_minimal() +
        geom_vline(xintercept = 0, linetype = 2) +
        xlab('Difference in population size at equilibrium') +
        ylab('') +
        theme(legend.position = 'none') +
        labs(
          title = paste0(
            'Competition between focal species: ', Sp
          )
        )
  )
}
dev.off()


pdf(
  file = 'simulations/species_pair_randomEffect/N_diff_re_vs_nore.pdf',
  width = 7,
  height = 10
)
for(Sp in unique(out$sp1)) {
  print(
    out_Ndiff_pair_re |>
      filter(sp1 == Sp) |>
      bind_cols(re = 'Yes') |>
      bind_rows(
        out_Ndiff_pair |>
          filter(sp1 == Sp) |>
          bind_cols(re = 'No')
      ) |>
      ggplot(aes(diffN, y = fct_reorder(sp2, diffN), fill = re)) +
        ggridges::geom_density_ridges2(alpha = 0.7, color = NA) + 
        #scale_fill_gradient2() +
        theme_minimal() +
        geom_vline(xintercept = 0, linetype = 2) +
        xlab('Difference in population size at equilibrium') +
        ylab('') +
        labs(
          title = paste0(
            'Competition between focal species: ', Sp
          ),
          fill = 'Rand effect'
        )
  )
}
dev.off()




parse_number(dir('simulations/species_pair_randomEffect/output/', pattern = 'RDS')) |>
  sort() |>
  map_dbl(
    ~ readRDS(
      paste0('simulations/species_pair/output/sim_', .x, '.RDS')
    )[['sim_time']] |>
    as.numeric(units = 'hours')
  ) |>
  hist()



# distribution across species
pdf( 
  file = 'simulations/species_pair_randomEffect/N_diff_summary.pdf', 
  width = 7, 
  height = 9 
)
out_Ndiff_pair_re |>
  bind_cols(re = 'Yes') |>
  bind_rows(
    out_Ndiff_pair |>
    bind_cols(re = 'No')
  ) |>
  select(!contains('BA_')) |>
  ggplot(aes(diffN, fct_reorder(sp1, diffN), fill = re)) +
    ggridges::geom_density_ridges2(color = NA, alpha = 0.7) +
    theme_minimal() + 
    geom_vline(xintercept = 0, linetype = 2) +
    ylab('') +
    xlab('Difference in population size between the focal and all other species') +
    labs(fill = 'Rand effect')
dev.off()
