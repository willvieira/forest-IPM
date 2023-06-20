
library(tidyverse)

pars_sim <- read_csv('simulations/species_pair/output/simulation_pars.csv')

out <- parse_number(dir('simulations/species_pair/output/', pattern = 'RDS')) |>
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

# evolution of pop N over time
pdf(
  file = 'simulations/species_pair/N_over_time.pdf',
  width = 20,
  height = 16
)
for(Sp in unique(out$sp1)) {
  print(
    out |>
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


out_Ndiff <- out |>
  group_by(sp_pair, replication) |>
  mutate(
    diffN = N_sp1[year == max(year)] - N_sp2[year == max(year)]
  ) |>
  slice_head(n = 1)

pdf(
  file = 'simulations/species_pair/N_diff.pdf',
  width = 6,
  height = 10
)
for(Sp in unique(out$sp1)) {
  print(
    out_Ndiff |>
      filter(sp1 == Sp) |>
      ggplot(aes(diffN, y = fct_reorder(sp2, diffN), fill = stat(x))) +
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
  file = 'simulations/species_pair/N_diff_v2.pdf', 
  width = 30, 
  height = 10 
)
out_Ndiff |>
  mutate(sp1 = gsub('[0-9]', '', sp1), sp2 = gsub('[0-9]', '', sp2)) |>
  ggplot(aes(diffN, y = sp2, fill = stat(x))) +
    ggridges::geom_density_ridges_gradient(color = NA) +  
    scale_fill_gradient2() + 
    theme_minimal() + 
    geom_vline(xintercept = 0, linetype = 2) + 
    facet_wrap(~sp1, nrow = 1) +
    theme(
      legend.position = 'none',
      panel.spacing.x = unit(0, "lines")
    )
dev.off()



parse_number(dir('simulations/species_pair/output/', pattern = 'RDS')) |>
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
  file = 'simulations/species_pair/N_diff_summary.pdf', 
  width = 7, 
  height = 9 
)
out_Ndiff |>
  select(!contains('BA_')) |>
  ggplot(aes(diffN, fct_reorder(sp1, diffN), fill = stat(x))) +
    ggridges::geom_density_ridges_gradient() +
    scale_fill_gradient2() + 
    theme_minimal() + 
    geom_vline(xintercept = 0, linetype = 2) +
    ylab('') +
    xlab('Difference in population size between the focal and all other species') + 
    theme(legend.position = 'none')
dev.off()
