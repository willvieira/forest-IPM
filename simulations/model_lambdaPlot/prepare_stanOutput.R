# Script to format stan parameter output

library(tidyverse)
library(cmdstanr)
library(posterior)

sim_path = 'simulations/model_lambdaPlot'

sim_pars <- readRDS(file.path(sim_path, 'simulation_pars.RDS'))

bad_chains <- read_csv(file.path(sim_path, 'chains_toRm.csv'))

out_pars = list()
out_diag = list()

arrays = 1:nrow(sim_pars)

for(i in arrays) {

  try({
  pars_i = sim_pars[i, ]

  bad_chains_i <- bad_chains |>
    filter(species_id == pars_i$species_id & border == pars_i$border & cond == pars_i$cond) |>
    pull(chain_rm) |>
    unique()
  
  md_out <- as_cmdstan_fit(
    dir(
      file.path(
        sim_path,
        'output',
        pars_i |> select(!array_id) |> paste0(collapse = '_')
      ),
      full.names = TRUE
    )
  )

  diag_out <- list(
    diag_summary = md_out$diagnostic_summary(),
    rhat = md_out$summary() |> select(variable, rhat),
    time = md_out$time()
  )

  out_diag[[i]] = diag_out

  md_out <- read_cmdstan_csv(
      dir(
      file.path(
        sim_path,
        'output',
        pars_i |> select(!array_id) |> paste0(collapse = '_')
      ),
      full.names = TRUE
    )
  )

  md_out$post_warmup_draws |>
    as_draws_df() |>
    as_tibble() |>
    filter(!.chain %in% bad_chains_i) |>
    filter(.iteration %in% 500:1000) |>
    slice_sample(n = 500) |>
    select(!contains('.')) |>
    mutate(iter = row_number()) |>
    pivot_longer(
      cols = !iter,
      names_to = 'par'
    ) |>
    bind_cols(array_id = i) ->
  params

  out_pars[[i]] = params
  })

  cat(' Progress ', round(i/length(arrays) * 100, 1), '%\r')
}
  

out_pars |>
  bind_rows() |>
  saveRDS(file.path(sim_path, 'param_posterior.RDS'))

saveRDS(out_diag, file.path(sim_path, 'out_diag.RDS'))
