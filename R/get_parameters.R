##############################################################
# Get parameters from cmdstanr output
# Will Vieira
# November 21, 2019
# last edited March 29, 2023
##############################################################



######################
  # /!!\ Parameters are available in local machine only /!!\
 # For each vital rate:
  # - get spIds
  # - load cmdstanr output
  # - save as iter x parameter matrix for each species x vital rate
 ######################


library(tidyverse)
library(data.table)
library(cmdstanr)
library(posterior)

simNames <- c(
  'growth' = 'bertalanffy_plotYear_BAspcTempPrecscl_fulldb',
  'mort' = 'mort_plotYear_sizeBAcompTempPrec_alldb',
  'rec' = 'rec_m_p_BA_plotYear',
  'sizeIngrowth' = 'sizeIng_time_truc_fulldb'
)

spIds <- list.dirs(
  file.path('..', 'TreesDemography', 'output', simNames[1]),
  full.names = FALSE,
  recursive = FALSE
)


dir.create(file.path('data', 'parameters'), recursive = TRUE)

# save all vital rate x species parameters as csv (population level parameters)
################################################################################

# list with parameter names for each vital rate
parNames <- list(
  'growth' = c(
    'r', 'sigma_plot', 'sigma_year', 'sigma_obs', 'Lmax', 'Beta',
    'theta', 'optimal_temp', 'tau_temp', 'optimal_prec', 'tau_prec'
  ),
  'mort' = c(
    'psi', 'sigma_plot', 'sigma_year', 'size_opt', 'size_var', 'Beta',
    'theta', 'optimal_temp', 'tau_temp', 'optimal_prec', 'tau_prec'
  ),
  'rec' = c(
    'mPop_log', 'sigma_plot', 'sigma_year', 'p_log', 'beta_p',
    'optimal_BA', 'sigma_BA'
  ),
  'sizeIngrowth' = c(
    'size_int', 'phi_time', 'sigma_size'
  )
)

# growth, mortality, and recruitment models
map(
  names(simNames),
  ~map2_dfr(
    spIds,
    .x,
    ~read_cmdstan_csv(
      dir(
        file.path(
          '..', 'TreesDemography', 'output', simNames[.y], .x
        ),
        full.names = TRUE,
        pattern = .y
      ),
      variables = parNames[[.y]]
    )$post_warmup_draws |>
    as_draws_df() |>
    mutate(
      species_id = .x,
      iter = row_number()
    ) |>
    select(!c(.chain, .iteration, .draw)) |>
    write_csv(
      file.path(
        'data', 'parameters',
        paste0(.y, '_', .x, '.csv'))
    )
  )
)
