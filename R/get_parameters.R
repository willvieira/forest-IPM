##############################################################
# Get parameters from cmdstanr output
# Will Vieira
# November 21, 2019
# last edited March 29, 2023
##############################################################



######################
 # For each vital rate:
  # - get spIds
  # - load cmdstanr output
  # - save as iter x parameter matrix for each species x vital rate
 ######################


library(tidyverse)


spIds <- c(
  '28728ACERUB', '28731ACESAC', '183302PICMAR', '18032ABIBAL', '505490THUOCC',
  '19489BETPAP', '195773POPTRE', '19290QUEALB', '19481BETALL', '19408QUERUB',
  '183397TSUCAN', '19462FAGGRA', '18037PINTAE', '183385PINSTR', '18034PICRUB',
  '183295PICGLA', '183319PINBAN'
)

simNames <- setNames(
  c('bertalanffy_plotInd_BAspcTempPrecscl', 'mort_plot_sizeBAcompTempPrec', 'rec_m_p_BA'),
  c('growth', 'mort', 'rec')
)


dir.create(file.path('data', 'parameters'))

map2(
  simNames,
  names(simNames),
  function(sim, simName)
    map_dfr(
        spIds,
        ~ readRDS(
          paste0(
            '../TreesDemography/output/',
            sim,
            '/posteriorPop_',
            .x,
            '.RDS'
          )
        ) |>
        pivot_wider(
          names_from = par,
          values_from = value
        ) |>
        write_csv(
          file.path(
            'data', 'parameters',
            paste0(simName, '_', .x, '.csv'))
        )
    )
)
