##############################################################
# Get parameters from rstan output
# Will Vieira
# November 21, 2019
# last edited October 10, 2021
##############################################################



######################
 # For each vital rate:
  # - get spIds
  # - load rstan output and save summary matrix of all chains
  # - create data frame with parameters mean, 2.5% and 97.5%
  # - save data frame as a text file
 ######################


library(tidyverse)


spIds <- c(
  '28728ACERUB', '28731ACESAC', '183302PICMAR', '18032ABIBAL', '505490THUOCC',
  '19489BETPAP', '195773POPTRE', '19290QUEALB', '19481BETALL', '19408QUERUB',
  '183397TSUCAN', '19462FAGGRA', '18037PINTAE', '183385PINSTR', '18034PICRUB',
  '183295PICGLA', '183319PINBAN'
)

simNames <- setNames(
  c('bertalanffy_plotInd_BAspcTempPrec', 'mort_plotYear_sizeBAcompTempPrec', 'rec_m_p_BA'),
  c('growth', 'mort', 'rec')
)


dir.create('data')

post_pars <- map2_dfr(
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
          bind_cols(species_id = .x)
      ) |>
      group_by(species_id, par) |>
      summarise(
        qLower = quantile(value, probs = 0.025),
        qMedian = quantile(value, probs = 0.5),
        qUpper = quantile(value, probs = 0.975)
      ) |>
      ungroup() |>
      bind_cols(vitalRate = simName)
) |>
saveRDS('data/pars_summary.RDS', )
