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
library(data.table)

spIds <- c(
  '28728ACERUB', '28731ACESAC', '183302PICMAR', '18032ABIBAL', '505490THUOCC',
  '19489BETPAP', '195773POPTRE', '19290QUEALB', '19481BETALL', '19408QUERUB',
  '183397TSUCAN', '19462FAGGRA', '18037PINTAE', '183385PINSTR', '18034PICRUB',
  '183295PICGLA', '183319PINBAN'
)

simNames <- setNames(
  c('bertalanffy_plotInd_BAspcTempPrecscl', 'mort_plot_sizeBAcompTempPrec', 'rec_m_p_BA_v2', 'sizeIng_time_truc'),
  c('growth', 'mort', 'rec', 'sizeIngrowth')
)


dir.create(file.path('data', 'parameters'))

# save all vital rate x species parameters as csv (population level parameters)
################################################################################
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



# save all individual-plot-year level parameters (slow code)
################################################################################

base_url = '../TreesDemography/output/'
treeData = readRDS('../TreesDemography/data/treeData.RDS')
paste0('data/parameters/randomEffects/', c('growth', 'mort', 'rec')) |>
  map(~ dir.create(.x, recursive = TRUE))

# Growth random effects (plot and individuals)
for(sp in spIds)
{
  trainData = readRDS(
    paste0(base_url, simNames['growth'], '/trainData_', sp, '.RDS')
  )
  
  # append plot_id info
  trainData[
    treeData,
    plot_id := i.plot_id,
    on = c('tree_id')
  ]

  # load plot parameters
  plot_pars <- readRDS(
      paste0(base_url, simNames['growth'], '/posteriorrPlot_', sp, '.RDS')
    ) |>
    mutate(plot_id_seq = parse_number(par)) |>
    select(!par) |>
    left_join(
      trainData[, .(plot_id_seq = unique(plot_id_seq)), by = plot_id]
    ) |>
    select(!plot_id_seq)
    
  # save raw posterior
  plot_pars |>
    saveRDS(
      paste0('data/parameters/randomEffects/growth/plot_', sp, '.RDS')
    )
  
  # save posterior mean
  plot_pars |>
    group_by(plot_id) |>
    reframe(value = mean(value)) |>
    saveRDS(
      paste0(
        'data/parameters/growth_plot_meanRandomEffect_',
        sp,
        '.RDS'
      )
    )

  # load tree parameters
  tree_pars <- readRDS(
      paste0(base_url, simNames['growth'], '/posteriorrTree_', sp, '.RDS')
    ) |>
    mutate(tree_id_seq = parse_number(par)) |>
    select(!par) |>
    left_join(
      trainData[, .(tree_id_seq = unique(tree_id_seq)), by = tree_id]
    ) |>
    select(!tree_id_seq)
    
  # save raw posterior
  tree_pars |>
    saveRDS(
      paste0('data/parameters/randomEffects/growth/tree_', sp, '.RDS')
    )

  # save posterior mean
  tree_pars |>
    group_by(tree_id) |>
    reframe(value = mean(value)) |>
    saveRDS(
      paste0(
        'data/parameters/growth_tree_meanRandomEffect_',
        sp,
        '.RDS'
      )
    )
}



# mort and recruitment plot random effects
for(vitalRate in c('mort', 'rec'))
{
  for(sp in spIds)
  {
    trainData = readRDS(
      paste0(
        base_url,
        simNames[vitalRate],
        ifelse(vitalRate == 'mort', '/trainData_', '/toSub_'),
        sp, '.RDS'
      )
    )
    
    # append plot_id info
    if(vitalRate == 'mort')
      trainData[
        treeData,
        plot_id := i.plot_id,
        on = c('tree_id')
      ]

    # load plot parameters
    plot_pars <- readRDS(
        paste0(
          base_url,
          simNames[vitalRate],
          '/posterior',
          ifelse(vitalRate == 'mort', 'psi', 'm'),
          'Plot_', sp, '.RDS')
      ) |>
      mutate(plot_id_seq = parse_number(par)) |>
      select(!par) |>
      left_join(
        trainData[, .(plot_id_seq = unique(plot_id_seq)), by = plot_id]
      ) |>
      select(!plot_id_seq)
      
    # save row samples
    plot_pars |>
      saveRDS(
        paste0(
          'data/parameters/randomEffects/',
          vitalRate,
          '/plot_', sp, '.RDS'
        )
      )

    # save posterior mean
    plot_pars |>
      group_by(plot_id) |>
      reframe(value = mean(value)) |>
      saveRDS(
        paste0(
          'data/parameters/',
          vitalRate,
          '_plot_meanRandomEffect_',
          sp,
          '.RDS'
        )
      )
  }
}
