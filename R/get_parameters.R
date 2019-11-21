##############################################################
# Get parameters from rstan output
# Will Vieira
# November 21, 2019
##############################################################



######################
 # For each vital rate:
  # - get spIds
  # - load rstan output and save summary matrix of all chains
  # - create data frame with parameters mean, 2.5% and 97.5%
  # - save data frame as a text file
 ######################



# Setup
  library(rstan)
  # species
  spIds <- '18032-ABI-BAL'
#



# Read simulation
  for(sp in spIds)
  {
    for(vRate in c('growth', 'mort', 'fec'))
    {
      # load simulation
      sim <- readRDS('')

      # Output
      # Get summary output for all chains
      simSummary <- summary(sim)[[1]]
      # drop `lp__` value
      simSummary <- simSummary[-nrow(simSummary), ]


      # Create parameters data frame
      params <- data.frame(
        simSummary[, 'mean'],
        simSummary[, '2.5%'],
        simSummary[, '97.5%']
      )

      # rename cols
      colnames(params) <- c('mean', 'lo2.5', 'up97.5')

      # save text  file
      if(!dir.exists('params')) dir.create('params')
      write.table(params, file = paste0('params/pars_', vRate, '_', sp, '.txt'))
    }
  }
#
