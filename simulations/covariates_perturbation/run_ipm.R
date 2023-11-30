# Script to run simulation on cluster

set.seed(0.0)

library(tidyverse)

dir.create('output')
dir.create('output_processed')

# create and export csv sim parameters
sim_pars <- readRDS('simulation_pars.RDS')

# Slurm returns an error for array above 10k
# loop over simulations to each batch is limited to 10k

i_max <- floor(nrow(sim_pars)/10000)
batch_max <- rep(10000, i_max)
batch_max[i_max + 1] <- nrow(sim_pars) - sum(batch_max)

# if # of sim > 100k, I cannot send all batches at the same time
if(length(batch_max) > 10) {

  number_of_sets <- ceiling(i_max/10)
  
  set = readline(
    prompt = paste0('The total number of simulations is larger than 100k,so you need to run this script ', number_of_sets, ' times. Tell me which time is this one? Valide values are [', paste0(seq_len(number_of_sets), collapse = ', '), ']')
  )

  if(!as.numeric(set) %in% seq_len(number_of_sets))
    stop('Answer is different than one of [', paste0(seq_len(number_of_sets), collapse = ', '),']')

  # first deal with past set of batch (merge and remove temp files)
  if(set > 1)
    source('merge_sim.R')

  # define sequence to run
  batch_seq = (0:9) + (10 * (as.numeric(set) - 1))
  if(as.numeric(set) == number_of_sets)
    batch_seq = batch_seq[1:which(batch_seq == length(batch_max))]

  # write slurm bash file
  for(i in batch_seq)
  {
    bash_file <- paste0('#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 23:50:00
#SBATCH --mem-per-cpu=1500M
#SBATCH --ntasks=1
#SBATCH --job-name=ipm', i, '
#SBATCH --mail-user=willian.vieira@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-', batch_max[i + 1], '

module load StdEnv/2020 r/4.1.0

BATCH=', i, ' R -f ipm_i.R
')

    # save bash file
    writeLines(bash_file, paste0('sub_', i, '.sh'))

    # run slurm
    system(paste0('sbatch sub_', i, '.sh'))
  }
}else{
  # write slurm bash file
  for(i in 0:i_max)
  {
    bash_file <- paste0('#!/bin/bash
#SBATCH --account=def-dgravel
#SBATCH -t 23:50:00
#SBATCH --mem-per-cpu=1500M
#SBATCH --ntasks=1
#SBATCH --job-name=ipm', i, '
#SBATCH --mail-user=willian.vieira@usherbrooke.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-', batch_max[i + 1], '

module load StdEnv/2020 r/4.1.0

BATCH=', i, ' R -f ipm_i.R
')

    # save bash file
    writeLines(bash_file, paste0('sub_', i, '.sh'))

    # run slurm
    system(paste0('sbatch sub_', i, '.sh'))
  }
}