# Script to check if chain finished `unexpectedly!`

library(tidyverse)

# load simulation pars
sim_pars = readRDS('simulation_pars.RDS')

# get file names
slurm_files = dir('.', pattern = 'slurm')

# test if the word `unexpectedly!` is present in the file 
slurm_files |>
  map(
    ~readLines(.x) |>
      grep(pattern = 'unexpectedly!')
  ) |>
  map(~ length(.x) > 0) |>
  unlist() ->
has_char

# sims to rerun
toRun = slurm_files[has_char]

# stop any running sim
for(i in 1:length(toRun)) {
  system(
    paste0(
      'scancel ',
      gsub('slurm-|.out', '', toRun[i])
    )
  )
  Sys.sleep(1)
}


# delete slurm files
file.remove(toRun)

# extract array_id
arrays_toRun = as.numeric(gsub('.out', '', gsub('.*_', '', toRun)))

# arrays_toRun = readRDS('toRun.RDS')[arrays_toRun]

# delete folders
for(i in arrays_toRun) {
  pars = sim_pars[i, ]
  sim_name_i = paste(pars$species_id, pars$border, pars$sim, pars$cond, sep = '_')
  unlink(
    file.path('output', sim_name_i),
    recursive = TRUE
  )
}

# save missing arrays
saveRDS(arrays_toRun, 'toRun.RDS')
