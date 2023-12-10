# script to merge simulations for a set of batch

# turn off scientific notation
options(scipen = 999)


set = as.numeric(set)

cat('Merging set ', set, '\n')

# First check that all simulations from last set were present
arrays_from_last_set <- 1:1e5 + ((set - 1) * 1e5)
if(max(arrays_from_last_set) > nrow(sim_pars))
  arrays_from_last_set = arrays_from_last_set[1:which(arrays_from_last_set == nrow(sim_pars))]

# get array from output folder
arrays_sim <- parse_number(dir('output'))
missing_arrays <- setdiff(arrays_from_last_set, arrays_sim)

stop if there are missing arrays
if(length(missing_arrays)) {
  saveRDS(missing_arrays, paste0('missingArrays_set_', set, '.RDS'))
  stop(paste0('There are ', length(missing_arrays), 'missing arrays. Their IDs are saved in the `missingArrays_set_', set, '.RDS` file.'))
}

# merge output arrays into single file
arrays_sim |>
  map(
    ~readRDS(paste0('output/sim_', .x, '.RDS')) |>
      bind_cols(array_id = .x)
  ) |>
  bind_rows() |>
  saveRDS(paste0('output_processed/set_', set, '.RDS'))


cat('Deleting set slurm files', set, '\n')
files_to_rm <- list.files('.', pattern = 'slurm-')
file.remove(files_to_rm)
files_to_rm <- list.files('.', pattern = 'sub_')
file.remove(files_to_rm)
