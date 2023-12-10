# script to set climate, competition, and initial conditions for sensitivity analysis
# Will Vieira
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# /!\ This should be ran locally (pre server run) /!\

library(tidyverse)
library(data.table)
set.seed(0.0)


# Prepare data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sim_path <- file.path('simulations', 'lambda_plot')

spIds <- read_csv(
  paste0(readLines('_data.path'), 'species_id.csv')
) |>
filter(sp_to_analyze)

treeData <- readRDS(
    paste0(readLines('_data.path'), 'treeData.RDS')
  )

treeData |>
  filter(status == 1) |>
  filter(dbh >= 127) |>
  select(!contains('_sc')) |>
  select(species_id, plot_id, year_measured, plot_size, dbh, contains('bio_')) |>
  drop_na() ->
treeData


# for the competition metric, I will use the variation in the size of each ind
# Then, to avoid loading the dataset for each simulation in the server,
# I will load the vector of conspecific and heterospecific individual sizes
load_sizeVec <- function(species_id, dbh)
{
  # Split dbh into conspecific and heterospecific groups
  con_dbh <- split(dbh, species_id)
  
  # get size of heterospecific individuals for each focal species
  het_dbh <- con_dbh |>
    imap(
      ~con_dbh[setdiff(names(con_dbh), .y)] |>
      unlist() |>
      unname()
    )
  
  tibble(
    con_dbh = con_dbh[species_id] |> setNames(nm = NULL),
    het_dbh = het_dbh[species_id] |> setNames(nm = NULL)
  )
}


treeData |>
  group_by(plot_id, year_measured) |>
  # filter plot that have at least one of the 31 modeled species
  mutate(has_modeled_species = any(species_id %in% spIds$species_id_old)) |>
  filter(has_modeled_species) |>
  select(!has_modeled_species) |>
  # get conspecific and heterospecific size per individual
  mutate(load_sizeVec(species_id, dbh)) |>
  group_by(species_id, plot_id, year_measured) |>
  # keep only one observation per species-plot-year
  slice_head(n = 1) |>
  filter(species_id %in% spIds$species_id_old) ->
input_data



# load parameters (population level)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pars_dir <- paste0(readLines('_data.path'), 'output_sim_processed')

models <- setNames(
  c(rep('intcpt_plot_comp_clim', 3), 'time_truc'),
  dir(pars_dir)
)

# function to load parameters from all species
load_all_pars <- function(vitalRate)
  spIds$species_id_old |>
    map(
      ~ readRDS(
          paste0(
            pars_dir, '/', vitalRate, '/', models[vitalRate], '/posterior_pop_', .x, '.RDS'
          )
        ) |>
        bind_cols(species_id = .x) |>
        filter(par != 'lp__')
    ) |>
    bind_rows() |>
    # add vital rate to the parameter name
    mutate(par = paste0(vitalRate, '.', par)) |>
    pivot_wider(names_from = par, values_from = value)


dir(pars_dir) |>
  map(load_all_pars) |>
  reduce(left_join, by = c('iter', 'species_id')) |>
  select(!iter) ->
population_pars

# # function to sample parameters from a given species
sample_pars <- function(sp, N = 100, dat = population_pars)
  dat |>
    filter(species_id == sp) |>
    select(!species_id) |>
    slice_sample(n = N) |>
    nest()

# input_data |>
#   mutate(pop_pars = sample_pars(species_id)) ->
# input_pars

# Code above is too slow, just saving the parameters and loading at every iteration
saveRDS(population_pars, file.path(sim_path, 'pop_pars.RDS'))

# load parameters (plot level)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to load random effects parameters for a single species
load_plot_pars <- function(sp)
  c('growth', 'mort', 'recruit') |>
    map(
      ~readRDS(
        paste0(
          pars_dir, '/', .x, '/', models[.x], '/posterior_plot_', sp, '.RDS'
        )
      ) |>
      select(!plot_id_seq) |>
      rename({{.x}} := value)
    ) |>
    reduce(left_join, by = c('iter', 'plot_id')) |>
    select(!iter)

sample_pars <- function(plot, N = 100, dat)
  dat |>
    filter(plot_id == plot) |>
    slice_sample(n = N) |>
    select(!plot_id)


plotYear_ls <- list()
for(Sp in spIds$species_id_old)
{
  plot_pars <- load_plot_pars(Sp)

  input_data |>
    filter(species_id == Sp) |>
    ungroup() |>
    select(plot_id, year_measured) |>
    filter(plot_id %in% unique(plot_pars$plot_id)) ->
  plot_year

  plot_year |>
    group_by(plot_id, year_measured) |>
    mutate(plot_pars = list(sample_pars(plot_id, dat = plot_pars))) ->
  plot_year

  plotYear_ls[[Sp]] = plot_year

  cat('Species', which(Sp == spIds$species_id_old), 'of', length(spIds$species_id_old))
}



# Save outputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Adding all plot random effects parameters to the input data will increase
# the time of loading the parameters RDS file for each batch file
# So here I will split the species level plot random effects in diff files
dir.create(file.path(sim_path, 'plot_parameters'))

plotYear_ls |>
  iwalk(
    ~ .x |>
      saveRDS(paste0(sim_path, '/plot_parameters/', .y, '.RDS'))
  )


# get all combination of species-plot-year from `plotYear_ls`
# to filter input_data with only plots that have random effects
plotYear_ls |>
  imap(~ .x |> select(plot_id, year_measured) |> bind_cols(species_id = .y)) |>
  bind_rows() ->
sp_plot_year


# Input parameters
input_data |>
  left_join(
    sp_plot_year |>
      bind_cols(toKeep = TRUE)
  ) |>
  filter(toKeep) |>
  select(!toKeep) |>
  saveRDS(file.path(sim_path, 'simulation_pars.RDS'))
