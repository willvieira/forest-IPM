# Prepare data for running Stan linear model of lambda against temperature
# Will Vieira
# Dec 2023


library(tidyverse)


# Load parameters list
readRDS('simulations/lambda_plot/simulation_pars.RDS') |>
  ungroup() |>
  mutate(array_id = row_number()) |>
  # determine if plot is upper or below mid temperature range
  group_by(species_id) |>  
  mutate( 
    mid_temp = (max(bio_01_mean) + min(bio_01_mean))/2, 
    border = if_else(bio_01_mean <= mid_temp, 'cold', 'hot') 
  ) |>
  ungroup() |>
  # determine if FIA or quebec plot
  # 1: FIA
  # 2: Quebec
  mutate(db_origin = if_else(plot_size < 200, 1, 2)) |>
  select(array_id, species_id, bio_01_mean, border, db_origin) ->
pars

# Load output of lambda simulation-condition (12 total per species-plot-year)
out <- readRDS('simulations/lambda_plot/final_outputv2.RDS')

# Compute number of plots per species-boder group
pars |>
  group_by(species_id, border, db_origin) |> 
  count() |>
  group_by(species_id, border) |>
  filter(db_origin == db_origin[which.max(n)]) |>
  arrange(desc(n)) |>
  as.data.frame()

# sample plots per species to avoid long Stan run
# First run using 5k plots
# (this only filters 7 of the 62 total species-border combinations)
nb_plots <- 5e3


# for each species-boder, sample `nb_plots` plots, and filter out
# object to keep only necessary info for the specific species-boder Stan run
out_dir <- 'simulations/model_lambdaPlot/out_species/'
dir.create(out_dir, recursive = TRUE)

count = 1
for(Sp in unique(pars$species_id))
{
  for(brd in c('cold', 'hot'))
  {    
    # get arrays to export
    pars |>
      filter(species_id == Sp) |>
      filter(border == brd) |>
      # When species is present at both FIA and QC database
      # Keep only the plots from the dataset that has most plots
      mutate(max_origin = as.numeric(names(which.max(table(db_origin))))) |>
      filter(db_origin == max_origin) |>
      slice_sample(n = nb_plots) ->
    pars_spBorder

    pars_spBorder |>
      pull(array_id) ->
    arrays_to_load

    # filter output object to specific species-boder
    out |>
      filter(array_id %in% arrays_to_load) ->
    out_spBorder
  
    # loop over simulations [0 to 3] and conditions [1 to 3]
    for(Sim in 1:2)
    {
      for(Cond in 1:3)
      {
        out_spBorder |>
          filter(sim == Sim) |>
          filter(cond == Cond) |>
          left_join(
            pars_spBorder |>
              select(array_id, bio_01_mean),
            by = 'array_id'
          ) |>
          saveRDS(
            paste0(out_dir, paste(Sp, brd, Sim, Cond, sep = '_'), '.RDS')
          )
        
        # progress
        cat(' Saving ', round(count/(31 * 2 * 2 * 3) * 100, 1), '%\r')
        count = count + 1
      }
    }
    rm(out_spBorder)
  }
}
