## code to prepare internal species parameter datasets for the forestIPM package
## Run this script once to regenerate inst/extdata/parameters/ after updating posteriors.
##
## Output: one .rds file per species in inst/extdata/parameters/
## Each file contains:
##   $draws — full set of 1000 posterior draws (iterations 2001-4000), per vital rate
##   $mean  — posterior mean across those draws, per vital rate (named numeric vector)
##
## Files are loaded on-demand inside species_model() via system.file(), not at library() time.
##
## Requirements: dplyr must be installed.
## Run data-raw/species_list.R before this script.
## This script reads from the local data path specified in _data.path.

library(dplyr)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
data_path  <- trimws(readLines("_data.path"))
output_dir <- file.path(data_path, "output_sim_processed")

vital_rate_dirs <- list(
  growth       = file.path(output_dir, "growth",       "intcpt_plot_comp_clim"),
  mort         = file.path(output_dir, "mort",         "intcpt_plot_comp_clim"),
  rec          = file.path(output_dir, "recruit",      "intcpt_plot_comp_clim"),
  sizeIngrowth = file.path(output_dir, "sizeIngrowth", "time_truc")
)

# ---------------------------------------------------------------------------
# Species list: canonical IDs and their posterior file stems
# ---------------------------------------------------------------------------
# species_list.csv is the source of truth for which species to include.
# species_id_new (strip prefix) = output file name; species_id_old = file stem.
sp_source <- read.csv(file.path(data_path, "species_id.csv"))
sp_source <- sp_source[sp_source$sp_to_analyze == TRUE, ]

sp_ids <- sub("^[0-9]+", "", sp_source$species_id_new)  # clean output IDs
stems  <- sp_source$species_id_old                       # posterior file stems

message("Processing ", length(sp_ids), " species: ", paste(sp_ids, collapse = ", "))

# ---------------------------------------------------------------------------
# Build and save one .rds per species to inst/extdata/parameters/
# ---------------------------------------------------------------------------
out_dir <- file.path("inst", "extdata", "parameters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(sp_ids)) {
  sp   <- sp_ids[[i]]
  stem <- stems[[i]]

  # Load each vital rate and retain draws 2001-4000 (1000 draws), renumbered 1-1000
  draws <- lapply(vital_rate_dirs, function(dir) {
    file.path(dir, paste0("posterior_pop_", stem, ".RDS")) |>
      readRDS() |>
      filter(.data$par != "lp__") |>
      filter(iter %in% 2001:4000) |>
      mutate(
        draw = iter - 2000L,
        iter = NULL
      ) |>
      relocate(draw, .before = par)
  })

  # Pre-compute the posterior mean across those 1000 draws for each vital rate
  mean_pars <- lapply(draws, function(vr_tbl) {
    vr_tbl |>
      group_by(.data$par) |>
      summarise(value = mean(.data$value), .groups = "drop") |>
      tidyr::pivot_wider(names_from = "par") |>
      unlist(use.names = TRUE)
  })

  saveRDS(
    list(draws = draws, mean = mean_pars),
    file     = file.path(out_dir, paste0(sp, "_pars.rds")),
    compress = "bzip2"
  )
  message("  Saved: ", sp, "_pars.rds")
}

message("\nDone. ", length(sp_ids), " species datasets written to ", out_dir, "/")
