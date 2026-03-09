## code to prepare inst/extdata/species_list.csv for the forestIPM package
## Run this script first, before species_parameters.R.
##
## Reads species_id.csv from the local data path, filters to species with
## sp_to_analyze == TRUE, derives the clean species_id from species_id_new
## (strip leading digits), and rebuilds species_list.csv with the structure:
##   species_id, common_name, nom_commun, growth_model, surv_model, recruit_model

library(dplyr)

data_path <- trimws(readLines("_data.path"))

# ---------------------------------------------------------------------------
# Common name lookup (English and French) keyed by clean species_id
# ---------------------------------------------------------------------------
common_names <- data.frame(
  species_id  = c(
    "ABIBAL", "ACERUB", "ACESAC", "BETALL", "BETPAP",
    "CARGLA", "CARTOM", "FAGGRA", "FRAAME", "FRANIG",
    "FRAPEN", "JUNVIR", "LIQSTY", "LIRTUL", "NYSSYL",
    "PICGLA", "PICMAR", "PICRUB", "PINBAN", "PINSTR",
    "POPGRA", "POPTRE", "PRUSER", "QUEALB", "QUEPRI",
    "QUERUB", "QUESTE", "QUEVEL", "THUOCC", "TILAME",
    "TSUCAN"
  ),
  common_name = c(
    "Balsam fir", "Red maple", "Sugar maple", "Yellow birch", "White birch",
    "Pignut hickory", "Mockernut hickory", "American beech", "White ash", "Black ash",
    "Green ash", "Eastern red cedar", "Sweetgum", "Tulip tree", "Black gum",
    "White spruce", "Black spruce", "Red spruce", "Jack pine", "Eastern white pine",
    "Bigtooth aspen", "Trembling aspen", "Black cherry", "White oak", "Chestnut oak",
    "Northern red oak", "Post oak", "Black oak", "Eastern white cedar", "American basswood",
    "Eastern hemlock"
  ),
  nom_commun  = c(
    "Sapin baumier", "Erable rouge", "Erable a sucre", "Bouleau jaune", "Bouleau blanc",
    "Caryer glabre", "Caryer tomenteux", "Hetre a grandes feuilles", "Frene blanc", "Frene noir",
    "Frene rouge", "Genevrier de Virginie", "Copalme d'Amerique", "Tulipier de Virginie", "Nysse sylvestre",
    "Epicea blanc", "Epicea noir", "Epicea rouge", "Pin gris", "Pin blanc",
    "Peuplier a grandes dents", "Peuplier faux-tremble", "Cerisier tardif", "Chene blanc", "Chene chataignier",
    "Chene rouge", "Chene etoile", "Chene velutineux", "Thuya occidental", "Tilleul d'Amerique",
    "Tsuga du Canada"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Load and filter source
# ---------------------------------------------------------------------------
sp_source   <- read.csv(file.path(data_path, "species_id.csv"))
sp_filtered <- sp_source[sp_source$sp_to_analyze == TRUE, ]

sp_filtered$species_id <- sub("^[0-9]+", "", sp_filtered$species_id_new)

# ---------------------------------------------------------------------------
# Extract environmental range per species
# ---------------------------------------------------------------------------
env_summary <- readRDS(file.path(data_path, "treeData.RDS")) |>
  filter(species_id %in% sp_source$species_id_old) |>
  group_by(species_id) |>
  reframe(
    MAT_min = min(bio_01_mean, na.rm = TRUE),
    MAT_max = max(bio_01_mean, na.rm = TRUE),
    MAP_min = min(bio_12_mean, na.rm = TRUE),
    MAP_max = max(bio_12_mean, na.rm = TRUE)
  ) |>
  mutate(across(starts_with("MA"), ~ round(.x, 1))) |>
  left_join(
    sp_source |>
      select(species_id = species_id_old, species_id_new, sp_to_analyze)
  ) |>
  filter(sp_to_analyze == TRUE) |>
  mutate(species_id = sub("^[0-9]+", "", species_id_new)) |>
  select(!c(species_id_new, sp_to_analyze))

# ---------------------------------------------------------------------------
# Build output
# ---------------------------------------------------------------------------
output <- sp_filtered |>
  select(species_id) |>
  left_join(common_names) |>
  left_join(env_summary) |>
  mutate(
    growth_model = "intcpt_plot_comp_cl",
    surv_model = "intcpt_plot_comp_cl",
    recruit_model = "intcpt_plot_comp_clim"
  ) |>
  arrange(species_id)

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
write.csv(output, file.path("inst", "extdata", "species_list.csv"), row.names = FALSE)

message("Written ", nrow(output), " species to inst/extdata/species_list.csv")
