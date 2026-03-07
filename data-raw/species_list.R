## code to prepare inst/extdata/species_list.csv for the forestIPM package
## Run this script first, before species_parameters.R.
##
## Reads species_id.csv from the local data path, filters to species with
## sp_to_analyze == TRUE, derives the clean species_id from species_id_new
## (strip leading digits), and rebuilds species_list.csv with the structure:
##   species_id, common_name, nom_commun, growth_model, surv_model, recruit_model

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
# Build output
# ---------------------------------------------------------------------------
output <- merge(
  sp_filtered[, "species_id", drop = FALSE],
  common_names,
  by = "species_id",
  all.x = TRUE
)

output$growth_model  <- "intcpt_plot_comp_clim"
output$surv_model    <- "intcpt_plot_comp_clim"
output$recruit_model <- "intcpt_plot_comp_clim"

output <- output[order(output$species_id), ]

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
write.csv(output, file.path("inst", "extdata", "species_list.csv"), row.names = FALSE)

message("Written ", nrow(output), " species to inst/extdata/species_list.csv")
