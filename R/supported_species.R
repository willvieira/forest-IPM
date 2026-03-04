#' Return the tibble of species supported by forestIPM
#'
#' @return A tibble with six columns: species_id, common_name, nom_commun,
#'   growth_model, surv_model, recruit_model. One row per supported species.
#' @export
supported_species <- function() {
  path <- system.file("extdata", "species_list.csv", package = "forestIPM")
  if (nchar(path) == 0) {
    # Fallback for devtools::load_all() during development
    path <- file.path("inst", "extdata", "species_list.csv")
  }
  tibble::as_tibble(utils::read.csv(path, stringsAsFactors = FALSE))
}
