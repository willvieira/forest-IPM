# new_ipm_spModel: low-level constructor
new_ipm_spModel <- function(species, params, on_missing) {
  structure(
    list(species = species, params = params, on_missing = on_missing),
    class = "ipm_spModel"
  )
}

# validate_ipm_spModel: checks species against supported_species()
validate_ipm_spModel <- function(x, on_missing) {
  supported <- supported_species()$species_id
  bad_ids <- setdiff(x$species, supported)

  if (length(bad_ids) > 0) {
    if (on_missing == "error") {
      # Closest-match suggestion using stringdist
      closest <- vapply(bad_ids, function(id) {
        idx <- stringdist::amatch(id, supported, maxDist = Inf)
        if (is.na(idx)) "<no match>" else supported[idx]
      }, character(1))
      cli::cli_abort(c(
        "{.arg x} contains species IDs not found in {.run supported_species()}.",
        "x" = "Unknown: {.val {bad_ids}}.",
        "i" = "Did you mean {.val {closest}}?",
        "i" = "Run {.run supported_species()} to see all valid IDs."
      ))
    } else if (on_missing == "drop") {
      x$species <- setdiff(x$species, bad_ids)
      x$params  <- x$params[x$species]
    }
    # on_missing == "static": keep all species in $species; no error
  }
  x
}

#' Create an ipm_spModel object defining which species to model
#'
#' @param x An \code{ipm_stand} object or a character vector of species IDs.
#' @param on_missing Character. Behavior when a species ID is not in
#'   \code{supported_species()}. One of \code{"error"} (default),
#'   \code{"drop"}, or \code{"static"}.
#' @return An object of S3 class \code{"ipm_spModel"} with fields
#'   \code{$species}, \code{$params}, and \code{$on_missing}.
#' @export
species_model <- function(x, on_missing = "error") {
  # Resolve species IDs from stand or character vector
  if (inherits(x, "ipm_stand")) {
    sp_ids <- x$species
  } else {
    sp_ids <- as.character(x)
  }

  if (!on_missing %in% c("error", "drop", "static")) {
    cli::cli_abort(c(
      "{.arg on_missing} must be one of {.val error}, {.val drop}, or {.val static}.",
      "x" = "Got {.val {on_missing}}."
    ))
  }

  # Load parameters eagerly for each species using the internal getPars_sp()
  # In Phase 3 this will be replaced by cloud Parquet fetching.
  # For now: attempt to load; if path does not exist, store NULL (Phase 3 will fill this).
  params <- lapply(stats::setNames(sp_ids, sp_ids), function(sp) {
    tryCatch(
      getPars_sp(sp = sp, method = "random"),  # load full posterior draws
      error = function(e) NULL  # Phase 2: local RDS may not exist yet; Phase 3 adds cloud fetch
    )
  })

  obj <- new_ipm_spModel(
    species    = sp_ids,
    params     = params,
    on_missing = on_missing
  )
  validate_ipm_spModel(obj, on_missing)
}

#' @export
print.ipm_spModel <- function(x, ...) {
  sp_str <- paste(x$species, collapse = ", ")
  if (nchar(sp_str) > 60) sp_str <- paste0(substr(sp_str, 1, 57), "...")
  cat(sprintf("<ipm_spModel>  %d species: %s\n", length(x$species), sp_str))
  invisible(x)
}

#' @export
summary.ipm_spModel <- function(object, ...) {
  cat("<ipm_spModel> summary\n")
  cat(sprintf("  Species (%d): %s\n", length(object$species),
              paste(object$species, collapse = ", ")))
  cat(sprintf("  on_missing: %s\n", object$on_missing))
  cat(sprintf("  Parameters loaded: %d / %d\n",
              sum(!vapply(object$params, is.null, logical(1))),
              length(object$species)))
  invisible(object)
}
