# new_ipm_stand: low-level constructor — no validation
new_ipm_stand <- function(trees, species, plot_size) {
  structure(
    list(trees = trees, species = species, plot_size = plot_size),
    class = "ipm_stand"
  )
}

# validate_ipm_stand: thorough value checks via cli_abort()
validate_ipm_stand <- function(x) {
  trees <- x$trees

  # Size column: accept "size_mm" (already normalized) — must be present
  if (!"size_mm" %in% names(trees)) {
    cli::cli_abort(c(
      "data must contain a column named {.field size} or {.field dbh} (size in mm).",
      "i" = "Provide individual tree sizes in mm."
    ))
  }
  if (any(trees[["size_mm"]] < 127, na.rm = TRUE)) {
    n <- sum(trees[["size_mm"]] < 127, na.rm = TRUE)
    cli::cli_abort(c(
      "All tree sizes must be >= 127 mm. Found {n} value{?s} below 127 mm.",
      "i" = "The 127 mm threshold matches the minimum DBH for inventory data."
    ))
  }

  # Species column: must be present after normalization
  if (!"species_id" %in% names(trees)) {
    cli::cli_abort(c(
      "data must contain a column named {.field sp}, {.field species} or {.field species_id}.",
      "i" = "Provide species identifiers as a character column."
    ))
  }

  # plot_size column
  if (!"plot_size" %in% names(trees)) {
    cli::cli_abort(c(
      "data must contain a {.field plot_size} column (plot area in m2).",
      "i" = "Add a plot_size column with a single positive numeric value for all rows."
    ))
  }

  x
}

# stand: user-facing helper — normalizes column names, then constructs + validates
#' Create an ipm_stand object representing a forest plot
#'
#' @param data A data.frame or tibble with individual-level tree records.
#'   Required columns: a size column (named \code{size} or \code{dbh}, in mm,
#'   >= 127 mm), a species column (named \code{sp}, \code{species}, or
#'   \code{species_id}), and a \code{plot_size} column (plot area in m2,
#'   same value for all rows).
#' @return An object of S3 class \code{"ipm_stand"} with fields \code{$trees}
#'   (data.frame with standardized columns \code{size_mm} and \code{species_id}),
#'   \code{$species} (character vector of unique species IDs), and
#'   \code{$plot_size} (numeric scalar).
#' @export
stand <- function(data) {
  data <- as.data.frame(data)

  # Normalize size column name to size_mm
  if ("dbh" %in% names(data)) {
    names(data)[names(data) == "dbh"] <- "size_mm"
  } else if ("size" %in% names(data)) {
    names(data)[names(data) == "size"] <- "size_mm"
  }

  # Normalize species column name to species_id
  for (col in c("sp", "species")) {
    if (col %in% names(data)) {
      names(data)[names(data) == col] <- "species_id"
      break
    }
  }

  plot_sz <- if ("plot_size" %in% names(data)) data$plot_size[[1]] else NA_real_

  # Keep only standardized columns in $trees
  keep_cols <- intersect(c("size_mm", "species_id", "plot_size"), names(data))
  trees_df <- data[, keep_cols, drop = FALSE]

  obj <- new_ipm_stand(
    trees     = trees_df,
    species   = unique(as.character(trees_df$species_id)),
    plot_size = plot_sz
  )
  validate_ipm_stand(obj)
}

#' @export
print.ipm_stand <- function(x, ...) {
  n_trees   <- nrow(x$trees)
  n_species <- length(x$species)
  cat(sprintf("<ipm_stand>  %d tree%s | %d species | %.0f m2 plot\n",
              n_trees, if (n_trees == 1) "" else "s",
              n_species, x$plot_size))
  invisible(x)
}

#' @export
summary.ipm_stand <- function(object, ...) {
  cat("<ipm_stand> summary\n")
  cat(sprintf("  Plot size: %.0f m2\n", object$plot_size))
  cat(sprintf("  Species (%d):\n", length(object$species)))
  for (sp in object$species) {
    sp_trees <- object$trees[object$trees$species_id == sp, ]
    cat(sprintf("    %s: %d trees, size range %.0f-%.0f mm\n",
                sp, nrow(sp_trees),
                min(sp_trees$size_mm), max(sp_trees$size_mm)))
  }
  invisible(object)
}
