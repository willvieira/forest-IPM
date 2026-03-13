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
#' @examples
#' df <- data.frame(size_mm = c(150, 200, 350),
#'                  species_id = "ABIBAL",
#'                  plot_size  = 1000)
#' s <- stand(df)
#' print(s)
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

# Internal: build nvec_list (species -> list(N_con, N_het)) from a stand object.
# Mirrors the approach in lambda(): uses dbh_to_sizeDist to bin observed tree
# sizes into the species mesh. N_het pools all trees that are not conspecific.
.stand_to_nvec <- function(stand, species, pars, bin_w) {
  lapply(stats::setNames(species, species), function(sp) {
    sp_pars     <- pars$species_params[[sp]]$fixed
    sp_trees    <- stand$trees$size_mm[stand$trees$species_id == sp]
    other_trees <- stand$trees$size_mm[stand$trees$species_id != sp]

    lmax <- if (!is.null(sp_pars) && !is.null(sp_pars$growth)) {
      ceiling(sp_pars$growth[["Lmax"]])
    } else if (length(sp_trees) > 0) {
      ceiling(max(sp_trees)) + bin_w * 5L
    } else {
      500L
    }

    m      <- max(1L, as.integer(ceiling((lmax - 127) / bin_w)))
    msh    <- 127 + ((seq_len(m)) - 0.5) * bin_w
    N_zero <- list(meshpts = msh, Nvec = rep(0.0, m), h = bin_w)

    N_con <- if (length(sp_trees)    > 0) dbh_to_sizeDist(sp_trees,    N_zero) else N_zero
    N_het <- if (length(other_trees) > 0) dbh_to_sizeDist(other_trees, N_zero) else N_zero

    list(N_con = N_con, N_het = N_het)
  })
}

# Internal: snapshot nvec_list as continuous size distributions for $stand_series.
# Returns an ipm_dist_snapshot: a named list (by species) of data frames with
# columns size_mm (mesh midpoints) and density (continuous Nvec, fractional individuals).
.nvec_to_stand <- function(nvec_list, species, plot_size) {
  dists <- stats::setNames(lapply(species, function(sp) {
    nvec <- nvec_list[[sp]]$N_con
    data.frame(size_mm = nvec$meshpts, density = nvec$Nvec, stringsAsFactors = FALSE)
  }), species)
  structure(
    list(distributions = dists, species = species, plot_size = plot_size),
    class = "ipm_dist_snapshot"
  )
}

#' @export
print.ipm_projection <- function(x, ...) {
  store_every <- if (length(x$years) > 1) x$years[2] - x$years[1] else 1
  cond <- x$conditions
  draw_str <- if (!is.null(cond)) {
    switch(cond$draw_type,
      mean         = "mean",
      random       = sprintf("random (id=%d, seed=%d)", cond$draw, cond$seed),
      user_defined = sprintf("draw=%d", cond$draw),
      "unknown"
    )
  } else "unknown"
  cat(sprintf("<ipm_projection>  %d species | %d years | store_every=%d | draw=%s\n",
              length(x$species),
              if (length(x$years) > 0) max(x$years) else 0,
              store_every,
              draw_str))
  invisible(x)
}

# Internal: recompute N_het for all focal_species from the current N_con Nvec
# objects. Bins each competitor's distribution onto the focal mesh using
# findInterval on the focal bin boundaries — sizes beyond the focal lmax are
# accumulated into the last bin so no individuals are lost.
# Returns the full nvec_list with N_het updated for focal_species.
.update_N_het <- function(focal_species, nvec_list) {
  updated <- lapply(stats::setNames(focal_species, focal_species), function(sp) {
    other_sp <- setdiff(names(nvec_list), sp)
    N_con    <- nvec_list[[sp]]$N_con
    m        <- length(N_con$meshpts)
    breaks   <- N_con$meshpts[1L] - N_con$h / 2 + (0:m) * N_con$h

    N_het      <- N_con
    N_het$Nvec <- if (length(other_sp) == 0) {
      rep(0.0, m)
    } else {
      Reduce(`+`, lapply(other_sp, function(s) {
        comp    <- nvec_list[[s]]$N_con
        bins    <- pmax(1L, pmin(findInterval(comp$meshpts, breaks), m))
        out     <- numeric(m)
        for (i in seq_along(bins)) out[bins[i]] <- out[bins[i]] + comp$Nvec[i]
        out
      }))
    }
    modifyList(nvec_list[[sp]], list(N_het = N_het))
  })
  modifyList(nvec_list, updated)
}
