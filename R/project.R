#' Project population dynamics over time
#'
#' @param mod An \code{ipm_spModel} object.
#' @param pars An \code{ipm_parameters} object.
#' @param stand An \code{ipm_stand} object. Initial stand state (t = 0).
#' @param env An \code{ipm_env} object. Climate drivers (may be time-varying).
#' @param ctrl An \code{ipm_control} object. Required.
#' @return An object of S3 class \code{"ipm_projection"} with fields:
#'   \code{$species}, \code{$years}, \code{$lambda}, \code{$stand_series}
#'   (list of \code{ipm_dist_snapshot} objects, each holding the continuous
#'   size distribution per species at that timestep), \code{$summary}.
#' @examples
#' df <- data.frame(size_mm = c(150, 200, 350),
#'                  species_id = "ABIBAL",
#'                  plot_size  = 1000)
#' s    <- stand(df)
#' mod  <- species_model(s)
#' pars <- parameters(mod, draw = "mean")
#' env  <- env_condition(MAT = 8, MAP = 1200)
#' ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
#' proj <- project(mod, pars, s, env, ctrl)
#' print(proj)
#' @export
project <- function(mod, pars, stand, env, ctrl) {
  if (!inherits(mod,   "ipm_spModel"))    cli::cli_abort("{.arg mod} must be {.cls ipm_spModel}.")
  if (!inherits(pars,  "ipm_parameters")) cli::cli_abort("{.arg pars} must be {.cls ipm_parameters}.")
  if (!inherits(stand, "ipm_stand"))      cli::cli_abort("{.arg stand} must be {.cls ipm_stand}.")
  if (!inherits(env,   "ipm_env"))        cli::cli_abort("{.arg env} must be {.cls ipm_env}.")
  if (!inherits(ctrl,  "ipm_control"))    cli::cli_abort("{.arg ctrl} must be {.cls ipm_control}.")

  # All mod species must have non-NULL parameters in pars
  missing_focal <- setdiff(mod$species, names(pars$species_params))
  null_focal    <- mod$species[purrr::map_lgl(mod$species, function(sp) {
    is.null(pars$species_params[[sp]]$fixed)
  })]
  bad_focal <- union(missing_focal, null_focal)
  if (length(bad_focal) > 0) {
    cli::cli_abort(c(
      "All {.arg mod} species must have parameters in {.arg pars}.",
      "x" = "Missing or NULL parameters for: {.val {bad_focal}}."
    ))
  }

  # Stand-only species (competitors not in mod) — handle via on_missing
  stand_only  <- setdiff(stand$species, mod$species)
  static_comp <- character(0)  # competitors whose Nvec will not update

  if (length(stand_only) > 0) {
    no_pars <- stand_only[!stand_only %in% names(pars$species_params) |
                            purrr::map_lgl(stand_only, function(sp) {
                              is.null(pars$species_params[[sp]]$fixed)
                            })]

    if (length(no_pars) > 0) {
      if (mod$on_missing == "error") {
        cli::cli_abort(c(
          "Stand contains competitor species without parameters in {.arg pars}.",
          "x" = "No parameters for: {.val {no_pars}}.",
          "i" = "Set {.code on_missing = \"drop\"} or {.code \"static\"} in {.fn species_model} to handle this."
        ))
      } else if (mod$on_missing == "drop") {
        stand_only <- setdiff(stand_only, no_pars)
      } else {  # "static"
        static_comp <- no_pars
      }
    }
  }

  warn_env_range(env, mod$species)

  all_species <- c(mod$species, stand_only)

  bin_w      <- ctrl$bin_width
  delta_time <- ctrl$delta_time

  # Build size distributions from stand data
  nvec_list <- .stand_to_nvec(stand, all_species, pars, bin_w)

  # Static competitors: save initial Nvec, never update
  static_nvec <- nvec_list[static_comp]

  # Storage
  stored_t     <- integer(0)
  lambda_store <- lapply(stats::setNames(mod$species, mod$species), function(sp) numeric(0))
  stand_series <- list()
  summary_rows <- list()

  # Projection loop
  if (ctrl$progress) {
    pb <- cli::cli_progress_bar(
      name   = "Projecting",
      total  = ctrl$years,
      format = "{cli::pb_name} {cli::pb_bar} {cli::pb_current}/{cli::pb_total} yr | ETA: {cli::pb_eta}"
    )
  }

  for (t in seq_len(ctrl$years)) {
    Temp <- if (is.function(env$.MAT_scl)) env$.MAT_scl(t) else env$.MAT_scl
    Prec <- if (is.function(env$.MAP_scl)) env$.MAP_scl(t) else env$.MAP_scl

    # Restore static competitor Nvec (unchanged)
    for (sp in static_comp) nvec_list[[sp]] <- static_nvec[[sp]]

    new_nvec_list <- nvec_list
    lambdas_t     <- stats::setNames(rep(NA_real_, length(mod$species)), mod$species)

    for (sp in mod$species) {
      sp_pars    <- pars$species_params[[sp]]$fixed
      plot_random <- if(is.null(pars$species_params[[sp]]$random_effects)) {
        c(0, 0, 0)
      }else {
        pars$species_params[[sp]]$random_effects
      }
      
      K_list <- mkKernel(
        Nvec_intra  = nvec_list[[sp]]$N_con,
        Nvec_inter  = nvec_list[[sp]]$N_het,
        delta_time  = delta_time,
        plotSize    = stand$plot_size,
        Temp        = Temp,
        Prec        = Prec,
        pars        = sp_pars,
        plot_random = plot_random
      )

      new_nvec_list[[sp]]$N_con$Nvec <- as.numeric(K_list$K %*% nvec_list[[sp]]$N_con$Nvec)

      if (ctrl$compute_lambda) {
        lambdas_t[[sp]] <- max(getEigenValues(K_list$K))
      }
    }

    nvec_list <- .update_N_het(mod$species, new_nvec_list)

    if (ctrl$progress) cli::cli_progress_update(id = pb)

    if (t %% ctrl$store_every == 0) {
      stored_t <- c(stored_t, t)

      for (sp in mod$species) {
        lambda_store[[sp]] <- c(lambda_store[[sp]], lambdas_t[[sp]])
      }

      t_stand <- .nvec_to_stand(nvec_list, mod$species, stand$plot_size)
      stand_series <- c(stand_series, list(t_stand))

      for (sp in mod$species) {
        summary_rows <- c(summary_rows, list(data.frame(
          timestep   = t,
          species_id = sp,
          lambda     = lambdas_t[[sp]],
          n_trees    = as.integer(round(sum(nvec_list[[sp]]$N_con$Nvec))),
          stringsAsFactors = FALSE
        )))
      }
    }
  }

  if (ctrl$progress) cli::cli_progress_done(id = pb)

  summary_tbl <- do.call(rbind, summary_rows)

  structure(
    list(
      species      = mod$species,
      years        = stored_t,
      lambda       = lambda_store,
      stand_series = stand_series,
      summary      = tibble::as_tibble(summary_tbl)
    ),
    class = "ipm_projection"
  )
}

#' @export
summary.ipm_projection <- function(object, ...) {
  cat("<ipm_projection> summary\n")
  cat(sprintf("  Species (%d): %s\n", length(object$species),
              paste(object$species, collapse = ", ")))
  if (length(object$years) > 0) {
    cat(sprintf("  Timesteps stored: %d (years %d to %d)\n",
                length(object$years), min(object$years), max(object$years)))
    cat("  Final timestep per species:\n")
    final_t <- max(object$years)
    final_rows <- object$summary[object$summary$timestep == final_t, ]
    for (i in seq_len(nrow(final_rows))) {
      row <- final_rows[i, ]
      lam_str <- if (is.na(row$lambda)) "NA" else sprintf("%.4f", row$lambda)
      cat(sprintf("    %s: lambda=%s, n_trees=%d\n",
                  row$species_id, lam_str, row$n_trees))
    }
  } else {
    cat("  No timesteps stored.\n")
  }
  invisible(object)
}

#' Plot an ipm_projection object
#'
#' @param x An \code{ipm_projection} object returned by \code{\link{project}}.
#' @param type Character or NULL. One of \code{"lambda"}, \code{"size_dist"},
#'   \code{"lambda_vs_n"}. If NULL (default), all three figures are rendered
#'   in sequence.
#' @param ... Additional arguments (currently unused).
#' @return \code{x}, invisibly.
#' @export
#' @examples
#' df <- data.frame(size_mm = seq(130, 600, by = 50),
#'                  species_id = "ABIBAL", plot_size = 400)
#' s    <- stand(df)
#' mod  <- species_model(s)
#' pars <- parameters(mod, draw = "mean")
#' env  <- env_condition(MAT = 8, MAP = 1200)
#' ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
#' proj <- project(mod, pars, s, env, ctrl)
#' plot(proj, type = "lambda")
plot.ipm_projection <- function(x, type = NULL, ...) {
  types <- c("lambda", "size_dist", "lambda_vs_n")
  if (!is.null(type) && !type %in% types) {
    cli::cli_abort(
      "{.arg type} must be one of {.val {types}} or NULL. Got {.val {type}}."
    )
  }
  to_plot <- if (is.null(type)) types else type
  for (t in to_plot) {
    switch(t,
      lambda      = .plot_lambda(x, ...),
      size_dist   = .plot_size_dist(x, ...),
      lambda_vs_n = .plot_lambda_vs_n(x, ...)
    )
  }
  invisible(x)
}

.plot_lambda <- function(x, ...) {
  sp_list <- x$species
  n_sp    <- length(sp_list)

  # Check if all lambda values are NA
  all_na <- all(vapply(x$lambda, function(l) all(is.na(l)), logical(1)))
  if (all_na) {
    message("Lambda not computed. Use control(compute_lambda = TRUE).")
    return(invisible(NULL))
  }

  cols <- if (n_sp <= 8) c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
                           "#ff7f00", "#a65628", "#f781bf", "#999999") else
          grDevices::rainbow(n_sp)

  for (i in seq_along(sp_list)) {
    sp  <- sp_list[[i]]
    lam <- x$lambda[[sp]]
    yrs <- seq_along(lam)
    if (i == 1) {
      plot(yrs, lam, type = "l", col = cols[i], lwd = 2,
           xlab = "Year", ylab = "Lambda", main = "Lambda over time", ...)
    } else {
      graphics::lines(yrs, lam, col = cols[i], lwd = 2)
    }
  }

  graphics::abline(h = 1, lty = 2, col = "grey50")
  if (n_sp > 1) {
    graphics::legend("topright", legend = sp_list,
                     col = cols[seq_along(sp_list)], lwd = 2, bty = "n")
  }
  invisible(NULL)
}

.plot_size_dist <- function(x, ...) {
  if (is.null(x$stand_series) || length(x$stand_series) == 0) {
    message("No population data available.")
    return(invisible(NULL))
  }

  sp_list <- x$species
  last_snap <- x$stand_series[[length(x$stand_series)]]

  for (sp in sp_list) {
    dist_df <- last_snap$distributions[[sp]]
    if (is.null(dist_df) || nrow(dist_df) == 0) next

    # Aggregate density by binned size class for barplot display
    breaks <- seq(min(dist_df$size_mm), max(dist_df$size_mm),
                  length.out = min(30L, nrow(dist_df)))
    bins    <- cut(dist_df$size_mm, breaks = breaks, include.lowest = TRUE)
    agg     <- tapply(dist_df$density, bins, mean, na.rm = TRUE)
    agg[is.na(agg)] <- 0

    graphics::barplot(agg,
      main = paste("Size distribution:", sp),
      xlab = "Size class",
      ylab = "Density",
      names.arg = rep("", length(agg)))
  }
  invisible(NULL)
}

.plot_lambda_vs_n <- function(x, ...) {
  if (is.null(x$stand_series) || length(x$stand_series) == 0) {
    return(invisible(NULL))
  }

  sp_list <- x$species
  smry    <- x$summary

  for (sp in sp_list) {
    sp_rows <- smry[smry$species_id == sp, ]
    if (nrow(sp_rows) == 0) next

    total_n <- sp_rows$n_trees
    lam     <- sp_rows$lambda

    if (all(is.na(lam))) {
      message(paste("Lambda not computed for", sp, ". Skipping lambda_vs_n plot."))
      next
    }

    plot(total_n, lam,
         xlab = "Total N",
         ylab = "Lambda",
         main = paste("Lambda vs N:", sp), ...)
    graphics::abline(h = 1, lty = 2)
  }
  invisible(NULL)
}