#' Project population dynamics over time
#'
#' @param mod An \code{ipm_spModel} object.
#' @param pars An \code{ipm_parameters} object.
#' @param stand An \code{ipm_stand} object. Initial stand state (t = 0).
#' @param env An \code{ipm_env} object. Climate drivers (may be time-varying).
#' @param ctrl An \code{ipm_control} object. Required.
#' @return An object of S3 class \code{"ipm_projection"} with fields:
#'   \code{$species}, \code{$years}, \code{$lambda}, \code{$stand_series},
#'   \code{$summary}.
#' @export
project <- function(mod, pars, stand, env, ctrl) {
  if (!inherits(mod,   "ipm_spModel"))    cli::cli_abort("{.arg mod} must be {.cls ipm_spModel}.")
  if (!inherits(pars,  "ipm_parameters")) cli::cli_abort("{.arg pars} must be {.cls ipm_parameters}.")
  if (!inherits(stand, "ipm_stand"))      cli::cli_abort("{.arg stand} must be {.cls ipm_stand}.")
  if (!inherits(env,   "ipm_env"))        cli::cli_abort("{.arg env} must be {.cls ipm_env}.")
  if (!inherits(ctrl,  "ipm_control"))    cli::cli_abort("{.arg ctrl} must be {.cls ipm_control}.")

  # All mod species must have non-NULL parameters in pars
  missing_focal <- setdiff(mod$species, names(pars$species_params))
  null_focal    <- mod$species[vapply(mod$species, function(sp) {
    is.null(pars$species_params[[sp]]$fixed)
  }, logical(1))]
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
                            vapply(stand_only, function(sp) {
                              is.null(pars$species_params[[sp]]$fixed)
                            }, logical(1))]

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
  for (t in seq_len(ctrl$years)) {
    Temp <- if (is.function(env$MAT)) env$MAT(t) else env$MAT
    Prec <- if (is.function(env$MAP)) env$MAP(t) else env$MAP

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

#' Plot lambda trajectory from an IPM projection
#'
#' @param x An \code{ipm_projection} object.
#' @param y Ignored.
#' @param ... Additional arguments passed to \code{plot()}.
#' @export
plot.ipm_projection <- function(x, y, ...) {
  if (all(vapply(x$lambda, function(l) all(is.na(l)), logical(1)))) {
    message("No lambda values to plot (compute_lambda = FALSE in ctrl).")
    return(invisible(x))
  }

  sp_list <- x$species
  n_sp    <- length(sp_list)

  all_lambdas <- unlist(x$lambda)
  ylim <- range(all_lambdas, na.rm = TRUE)
  xlim <- range(x$years)

  plot(xlim, ylim, type = "n",
       xlab = "Timestep (years)", ylab = "Lambda",
       main = "Population growth rate over time", ...)

  cols <- if (n_sp <= 8) c("#e41a1c","#377eb8","#4daf4a","#984ea3",
                           "#ff7f00","#a65628","#f781bf","#999999") else
          grDevices::rainbow(n_sp)

  for (i in seq_along(sp_list)) {
    sp  <- sp_list[[i]]
    lam <- x$lambda[[sp]]
    graphics::lines(x$years, lam, col = cols[i], lwd = 2)
  }

  graphics::abline(h = 1, lty = 2, col = "grey50")
  graphics::legend("topright", legend = sp_list, col = cols[seq_along(sp_list)],
                   lwd = 2, bty = "n")
  invisible(x)
}