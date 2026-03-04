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
  # Input validation
  if (!inherits(mod,   "ipm_spModel"))    cli::cli_abort("{.arg mod} must be {.cls ipm_spModel}.")
  if (!inherits(pars,  "ipm_parameters")) cli::cli_abort("{.arg pars} must be {.cls ipm_parameters}.")
  if (!inherits(stand, "ipm_stand"))      cli::cli_abort("{.arg stand} must be {.cls ipm_stand}.")
  if (!inherits(env,   "ipm_env"))        cli::cli_abort("{.arg env} must be {.cls ipm_env}.")
  if (!inherits(ctrl,  "ipm_control"))    cli::cli_abort("{.arg ctrl} must be {.cls ipm_control}.")

  # Cross-constructor checks (same as lambda())
  missing_sp <- setdiff(mod$species, names(pars$species_params))
  if (length(missing_sp) > 0) {
    cli::cli_abort("pars does not contain parameters for all species in mod. Missing: {.val {missing_sp}}.")
  }
  extra_sp <- setdiff(stand$species, mod$species)
  if (length(extra_sp) > 0) {
    cli::cli_abort("stand contains species not in mod: {.val {extra_sp}}.")
  }

  bin_w      <- ctrl$bin_width
  delta_time <- ctrl$delta_time
  plot_random <- c(0, 0, 0)  # Phase 2: no plot random effects

  # Initialize size distributions from stand
  nvec_list <- lapply(stats::setNames(mod$species, mod$species), function(sp) {
    sp_pars <- pars$species_params[[sp]]$fixed
    n_trees <- sum(stand$trees$species_id == sp)

    lmax <- if (!is.null(sp_pars) && !is.null(sp_pars$growth)) sp_pars$growth["Lmax"] else 500
    pop  <- init_pop(params = list(growth = c(Lmax = lmax)), L = 127, h = bin_w, N = max(n_trees, 1))

    if (n_trees > 0) {
      sp_dbh <- stand$trees$size_mm[stand$trees$species_id == sp]
      pop    <- dbh_to_sizeDist(dbh = sp_dbh, N_intra = pop)
    }
    pop
  })

  # Helper: build competitor Nvec by aggregating other species onto focal mesh
  .make_inter <- function(focal_sp, nvec_list) {
    other_sp <- setdiff(names(nvec_list), focal_sp)
    if (length(other_sp) == 0) return(nvec_list[[focal_sp]])
    focal_mesh <- nvec_list[[focal_sp]]
    inter_nvec <- rowSums(sapply(other_sp, function(sp) {
      stats::approx(nvec_list[[sp]]$meshpts, nvec_list[[sp]]$Nvec,
                    xout = focal_mesh$meshpts, rule = 2)$y
    }))
    focal_mesh$Nvec <- inter_nvec
    focal_mesh
  }

  # Storage
  stored_t      <- integer(0)
  lambda_store  <- lapply(stats::setNames(mod$species, mod$species), function(sp) numeric(0))
  stand_series  <- list()
  summary_rows  <- list()

  # Projection loop
  for (t in seq_len(ctrl$years)) {
    Temp <- if (is.function(env$MAT)) env$MAT(t) else env$MAT
    Prec <- if (is.function(env$MAP)) env$MAP(t) else env$MAP

    new_nvec_list <- nvec_list  # will be updated each species this step
    lambdas_t     <- numeric(length(mod$species))
    names(lambdas_t) <- mod$species

    for (sp in mod$species) {
      sp_pars <- pars$species_params[[sp]]$fixed
      if (is.null(sp_pars)) {
        lambdas_t[[sp]] <- NA_real_
        next
      }

      Nvec_intra <- nvec_list[[sp]]
      Nvec_inter <- .make_inter(sp, nvec_list)

      K_list <- mkKernel(
        Nvec_intra  = Nvec_intra,
        Nvec_inter  = Nvec_inter,
        delta_time  = delta_time,
        plotSize    = stand$plot_size,
        Temp        = Temp,
        Prec        = Prec,
        pars        = sp_pars,
        plot_random = plot_random
      )

      # Project one timestep: new Nvec = K %*% old Nvec
      new_nvec                   <- Nvec_intra
      new_nvec$Nvec              <- as.numeric(K_list$K %*% Nvec_intra$Nvec)
      new_nvec_list[[sp]]        <- new_nvec

      lambdas_t[[sp]] <- max(getEigenValues(K_list$K))
    }

    # Update nvec_list for next timestep (fully coupled for multi-species)
    nvec_list <- new_nvec_list

    # Store this timestep?
    if (t %% ctrl$store_every == 0) {
      stored_t <- c(stored_t, t)

      for (sp in mod$species) {
        lambda_store[[sp]] <- c(lambda_store[[sp]], lambdas_t[[sp]])
      }

      # Reconstruct ipm_stand for this timestep
      t_stand <- .nvec_to_stand(nvec_list, mod$species, stand$plot_size)
      stand_series <- c(stand_series, list(t_stand))

      # Summary rows
      for (sp in mod$species) {
        summary_rows <- c(summary_rows, list(data.frame(
          timestep   = t,
          species_id = sp,
          lambda     = lambdas_t[[sp]],
          n_trees    = as.integer(round(sum(nvec_list[[sp]]$Nvec))),
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

# Internal: convert nvec_list back to an ipm_stand for storage in $stand_series
.nvec_to_stand <- function(nvec_list, species, plot_size) {
  tree_rows <- lapply(species, function(sp) {
    nvec  <- nvec_list[[sp]]
    n_ind <- round(sum(nvec$Nvec))
    if (n_ind < 1) return(NULL)
    # Represent each "individual" by its mesh midpoint (aggregate representation)
    sizes <- rep(nvec$meshpts, times = round(nvec$Nvec))
    data.frame(
      size_mm    = sizes,
      species_id = sp,
      plot_size  = plot_size,
      stringsAsFactors = FALSE
    )
  })
  tree_rows <- Filter(Negate(is.null), tree_rows)
  if (length(tree_rows) == 0) {
    trees_df <- data.frame(size_mm = numeric(0), species_id = character(0),
                           plot_size = numeric(0), stringsAsFactors = FALSE)
  } else {
    trees_df <- do.call(rbind, tree_rows)
  }
  structure(
    list(
      trees     = trees_df,
      species   = species,
      plot_size = plot_size
    ),
    class = "ipm_stand"
  )
}

#' @export
print.ipm_projection <- function(x, ...) {
  store_every <- if (length(x$years) > 1) x$years[2] - x$years[1] else 1
  cat(sprintf("<ipm_projection>  %d species | %d years | store_every=%d\n",
              length(x$species),
              if (length(x$years) > 0) max(x$years) else 0,
              store_every))
  invisible(x)
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
      cat(sprintf("    %s: lambda=%.4f, n_trees=%d\n",
                  row$species_id, row$lambda, row$n_trees))
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
  sp_list <- x$species
  n_sp    <- length(sp_list)

  # Set up plot area
  all_lambdas <- unlist(x$lambda)
  ylim <- range(all_lambdas, na.rm = TRUE)
  xlim <- range(x$years)

  plot(xlim, ylim, type = "n",
       xlab = "Timestep (years)", ylab = "Lambda",
       main = "Population growth rate over time", ...)

  # One line per species
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
