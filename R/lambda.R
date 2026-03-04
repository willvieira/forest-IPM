#' Compute the asymptotic population growth rate (lambda) per species
#'
#' @param mod An \code{ipm_spModel} object.
#' @param pars An \code{ipm_parameters} object. Must contain parameters for all
#'   species in \code{mod}.
#' @param stand An \code{ipm_stand} object. Initial size distribution.
#' @param env An \code{ipm_env} object. Climate drivers.
#' @param ctrl An \code{ipm_control} object or NULL. Only \code{bin_width} and
#'   \code{delta_time} are used; \code{years} and \code{store_every} are ignored.
#'   If NULL, defaults are used (bin_width = 1, delta_time = 1).
#' @return An object of S3 class \code{"ipm_lambda"} — a named numeric vector
#'   with one element per species in \code{mod}.
#' @export
lambda <- function(mod, pars, stand, env, ctrl = NULL) {
  # Input validation
  if (!inherits(mod,   "ipm_spModel"))    cli::cli_abort("{.arg mod} must be {.cls ipm_spModel}.")
  if (!inherits(pars,  "ipm_parameters")) cli::cli_abort("{.arg pars} must be {.cls ipm_parameters}.")
  if (!inherits(stand, "ipm_stand"))      cli::cli_abort("{.arg stand} must be {.cls ipm_stand}.")
  if (!inherits(env,   "ipm_env"))        cli::cli_abort("{.arg env} must be {.cls ipm_env}.")

  # Cross-constructor checks
  missing_sp <- setdiff(mod$species, names(pars$species_params))
  if (length(missing_sp) > 0) {
    cli::cli_abort("pars does not contain parameters for all species in mod. Missing: {.val {missing_sp}}.")
  }
  extra_sp <- setdiff(stand$species, mod$species)
  if (length(extra_sp) > 0) {
    cli::cli_abort("stand contains species not in mod: {.val {extra_sp}}.")
  }

  if (is.null(ctrl)) ctrl <- control()
  bin_w      <- ctrl$bin_width
  delta_time <- ctrl$delta_time

  # Resolve static climate at t = 0
  Temp <- if (is.function(env$MAT)) env$MAT(0) else env$MAT
  Prec <- if (is.function(env$MAP)) env$MAP(0) else env$MAP
  plot_random <- c(0, 0, 0)  # Phase 2: no plot random effects; Phase 3 will use pars$species_params[[sp]]$random_effects

  # Build size distributions for all species
  nvec_list <- lapply(stats::setNames(mod$species, mod$species), function(sp) {
    sp_pars  <- pars$species_params[[sp]]$fixed
    n_trees  <- sum(stand$trees$species_id == sp)

    if (is.null(sp_pars)) {
      # Phase 2 fallback: parameters not loaded; use a minimal default
      lmax <- 500
      sp_pars_stub <- list(growth = c(Lmax = lmax))
      pop <- init_pop(params = sp_pars_stub, L = 127, h = bin_w, N = max(n_trees, 1))
    } else {
      pop <- init_pop(params = sp_pars, L = 127, h = bin_w, N = max(n_trees, 1))
    }

    # If species has actual observed trees, replace smooth init with observed distribution
    if (n_trees > 0) {
      sp_dbh <- stand$trees$size_mm[stand$trees$species_id == sp]
      pop    <- dbh_to_sizeDist(dbh = sp_dbh, N_intra = pop)
    }
    pop
  })

  # Compute competitor Nvec (sum of all other species) for background competition
  .make_inter <- function(focal_sp, nvec_list) {
    other_sp <- setdiff(names(nvec_list), focal_sp)
    if (length(other_sp) == 0) {
      return(nvec_list[[focal_sp]])  # single species: intra = inter
    }
    # Aggregate competitors onto focal species mesh
    focal_mesh <- nvec_list[[focal_sp]]
    inter_nvec <- rowSums(sapply(other_sp, function(sp) {
      # Project competitor Nvec onto focal mesh via linear interpolation
      stats::approx(nvec_list[[sp]]$meshpts, nvec_list[[sp]]$Nvec,
                    xout = focal_mesh$meshpts, rule = 2)$y
    }))
    focal_mesh$Nvec <- inter_nvec
    focal_mesh
  }

  # Compute lambda per species
  lambdas <- vapply(mod$species, function(sp) {
    sp_pars <- pars$species_params[[sp]]$fixed
    if (is.null(sp_pars)) return(NA_real_)

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

    max(getEigenValues(K_list$K))
  }, numeric(1))

  structure(lambdas, class = "ipm_lambda")
}

#' @export
print.ipm_lambda <- function(x, ...) {
  cat("<ipm_lambda>\n")
  for (sp in names(x)) {
    cat(sprintf("  %s: %.4f\n", sp, x[[sp]]))
  }
  invisible(x)
}

#' @export
summary.ipm_lambda <- function(object, ...) {
  cat("<ipm_lambda> summary\n")
  cat(sprintf("  Species: %d\n", length(object)))
  for (sp in names(object)) {
    status <- if (is.na(object[[sp]])) "NA (parameters not loaded)"
              else if (object[[sp]] > 1) "growing (lambda > 1)"
              else if (object[[sp]] < 1) "declining (lambda < 1)"
              else "stable (lambda = 1)"
    cat(sprintf("  %s: %.4f  [%s]\n", sp, object[[sp]], status))
  }
  invisible(object)
}
