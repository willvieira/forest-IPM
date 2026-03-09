#' Compute the asymptotic population growth rate (lambda) per species
#'
#' @param mod An \code{ipm_spModel} object.
#' @param pars An \code{ipm_parameters} object. Must contain parameters for at
#'   least one species in \code{mod}. Lambda is computed only for species present
#'   in both \code{mod} and \code{pars} — other species in \code{stand} are used
#'   as competitors but do not appear in the output.
#' @param stand An \code{ipm_stand} object. Provides size distributions for all
#'   species (focal and competitors).
#' @param env An \code{ipm_env} object. Climate drivers.
#' @param ctrl An \code{ipm_control} object or NULL. Only \code{bin_width} and
#'   \code{delta_time} are used. If NULL, defaults are used.
#' @return An object of S3 class \code{"ipm_lambda"} — a named numeric vector
#'   with one element per focal species (species with available parameters).
#' @export
lambda <- function(mod, pars, stand, env, ctrl = NULL) {
  if (!inherits(mod,   "ipm_spModel"))    cli::cli_abort("{.arg mod} must be {.cls ipm_spModel}.")
  if (!inherits(pars,  "ipm_parameters")) cli::cli_abort("{.arg pars} must be {.cls ipm_parameters}.")
  if (!inherits(stand, "ipm_stand"))      cli::cli_abort("{.arg stand} must be {.cls ipm_stand}.")
  if (!inherits(env,   "ipm_env"))        cli::cli_abort("{.arg env} must be {.cls ipm_env}.")

  # Focal species: mod species that have parameters in pars
  focal_species <- intersect(mod$species, names(pars$species_params))
  if (length(focal_species) == 0) {
    cli::cli_abort(
      "No overlap between {.arg mod} species and {.arg pars} species. \\
      {.arg pars} must contain parameters for at least one species in {.arg mod}."
    )
  }

  warn_env_range(env, focal_species)

  if (is.null(ctrl)) ctrl <- control()
  bin_w      <- ctrl$bin_width
  delta_time <- ctrl$delta_time

  Temp <- if (is.function(env$.MAT_scl)) env$.MAT_scl(0) else env$.MAT_scl
  Prec <- if (is.function(env$.MAP_scl)) env$.MAP_scl(0) else env$.MAP_scl

  # Build size distributions from stand data for ALL species in stand + focal species.
  # Observed tree sizes determine the distribution.
  all_species <- union(stand$species, focal_species)

  nvec_list <- .stand_to_nvec(stand, all_species, pars, bin_w)

  # Compute lambda only for focal species
  lambdas <- vapply(focal_species, function(sp) {
    sp_pars <- pars$species_params[[sp]]$fixed
    if (is.null(sp_pars)) {
      cli::cli_abort(
        "Parameters for focal species {.val {sp}} are NULL. \\
        Cannot compute lambda — ensure the RDS file for this species is present."
      )
    }

    # Per-species plot-level random effects; default to zero offsets if not set
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
    status <- if (object[[sp]] > 1) "growing (lambda > 1)"
              else if (object[[sp]] < 1) "declining (lambda < 1)"
              else "stable (lambda = 1)"
    cat(sprintf("  %s: %.4f  [%s]\n", sp, object[[sp]], status))
  }
  invisible(object)
}