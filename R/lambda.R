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

  if (is.null(ctrl)) ctrl <- control()
  bin_w      <- ctrl$bin_width
  delta_time <- ctrl$delta_time

  Temp <- if (is.function(env$MAT)) env$MAT(0) else env$MAT
  Prec <- if (is.function(env$MAP)) env$MAP(0) else env$MAP

  # Build size distributions from stand data for ALL species in stand + focal species.
  # Observed tree sizes determine the distribution; init_pop() is not used.
  all_species <- union(stand$species, focal_species)

  nvec_list <- lapply(stats::setNames(all_species, all_species), function(sp) {
    sp_trees <- stand$trees$size_mm[stand$trees$species_id == sp]
    sp_pars  <- pars$species_params[[sp]]$fixed

    # Upper mesh bound: Lmax from params if available, else from observed sizes
    lmax <- if (!is.null(sp_pars) && !is.null(sp_pars$growth)) {
      round(sp_pars$growth["Lmax"], 0)
    } else if (length(sp_trees) > 0) {
      ceiling(max(sp_trees)) + bin_w * 5L
    } else {
      500L
    }

    m   <- max(1L, as.integer(ceiling((lmax - 127) / bin_w)))
    msh <- 127 + ((seq_len(m)) - 0.5) * bin_w
    out <- list(meshpts = msh, Nvec = rep(0.0, m), h = bin_w)

    if (length(sp_trees) > 0) out <- dbh_to_sizeDist(dbh = sp_trees, N_intra = out)
    out
  })

  # Aggregate competitors onto focal mesh via linear interpolation
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
    status <- if (object[[sp]] > 1) "growing (lambda > 1)"
              else if (object[[sp]] < 1) "declining (lambda < 1)"
              else "stable (lambda = 1)"
    cat(sprintf("  %s: %.4f  [%s]\n", sp, object[[sp]], status))
  }
  invisible(object)
}