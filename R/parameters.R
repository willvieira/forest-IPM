# new_ipm_parameters: low-level constructor
new_ipm_parameters <- function(species_params, draw_type, seed) {
  structure(
    list(species_params = species_params, draw_type = draw_type, seed = seed),
    class = "ipm_parameters"
  )
}

# validate_ipm_parameters: checks mod and draw argument values
validate_ipm_parameters <- function(x) {
  x
}

#' Resolve a single parameter realization from Bayesian posteriors
#'
#' @param mod An \code{ipm_spModel} object.
#' @param draw Character \code{"mean"} or \code{"random"}, or a positive integer
#'   1-1000 selecting a specific posterior draw index.
#' @param seed Integer or NULL. Random seed used only when \code{draw = "random"}.
#' @return An object of S3 class \code{"ipm_parameters"}.
#' @export
parameters <- function(mod, draw = "random", seed = NULL) {
  if (!inherits(mod, "ipm_spModel")) {
    cli::cli_abort(c(
      "{.arg mod} must be an object of class {.cls ipm_spModel}. Got {.cls {class(mod)}}."
    ))
  }

  # Validate draw argument
  draw_is_valid <- isTRUE(draw == "mean") ||
                   isTRUE(draw == "random") ||
                   (is.numeric(draw) && length(draw) == 1 && draw >= 1 && draw <= 1000 && draw == floor(draw))

  if (!draw_is_valid) {
    cli::cli_abort(c(
      '{.arg draw} must be {.val "mean"}, {.val "random"}, or a positive integer (1-1000). Got {.val {draw}}.'
    ))
  }

  if (is.numeric(draw) && (draw < 1 || draw > 1000)) {
    cli::cli_abort(c(
      "{.arg draw} integer index must be between 1 and 1000. Got {draw}."
    ))
  }

  # Resolve draw type string for $draw_type field
  draw_type <- if (is.numeric(draw)) as.character(as.integer(draw)) else draw

  # Set seed if random draw
  if (isTRUE(draw == "random") && !is.null(seed)) set.seed(seed)

  # Extract one parameter realization per species from mod$params
  species_params <- lapply(stats::setNames(mod$species, mod$species), function(sp) {
    sp_pars <- mod$params[[sp]]

    if (is.null(sp_pars)) {
      # Phase 2: parameters not yet loaded (local RDS missing) — return NULL placeholder
      # Phase 3 will replace this with cloud fetch
      return(list(fixed = NULL, random_effects = NULL, draw_id = NA_integer_))
    }

    # sp_pars is the named list from getPars_sp() with keys: growth, mort, rec, sizeIngrowth
    # Each element is a named numeric vector of parameter values for that vital rate.
    # For draw = "mean": getPars_sp() was already called with method = "mean" in species_model()
    # For draw = "random" or integer: select the draw from stored posterior

    # draw_id: which posterior draw was selected
    if (isTRUE(draw == "mean")) {
      draw_id <- NA_integer_
      pars_resolved <- sp_pars  # already the mean from species_model()
    } else if (isTRUE(draw == "random")) {
      draw_id <- sample(1:1000, 1)
      pars_resolved <- sp_pars  # species_model() loaded random draw; use as-is
    } else {
      # integer draw index
      draw_id <- as.integer(draw)
      pars_resolved <- sp_pars  # Phase 3: will select exact draw from Parquet
    }

    list(
      fixed          = pars_resolved,   # named list of vital rate parameter vectors
      random_effects = NULL,            # Phase 3: plot-level random effects from Parquet
      draw_id        = draw_id
    )
  })

  obj <- new_ipm_parameters(
    species_params = species_params,
    draw_type      = draw_type,
    seed           = seed
  )
  validate_ipm_parameters(obj)
}

#' @export
print.ipm_parameters <- function(x, ...) {
  seed_str <- if (!is.null(x$seed)) paste0("seed=", x$seed) else "seed=NULL"
  cat(sprintf("<ipm_parameters>  draw=%s %s | %d species\n",
              x$draw_type, seed_str, length(x$species_params)))
  invisible(x)
}

#' @export
summary.ipm_parameters <- function(object, ...) {
  cat("<ipm_parameters> summary\n")
  cat(sprintf("  Draw type: %s\n", object$draw_type))
  cat(sprintf("  Seed: %s\n", if (!is.null(object$seed)) object$seed else "NULL"))
  cat(sprintf("  Species (%d):\n", length(object$species_params)))
  for (sp in names(object$species_params)) {
    draw_id <- object$species_params[[sp]]$draw_id
    cat(sprintf("    %s: draw_id=%s\n", sp,
                if (is.na(draw_id)) "mean" else as.character(draw_id)))
  }
  invisible(object)
}
