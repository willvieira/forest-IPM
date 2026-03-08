# new_ipm_parameters: low-level constructor
new_ipm_parameters <- function(species_params, draw_type, draw, seed) {
  structure(
    list(
      species_params = species_params,
      draw_type      = draw_type,   # "mean" | "random" | "user_defined"
      draw           = draw,        # integer draw ID, or NULL for "mean"
      seed           = seed
    ),
    class = "ipm_parameters"
  )
}

# validate_ipm_parameters: checks draw_type and draw consistency
validate_ipm_parameters <- function(x) {
  valid_types <- c("mean", "random", "user_defined")
  if (!x$draw_type %in% valid_types) {
    cli::cli_abort(
      '{.field draw_type} must be one of {.val {valid_types}}. Got {.val {x$draw_type}}.'
    )
  }
  if (x$draw_type == "mean") {
    if (!is.null(x$draw)) {
      cli::cli_abort('{.field draw} must be NULL when {.field draw_type} is "mean".')
    }
  } else {
    if (!is.integer(x$draw) || length(x$draw) != 1L ||
        is.na(x$draw) || x$draw < 1L || x$draw > 2000L) {
      cli::cli_abort(
        '{.field draw} must be an integer between 1 and 2000 when {.field draw_type} is {.val {x$draw_type}}.'
      )
    }
  }
  x
}

#' Resolve a single parameter realization from Bayesian posteriors
#'
#' @param mod An \code{ipm_spModel} object.
#' @param draw Character \code{"mean"} or \code{"random"}, or a positive integer
#'   1-1000 selecting a specific posterior draw index.
#' @param seed Integer or NULL. Random seed. When \code{draw = "random"} and
#'   \code{seed = NULL}, a seed is auto-generated so the draw is reproducible.
#'   Retrieve it from \code{$seed} on the returned object.
#' @return An object of S3 class \code{"ipm_parameters"}.
#' @export
parameters <- function(mod, draw = "random", seed = NULL) {
  if (!inherits(mod, "ipm_spModel")) {
    cli::cli_abort(
      "{.arg mod} must be an object of class {.cls ipm_spModel}. Got {.cls {class(mod)}}."
    )
  }

  # Auto-generate seed when random draw and no seed provided
  if (isTRUE(draw == "random") && is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  if (!is.null(seed)) set.seed(seed)

  # Resolve draw_type and draw_id
  if (isTRUE(draw == "mean")) {
    draw_type <- "mean"
    draw_id   <- NULL
  } else if (isTRUE(draw == "random")) {
    draw_type <- "random"
    draw_id   <- sample(1L:2000L, 1L)
  } else {
    draw_type <- "user_defined"
    draw_id   <- as.integer(draw)
  }

  # Extract one parameter realization per species from mod$params (loaded RDS)
  species_params <- lapply(stats::setNames(mod$species, mod$species), function(sp) {
    sp_data <- mod$params[[sp]]  # list(draws = ..., mean = ...) or NULL

    if (is.null(sp_data)) {
      # RDS not shipped yet; Phase 3 will add cloud fetch
      return(list(fixed = NULL, random_effects = NULL))
    }

    if (draw_type == "mean") {
      pars_resolved <- sp_data$mean  # named list of named numeric vectors per vital rate
    } else {
      # Filter one draw row per vital rate tibble → named numeric vector
      pars_resolved <- lapply(sp_data$draws, function(vr_tbl) {
        vr_tbl |>
          dplyr::filter(.data$draw == draw_id) |>
          dplyr::select(-"draw") |>
          tidyr::pivot_wider(names_from = "par", values_from = "value") |>
          unlist(use.names = TRUE)
      })
    }

    list(fixed = pars_resolved, random_effects = NULL)
  })

  obj <- new_ipm_parameters(
    species_params = species_params,
    draw_type      = draw_type,
    draw           = draw_id,
    seed           = seed
  )
  validate_ipm_parameters(obj)
}

#' Set plot-level random effects for one or more species
#'
#' Populates the \code{$random_effects} slot inside each species entry of an
#' \code{ipm_parameters} object. Engines (\code{lambda()}, \code{project()})
#' default to \code{c(0, 0, 0)} for any species whose slot is \code{NULL}.
#'
#' @param pars An \code{ipm_parameters} object.
#' @param values A numeric vector of length 3 (growth, survival, recruitment
#'   offsets), applied to all species named in \code{species}; \emph{or} a
#'   named list keyed by species ID, each element a numeric vector of length 3.
#' @param species Character vector of species IDs. Used only when \code{values}
#'   is a numeric vector; ignored when \code{values} is a named list.
#'   If \code{NULL} and \code{values} is a numeric vector, applies to all species.
#' @return The modified \code{ipm_parameters} object.
#' @export
set_random_effects <- function(pars, values, species = NULL) {
  if (!inherits(pars, "ipm_parameters")) {
    cli::cli_abort("{.arg pars} must be an {.cls ipm_parameters} object.")
  }

  if (is.list(values)) {
    # Named list keyed by species
    for (sp in names(values)) {
      if (!sp %in% names(pars$species_params)) {
        cli::cli_abort("{.val {sp}} is not a species in {.arg pars}.")
      }
      re <- values[[sp]]
      if (!is.numeric(re) || length(re) != 3L) {
        cli::cli_abort(
          "Random effects for {.val {sp}} must be a numeric vector of length 3 \\
          (growth, survival, recruitment)."
        )
      }
      pars$species_params[[sp]]$random_effects <- re
    }
  } else {
    # Single numeric(3) applied to target species
    if (!is.numeric(values) || length(values) != 3L) {
      cli::cli_abort(
        "{.arg values} must be a numeric vector of length 3 (growth, survival, recruitment)."
      )
    }
    target <- if (is.null(species)) names(pars$species_params) else species
    for (sp in target) {
      if (!sp %in% names(pars$species_params)) {
        cli::cli_abort("{.val {sp}} is not a species in {.arg pars}.")
      }
      pars$species_params[[sp]]$random_effects <- values
    }
  }
  pars
}

#' @export
print.ipm_parameters <- function(x, ...) {
  draw_str <- switch(x$draw_type,
    mean         = "mean",
    random       = sprintf("random (id=%d)", x$draw),
    user_defined = sprintf("user_defined (id=%d)", x$draw)
  )
  seed_str <- if (!is.null(x$seed)) paste0(" seed=", x$seed) else ""
  cat(sprintf("<ipm_parameters>  draw=%s%s | %d species\n",
              draw_str, seed_str, length(x$species_params)))
  invisible(x)
}

#' @export
summary.ipm_parameters <- function(object, ...) {
  cat("<ipm_parameters> summary\n")
  cat(sprintf("  Draw type: %s\n", object$draw_type))
  if (!is.null(object$draw)) cat(sprintf("  Draw ID: %d\n", object$draw))
  cat(sprintf("  Seed: %s\n", if (!is.null(object$seed)) object$seed else "NULL"))
  cat(sprintf("  Species (%d):\n", length(object$species_params)))
  for (sp in names(object$species_params)) {
    re  <- object$species_params[[sp]]$random_effects
    re_str <- if (is.null(re)) "re=NULL" else
      sprintf("re=c(%s)", paste(round(re, 3), collapse = ","))
    cat(sprintf("    %s: %s\n", sp, re_str))
  }
  invisible(object)
}