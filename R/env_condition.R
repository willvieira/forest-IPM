# new_ipm_env: low-level constructor
# MAT/MAP: natural-scale values (or functions)
# .MAT_scl/.MAP_scl: scaled [0,1] values (or functions) used by engines
new_ipm_env <- function(MAT, MAP) {
  if (is.function(MAT)) {
    MAT_scl <- function(t) scale_env(MAT(t), 0)$MAT
  } else {
    MAT_scl <- scale_env(MAT, 0)$MAT
  }
  if (is.function(MAP)) {
    MAP_scl <- function(t) scale_env(0, MAP(t))$MAP
  } else {
    MAP_scl <- scale_env(0, MAP)$MAP
  }
  structure(
    list(MAT = MAT, MAP = MAP, .MAT_scl = MAT_scl, .MAP_scl = MAP_scl),
    class = "ipm_env"
  )
}

# validate_ipm_env: checks MAT and MAP types
validate_ipm_env <- function(x) {
  if (!is.numeric(x$MAT) && !is.function(x$MAT)) {
    cli::cli_abort(c(
      '{.arg MAT} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(x$MAT)}}.'
    ))
  }
  if (!is.numeric(x$MAP) && !is.function(x$MAP)) {
    cli::cli_abort(c(
      '{.arg MAP} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(x$MAP)}}.'
    ))
  }
  x
}

#' Create an ipm_env object specifying climate drivers
#'
#' @param MAT Numeric scalar or \code{function(t)}. Mean Annual Temperature in
#'   Celsius. Use a function for time-varying climate scenarios.
#' @param MAP Numeric scalar or \code{function(t)}. Mean Annual Precipitation in
#'   mm/year.
#' @return An object of S3 class \code{"ipm_env"} with fields \code{$MAT} and
#'   \code{$MAP} (unscaled; internal scaling applied by engines).
#' @examples
#' env <- env_condition(MAT = 8, MAP = 1200)
#' print(env)
#' @export
env_condition <- function(MAT, MAP) {
  validate_ipm_env(new_ipm_env(MAT = MAT, MAP = MAP))
}

# warn_env_range: issues a warning for species whose static MAT/MAP falls
# outside the range observed in the training data.
# species: character vector of focal species IDs
warn_env_range <- function(env, species) {
  static_MAT <- if (is.numeric(env$MAT)) env$MAT else NULL
  static_MAP <- if (is.numeric(env$MAP)) env$MAP else NULL

  if (is.null(static_MAT) && is.null(static_MAP)) return(invisible(NULL))

  sp_tbl <- supported_species()
  sp_tbl <- sp_tbl[sp_tbl$species_id %in% species, ]
  if (nrow(sp_tbl) == 0L) return(invisible(NULL))

  out_msgs <- character(0)

  for (i in seq_len(nrow(sp_tbl))) {
    sp   <- sp_tbl$species_id[i]
    name <- sp_tbl$common_name[i]
    msgs <- character(0)

    if (!is.null(static_MAT)) {
      if (static_MAT < sp_tbl$MAT_min[i] || static_MAT > sp_tbl$MAT_max[i]) {
        msgs <- c(msgs, sprintf(
          "MAT = %.1f\u00b0C is outside the observed range [%.1f, %.1f]\u00b0C",
          static_MAT, sp_tbl$MAT_min[i], sp_tbl$MAT_max[i]
        ))
      }
    }

    if (!is.null(static_MAP)) {
      if (static_MAP < sp_tbl$MAP_min[i] || static_MAP > sp_tbl$MAP_max[i]) {
        msgs <- c(msgs, sprintf(
          "MAP = %.0f mm/yr is outside the observed range [%.0f, %.0f] mm/yr",
          static_MAP, sp_tbl$MAP_min[i], sp_tbl$MAP_max[i]
        ))
      }
    }

    if (length(msgs) > 0L) {
      out_msgs <- c(out_msgs, stats::setNames(
        paste0(name, " (", sp, "): ", paste(msgs, collapse = "; ")),
        "*"
      ))
    }
  }

  if (length(out_msgs) > 0L) {
    cli::cli_warn(c(
      "Environmental conditions are outside the range observed for {length(out_msgs)} species.",
      "!" = "Model extrapolation may be unreliable.",
      out_msgs
    ))
  }

  invisible(NULL)
}

#' @export
print.ipm_env <- function(x, ...) {
  mat_str <- if (is.function(x$MAT)) "<function>" else sprintf("%.1f\u00b0C", x$MAT)
  map_str <- if (is.function(x$MAP)) "<function>" else sprintf("%.0f mm/yr", x$MAP)
  cat(sprintf("<ipm_env>  MAT=%s  MAP=%s\n", mat_str, map_str))
  invisible(x)
}

#' @export
summary.ipm_env <- function(object, ...) {
  cat("<ipm_env> summary\n")
  cat(sprintf("  MAT: %s\n",
              if (is.function(object$MAT)) "time-varying function" else
                sprintf("%.2f C (static)", object$MAT)))
  cat(sprintf("  MAP: %s\n",
              if (is.function(object$MAP)) "time-varying function" else
                sprintf("%.0f mm/yr (static)", object$MAP)))
  invisible(object)
}
