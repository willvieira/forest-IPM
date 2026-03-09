# new_ipm_control: low-level constructor
new_ipm_control <- function(years, delta_time, store_every, bin_width, compute_lambda, progress) {
  structure(
    list(years = years, delta_time = delta_time,
         store_every = store_every, bin_width = bin_width,
         compute_lambda = compute_lambda, progress = progress),
    class = "ipm_control"
  )
}

# validate_ipm_control: checks all arguments
validate_ipm_control <- function(x) {
  is_pos_int <- function(v) is.numeric(v) && length(v) == 1 && v > 0 && v == floor(v)
  is_pos_num <- function(v) is.numeric(v) && length(v) == 1 && v > 0

  if (!is_pos_int(x$years)) {
    cli::cli_abort(
      "{.arg years} must be a positive integer. Got {.val {x$years}}."
    )
  }
  if (!is_pos_num(x$delta_time)) {
    cli::cli_abort(
      "{.arg delta_time} must be a positive numeric. Got {.val {x$delta_time}}."
    )
  }
  if (!is_pos_int(x$store_every)) {
    cli::cli_abort(
      "{.arg store_every} must be a positive integer. Got {.val {x$store_every}}."
    )
  }
  if (!is_pos_int(x$bin_width)) {
    cli::cli_abort(
      "{.arg bin_width} must be a positive integer. Got {.val {x$bin_width}}."
    )
  }
  if (!is.logical(x$compute_lambda) || length(x$compute_lambda) != 1 || is.na(x$compute_lambda)) {
    cli::cli_abort(
      "{.arg compute_lambda} must be TRUE or FALSE."
    )
  }
  if (!is.logical(x$progress) || length(x$progress) != 1 || is.na(x$progress)) {
    cli::cli_abort(
      "{.arg progress} must be TRUE or FALSE."
    )
  }
  x
}

#' Configure IPM projection settings
#'
#' @param years Positive integer. Number of simulation timesteps. Default 100.
#' @param delta_time Positive numeric. Duration of each timestep in years. Default 1.
#' @param store_every Positive integer. Store stand state every N timesteps. Default 1.
#' @param bin_width Positive integer. Bin width for IPM kernel discretization. Default 1.
#' @param compute_lambda Logical. Whether to compute the asymptotic lambda at each
#'   timestep via eigendecomposition. Set to FALSE to skip (faster projections when
#'   only population structure is needed). Default FALSE
#' @param progress Logical. Whether to display a progress bar during projection.
#'   Default TRUE.
#' @return An object of S3 class \code{"ipm_control"}.
#' @export
control <- function(years = 100, delta_time = 1, store_every = 1, bin_width = 1,
                    compute_lambda = FALSE, progress = TRUE) {
  validate_ipm_control(
    new_ipm_control(
      years          = years,
      delta_time     = delta_time,
      store_every    = store_every,
      bin_width      = bin_width,
      compute_lambda = compute_lambda,
      progress       = progress
    )
  )
}

#' @export
print.ipm_control <- function(x, ...) {
  cat(sprintf("<ipm_control>  %d years | dt=%.1f | store_every=%d | bin_width=%d | lambda=%s | progress=%s\n",
              x$years, x$delta_time, x$store_every, x$bin_width,
              if (x$compute_lambda) "yes" else "no",
              if (x$progress) "yes" else "no"))
  invisible(x)
}

#' @export
summary.ipm_control <- function(object, ...) {
  cat("<ipm_control> summary\n")
  stored_steps <- ceiling(object$years / object$store_every)
  cat(sprintf("  Total years: %d\n", object$years))
  cat(sprintf("  Timestep (delta_time): %.2f years\n", object$delta_time))
  cat(sprintf("  Store every: %d timestep(s) => %d stored states\n",
              object$store_every, stored_steps))
  cat(sprintf("  Kernel bin width: %d mm\n", object$bin_width))
  cat(sprintf("  Compute lambda: %s\n", if (object$compute_lambda) "yes" else "no"))
  cat(sprintf("  Progress bar:   %s\n", if (object$progress) "yes" else "no"))
  invisible(object)
}