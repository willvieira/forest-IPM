# Internal helper: load range values from package data
.env_ranges <- function() {
  path <- system.file("extdata", "env_scaling.RDS", package = "forestIPM")
  readRDS(path)
}

#' Scale environmental values to [0, 1]
#'
#' Scales mean annual temperature and precipitation to the [0, 1] range using
#' the observed ranges stored in the package. These scaled values are what the
#' IPM vital rate functions expect internally.
#'
#' @param MAT Numeric. Mean annual temperature in degrees Celsius.
#' @param MAP Numeric. Mean annual precipitation in mm/year.
#' @return A named list with elements \code{MAT} and \code{MAP}, each scaled
#'   to [0, 1].
#' @export
scale_env <- function(MAT, MAP) {
  rng <- .env_ranges()
  list(
    MAT = (MAT - rng$bio_01_mean[1]) / (rng$bio_01_mean[2] - rng$bio_01_mean[1]),
    MAP = (MAP - rng$bio_12_mean[1]) / (rng$bio_12_mean[2] - rng$bio_12_mean[1])
  )
}

#' Unscale environmental values from [0, 1] to natural units
#'
#' Reverses the scaling applied by \code{scale_env()}, converting scaled values
#' back to degrees Celsius and mm/year.
#'
#' @param MAT Numeric. Scaled mean annual temperature in [0, 1].
#' @param MAP Numeric. Scaled mean annual precipitation in [0, 1].
#' @return A named list with elements \code{MAT} (degrees Celsius) and
#'   \code{MAP} (mm/year).
#' @export
unscale_env <- function(MAT, MAP) {
  rng <- .env_ranges()
  list(
    MAT = MAT * (rng$bio_01_mean[2] - rng$bio_01_mean[1]) + rng$bio_01_mean[1],
    MAP = MAP * (rng$bio_12_mean[2] - rng$bio_12_mean[1]) + rng$bio_12_mean[1]
  )
}
