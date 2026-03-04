# globals.R
# Centralized package-level roxygen2 declarations:
# - Rcpp shared library registration (required for getEigenValues() to load)
# - All @importFrom declarations for symbols used in params.R without ::

#' @useDynLib forestIPM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom purrr map map_chr
#' @importFrom dplyr filter group_by reframe select mutate case_match group_split
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_replace
#' @importFrom rlang .data set_names
#' @importFrom magrittr %>%
NULL
