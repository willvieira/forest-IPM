# globals.R
# Centralized package-level roxygen2 declarations:
# - Rcpp shared library registration (required for getEigenValues() to load)
# - All @importFrom declarations for symbols used in params.R without ::

#' @useDynLib forestIPM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom purrr map map_chr map_dbl map_lgl
#' @importFrom dplyr filter group_by reframe select mutate case_match group_split left_join bind_rows
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_path geom_hline theme_classic labs theme element_blank element_text
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom stringr str_replace
#' @importFrom rlang .data set_names
NULL
