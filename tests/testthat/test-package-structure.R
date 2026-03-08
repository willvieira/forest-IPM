# Tests for Phase 02-03: Package infrastructure (DESCRIPTION, NAMESPACE, globals.R)

pkg <- pkgload::pkg_path()

test_that("DESCRIPTION has Package: forestIPM", {
  desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
  expect_equal(as.character(desc[, "Package"]), "forestIPM")
})

test_that("DESCRIPTION Imports includes all required runtime packages", {
  desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
  imports <- desc[, "Imports"]
  required <- c("Rcpp", "cli", "dplyr", "tidyr", "purrr", "stringr",
                 "stringdist", "rlang", "magrittr", "tibble", "truncnorm")
  for (pkg_name in required) {
    expect_true(grepl(pkg_name, imports), label = paste("Imports contains", pkg_name))
  }
})

test_that("DESCRIPTION LinkingTo includes Rcpp and RcppEigen", {
  desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
  linking <- desc[, "LinkingTo"]
  expect_true(grepl("Rcpp", linking))
  expect_true(grepl("RcppEigen", linking))
})

test_that("NAMESPACE has useDynLib(forestIPM)", {
  ns_lines <- readLines(file.path(pkg, "NAMESPACE"))
  expect_true(any(grepl("useDynLib\\(forestIPM", ns_lines)))
})

test_that("NAMESPACE has at least 10 importFrom declarations", {
  ns_lines <- readLines(file.path(pkg, "NAMESPACE"))
  import_lines <- grep("^importFrom", ns_lines, value = TRUE)
  expect_gte(length(import_lines), 10L)
})

test_that("R/globals.R has centralized @useDynLib forestIPM tag", {
  src_lines <- readLines(file.path(pkg, "R/globals.R"))
  expect_true(any(grepl("@useDynLib forestIPM", src_lines)))
})

test_that("R/globals.R has @importFrom Rcpp sourceCpp", {
  src_lines <- readLines(file.path(pkg, "R/globals.R"))
  expect_true(any(grepl("@importFrom Rcpp sourceCpp", src_lines)))
})

test_that(".Rbuildignore excludes simulations/ directory", {
  rb_lines <- readLines(file.path(pkg, ".Rbuildignore"))
  expect_true(any(grepl("simulations", rb_lines)))
})
