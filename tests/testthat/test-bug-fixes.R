# Tests for Phase 02-01: Bug fixes (BUG-01, BUG-02, BUG-03)

pkg <- pkgload::pkg_path()

# BUG-01: General EigenSolver for non-symmetric matrices
test_that("BUG-01: getEigenValues() matches base R eigen() within 1e-10 on asymmetric matrix", {
  K <- matrix(c(0.8, 0.1, 0.05, 0.9), nrow = 2)
  rcpp_lambda <- max(getEigenValues(K))
  base_lambda <- Re(eigen(K)$values[1])
  expect_true(isTRUE(all.equal(rcpp_lambda, base_lambda, tolerance = 1e-10)))
})

test_that("BUG-01: getEigenValues() returns a real numeric vector (not complex)", {
  K <- matrix(c(0.8, 0.1, 0.05, 0.9), nrow = 2)
  result <- getEigenValues(K)
  expect_type(result, "double")
  expect_false(is.complex(result))
  expect_equal(length(result), 2L)
})

test_that("BUG-01: SelfAdjointEigenSolver is absent from src/eigen.cpp", {
  src_lines <- readLines(file.path(pkg, "src/eigen.cpp"))
  expect_equal(sum(grepl("SelfAdjointEigenSolver", src_lines)), 0)
  expect_true(any(grepl("EigenSolver", src_lines)))
})

# BUG-02: init_pop() no longer uses exists('fct') global env check
test_that("BUG-02: exists('fct') is absent from R/kernel.R", {
  src_lines <- readLines(file.path(pkg, "R/kernel.R"))
  expect_equal(sum(grepl("exists\\('fct'\\)", src_lines)), 0)
})

test_that("BUG-02: fct <- 1 initializer is present in R/kernel.R before the while-loop", {
  src_lines <- readLines(file.path(pkg, "R/kernel.R"))
  expect_true(any(grepl("fct <- 1", src_lines)))
})

# BUG-03: purrr::as_vector() removed from R/params.R
test_that("BUG-03: zero as_vector() occurrences remain in R/params.R", {
  src_lines <- readLines(file.path(pkg, "R/params.R"))
  expect_equal(sum(grepl("as_vector", src_lines)), 0)
})

test_that("BUG-03: unlist(use.names = TRUE) is present in R/params.R", {
  src_lines <- readLines(file.path(pkg, "R/params.R"))
  expect_true(any(grepl("unlist\\(use\\.names\\s*=\\s*TRUE\\)", src_lines)))
})
