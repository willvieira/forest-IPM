# tests/testthat/test-regression-baselines.R
# Regression guard: any optimization must produce output identical to this baseline.
# Tolerance 1e-10 catches floating-point drift from refactoring; must fail for wrong math.

test_that("lambda() output matches pre-optimization baseline for QUERUB", {
  skip_on_ci()
  skip_if_not(file.exists(testthat::test_path("fixtures/lambda_baseline.rds")),
              "Baseline fixture not found — run Task 1 setup to generate it")

  df_single <- data.frame(
    size_mm    = seq(130, 600, by = 50),
    species_id = "QUERUB",
    plot_size  = 1000
  )
  s    <- stand(df_single)
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 10, MAP = 1100)
  ctrl <- control(years = 1)

  result   <- lambda(mod, pars, s, env, ctrl)
  baseline <- readRDS(testthat::test_path("fixtures/lambda_baseline.rds"))

  expect_equal(result, baseline, tolerance = 1e-10)
})

test_that("project() lambda trajectory matches pre-optimization baseline for QUERUB", {
  skip_on_ci()
  skip_if_not(file.exists(testthat::test_path("fixtures/project_baseline.rds")),
              "Baseline fixture not found — run Task 1 setup to generate it")

  df_single <- data.frame(
    size_mm    = seq(130, 600, by = 50),
    species_id = "QUERUB",
    plot_size  = 1000
  )
  s    <- stand(df_single)
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 10, MAP = 1100)
  ctrl <- control(years = 10, compute_lambda = TRUE)

  result   <- project(mod, pars, s, env, ctrl)$lambda[["QUERUB"]]
  baseline <- readRDS(testthat::test_path("fixtures/project_baseline.rds"))

  expect_equal(result, baseline, tolerance = 1e-10)
})
