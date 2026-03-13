# tests/testthat/test-workflow.R
# Full workflow integration tests: stand -> species_model -> parameters -> lambda -> project
# Species: ABIBAL (bundled RDS in inst/extdata/parameters/ — no network needed)
# All projections use years = 5 and progress = FALSE for speed

# ---------------------------------------------------------------------------
# Shared fixture helper
# ---------------------------------------------------------------------------
make_abibal_stand <- function(n_trees = 10) {
  stand(data.frame(
    size_mm    = seq(130, 600, length.out = n_trees),
    species_id = "ABIBAL",
    plot_size  = 1000
  ))
}

# ---------------------------------------------------------------------------
# Constructor chain
# ---------------------------------------------------------------------------
test_that("stand() returns ipm_stand with correct structure", {
  s <- make_abibal_stand()
  expect_s3_class(s, "ipm_stand")
  expect_true("size_mm" %in% names(s$trees))
  expect_true("species_id" %in% names(s$trees))
})

test_that("species_model() returns ipm_spModel for valid stand", {
  s   <- make_abibal_stand()
  mod <- species_model(s)
  expect_s3_class(mod, "ipm_spModel")
})

test_that("parameters() returns ipm_parameters with draw = 'mean' (offline)", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  expect_s3_class(pars, "ipm_parameters")
})

test_that("env_condition() returns ipm_env", {
  env <- env_condition(MAT = 8, MAP = 1200)
  expect_s3_class(env, "ipm_env")
})

test_that("control() returns ipm_control", {
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  expect_s3_class(ctrl, "ipm_control")
})

# ---------------------------------------------------------------------------
# lambda engine
# ---------------------------------------------------------------------------
test_that("lambda() returns ipm_lambda with finite values", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  lam  <- lambda(mod, pars, s, env, ctrl)
  expect_s3_class(lam, "ipm_lambda")
  lam_vals <- lam[["ABIBAL"]]
  expect_true(all(is.finite(lam_vals)))
  expect_true(all(lam_vals > 0))
})

# ---------------------------------------------------------------------------
# project engine
# ---------------------------------------------------------------------------
test_that("project() returns ipm_projection with lambda and stand_series fields", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_s3_class(proj, "ipm_projection")
  expect_true(!is.null(proj$lambda))
  expect_true(!is.null(proj$stand_series))
  lam_vals <- proj$lambda[["ABIBAL"]]
  expect_true(all(is.finite(lam_vals)))
})

test_that("project() with compute_lambda = FALSE returns NA lambda", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = FALSE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_s3_class(proj, "ipm_projection")
  lam_vals <- proj$lambda[["ABIBAL"]]
  expect_true(all(is.na(lam_vals)))
})

# ---------------------------------------------------------------------------
# plot.ipm_projection()
# ---------------------------------------------------------------------------
test_that("plot(proj, type = 'lambda') renders without error", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_no_error({
    png(tempfile())
    plot(proj, type = "lambda")
    dev.off()
  })
})

test_that("plot(proj, type = 'size_dist') renders without error", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_no_error({
    png(tempfile())
    plot(proj, type = "size_dist")
    dev.off()
  })
})

test_that("plot(proj, type = 'lambda_vs_n') renders without error", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_no_error({
    png(tempfile())
    plot(proj, type = "lambda_vs_n")
    dev.off()
  })
})

test_that("plot(proj) renders all three figures without error", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_no_error({
    png(tempfile())
    plot(proj)
    dev.off()
  })
})

test_that("plot(proj, type = 'invalid') raises error", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_error(plot(proj, type = "invalid"))
})

# ---------------------------------------------------------------------------
# supported_species()
# ---------------------------------------------------------------------------
test_that("supported_species() returns a data frame with species_id column", {
  sp <- supported_species()
  expect_s3_class(sp, "data.frame")
  expect_true("species_id" %in% names(sp))
  expect_true("ABIBAL" %in% sp$species_id)
})

# ---------------------------------------------------------------------------
# stand() input validation
# ---------------------------------------------------------------------------
test_that("stand() rejects trees below minimum DBH threshold (127mm)", {
  df_small <- data.frame(size_mm = c(100, 120), species_id = "ABIBAL", plot_size = 400)
  expect_error(stand(df_small))
})
