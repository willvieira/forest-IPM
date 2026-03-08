# Tests for Phase 02-02: API constructors and engines

# --- env_condition() ---

test_that("env_condition() accepts numeric scalars and returns ipm_env", {
  env <- env_condition(MAT = 8, MAP = 1200)
  expect_s3_class(env, "ipm_env")
  expect_equal(env$MAT, 8)
  expect_equal(env$MAP, 1200)
})

test_that("env_condition() accepts function(t) for MAT", {
  env <- env_condition(MAT = function(t) 8 + 0.02 * t, MAP = 1200)
  expect_s3_class(env, "ipm_env")
  expect_true(is.function(env$MAT))
})

# --- control() ---

test_that("control() returns ipm_control with correct defaults", {
  ctrl <- control()
  expect_s3_class(ctrl, "ipm_control")
  expect_equal(ctrl$years, 100)
  expect_equal(ctrl$delta_time, 1)
  expect_equal(ctrl$store_every, 1)
  expect_equal(ctrl$bin_width, 1)
})

test_that("control() accepts custom values", {
  ctrl <- control(years = 50, delta_time = 0.5)
  expect_equal(ctrl$years, 50)
  expect_equal(ctrl$delta_time, 0.5)
})

# --- supported_species() ---

test_that("supported_species() returns a tibble with required columns", {
  sp <- supported_species()
  expect_true(inherits(sp, "tbl_df"))
  expect_gt(nrow(sp), 0L)
  expect_true("species_id" %in% names(sp))
  expect_true("common_name" %in% names(sp))
})

# --- stand() ---

test_that("stand() returns ipm_stand with $trees, $species, $plot_size", {
  df <- data.frame(
    size_mm = c(150, 200, 300),
    species_id = "ABIBAL",
    plot_size = 1000
  )
  s <- stand(df)
  expect_s3_class(s, "ipm_stand")
  expect_true(!is.null(s$trees))
  expect_true(!is.null(s$species))
  expect_true(!is.null(s$plot_size))
})

test_that("stand() rejects trees smaller than 127 mm with an error", {
  df <- data.frame(
    size_mm = c(100, 200),
    species_id = "ABIBAL",
    plot_size = 1000
  )
  expect_error(stand(df))
})

# --- species_model() ---

test_that("species_model() returns ipm_spModel with $species, $params, $on_missing", {
  df <- data.frame(size_mm = c(150, 200, 300), species_id = "ABIBAL", plot_size = 1000)
  s <- stand(df)
  mod <- species_model(s)
  expect_s3_class(mod, "ipm_spModel")
  expect_true("species" %in% names(mod))
  expect_true("params" %in% names(mod))
  expect_true("on_missing" %in% names(mod))
})

test_that("species_model() errors on unknown species", {
  df <- data.frame(size_mm = c(150, 200), species_id = "FAKESPP", plot_size = 1000)
  s <- stand(df)
  expect_error(species_model(s))
})

# --- parameters() ---

test_that("parameters() returns ipm_parameters with required structure", {
  df <- data.frame(size_mm = c(150, 200, 300), species_id = "ABIBAL", plot_size = 1000)
  s <- stand(df)
  mod <- species_model(s)
  p <- parameters(mod, draw = "random", seed = 42L)
  expect_s3_class(p, "ipm_parameters")
  expect_true("species_params" %in% names(p))
  expect_true("draw_type" %in% names(p))
  expect_true("seed" %in% names(p))
})
