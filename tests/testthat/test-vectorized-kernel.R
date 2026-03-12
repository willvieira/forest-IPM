# tests/testthat/test-vectorized-kernel.R
# TDD RED phase: verify that vectorized rep()/matrix() pattern produces
# element-wise identical output to original outer() for a tiny 3x3 kernel.
# This test will FAIL until the vectorized implementation is in kernel.R.

test_that("vectorized P matrix equals outer() result for 3x3 micro-kernel", {
  # Requires package loaded
  s <- stand(data.frame(
    size_mm    = seq(130, 600, by = 50),
    species_id = "QUERUB",
    plot_size  = 1000
  ))
  pars <- parameters(species_model(s), draw = "mean")
  sp_pars <- pars$species_params[["QUERUB"]]$fixed

  test_meshpts <- c(130, 200, 300)
  n_test       <- length(test_meshpts)

  # Build tiny Nvec structure using 3 mesh points
  # Reuse size_to_BAcomp with a minimal Nvec structure
  h_test <- test_meshpts[2] - test_meshpts[1]  # 70 mm

  # Minimal Nvec list for BA computation
  Nvec_mini <- list(
    meshpts = test_meshpts,
    Nvec    = rep(10, n_test),
    h       = h_test
  )
  BA_intra <- size_to_BAcomp(N_intra = Nvec_mini, plot_size = 1000)
  BA_inter <- size_to_BAcomp(N_intra = Nvec_mini, N_inter = Nvec_mini, plot_size = 1000)

  # Original outer() result
  P_outer <- h_test * outer(
    test_meshpts, test_meshpts,
    P_xEC,
    1, BA_intra, BA_inter, 10, 1100,
    sp_pars[["growth"]], sp_pars[["mort"]], c(0, 0, 0)
  )

  # Vectorized result (rep + direct eval + matrix)
  # outer(X,Y,FUN): X_expanded=rep(X,times=n), Y_expanded=rep(Y,each=n)
  s1 <- rep(test_meshpts, times = n_test)   # outer's X_expanded -> 1st arg (size_t1)
  s0 <- rep(test_meshpts, each  = n_test)   # outer's Y_expanded -> 2nd arg (size_t0)
  P_vals <- P_xEC(
    s1, s0,
    1, BA_intra, BA_inter, 10, 1100,
    sp_pars[["growth"]], sp_pars[["mort"]], c(0, 0, 0)
  )
  P_vec <- h_test * matrix(P_vals, nrow = n_test, ncol = n_test)

  expect_equal(P_outer, P_vec, tolerance = .Machine$double.eps)
})

test_that("vectorized F matrix equals outer() result for 3x3 micro-kernel", {
  s <- stand(data.frame(
    size_mm    = seq(130, 600, by = 50),
    species_id = "QUERUB",
    plot_size  = 1000
  ))
  pars <- parameters(species_model(s), draw = "mean")
  sp_pars <- pars$species_params[["QUERUB"]]$fixed

  test_meshpts <- c(130, 200, 300)
  n_test       <- length(test_meshpts)
  h_test       <- test_meshpts[2] - test_meshpts[1]

  Nvec_mini <- list(
    meshpts = test_meshpts,
    Nvec    = rep(10, n_test),
    h       = h_test
  )
  BAplot_intra <- size_to_BAplot(N = Nvec_mini, plot_size = 1000)
  BAplot_inter <- size_to_BAplot(N = Nvec_mini, plot_size = 1000)
  BAplot_total <- BAplot_intra + BAplot_inter

  # Original outer() result
  F_outer <- h_test * outer(
    test_meshpts, test_meshpts,
    ingrowth_lk,
    1, 1000, BAplot_intra, BAplot_total, 10, 1100,
    sp_pars[["rec"]], sp_pars[["sizeIngrowth"]], 0
  )

  # Vectorized result
  # outer(X,Y,FUN): X_expanded=rep(X,times=n), Y_expanded=rep(Y,each=n)
  s1 <- rep(test_meshpts, times = n_test)   # outer's X_expanded -> 1st arg (size_t1)
  s0 <- rep(test_meshpts, each  = n_test)   # outer's Y_expanded -> 2nd arg (size_t0)
  F_vals <- ingrowth_lk(
    s1, s0,
    1, 1000, BAplot_intra, BAplot_total, 10, 1100,
    sp_pars[["rec"]], sp_pars[["sizeIngrowth"]], 0
  )
  F_vec <- h_test * matrix(F_vals, nrow = n_test, ncol = n_test)

  expect_equal(F_outer, F_vec, tolerance = .Machine$double.eps)
})
