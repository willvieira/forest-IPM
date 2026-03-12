# benchmark_before_after.R
# Before/after benchmark for Plan 06-03 optimizations.
# Measures mkKernel() performance after Target 1 (vectorize outer) + Target 2 (truncnorm->dnorm/pnorm).
# Run from the package root after devtools::load_all().
#
# Usage:
#   devtools::load_all()
#   source(".planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/benchmark_before_after.R")

if (!requireNamespace("microbenchmark", quietly = TRUE)) install.packages("microbenchmark")
library(microbenchmark)

devtools::load_all(".")

# -- Setup inputs (same as profile_benchmark.R) -----------------------------------
df_single <- data.frame(
  size_mm    = seq(130, 600, by = 20),
  species_id = "QUERUB",
  plot_size  = 400
)
s_single    <- stand(df_single)
mod_single  <- species_model(s_single)
pars_single <- parameters(mod_single, draw = "mean")
env         <- env_condition(MAT = 10, MAP = 1100)
ctrl_1      <- control(years = 1, compute_lambda = FALSE)

nvec_list <- forestIPM:::.stand_to_nvec(s_single, "QUERUB", pars_single, ctrl_1$bin_width)
sp_pars   <- pars_single$species_params[["QUERUB"]]$fixed
meshpts   <- nvec_list[["QUERUB"]]$N_con$meshpts
h         <- nvec_list[["QUERUB"]]$N_con$h
n         <- length(meshpts)

cat(sprintf("QUERUB mesh: %d points — kernel dimensions %d x %d = %d cells\n",
            n, n, n, n^2))

# Pre-compute BA metrics (isolate kernel build cost from BA cost)
BA_comp_intra <- size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
BA_comp_inter <- size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con,
                                N_inter  = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
BAplot_intra  <- size_to_BAplot(N = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
BAplot_inter  <- size_to_BAplot(N = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
BAplot_total  <- BAplot_intra + BAplot_inter

# -- Full mkKernel() after optimization (50 reps) --------------------------------
cat("\n=== Full mkKernel() — AFTER optimization (50 reps) ===\n")
mb_optimized <- microbenchmark(
  mkKernel_optimized = {
    mkKernel(
      Nvec_intra  = nvec_list[["QUERUB"]]$N_con,
      Nvec_inter  = nvec_list[["QUERUB"]]$N_het,
      delta_time  = ctrl_1$delta_time,
      plotSize    = s_single$plot_size,
      Temp        = env$.MAT_scl,
      Prec        = env$.MAP_scl,
      pars        = sp_pars,
      plot_random = c(0, 0, 0)
    )
  },
  times = 50L
)
print(mb_optimized)

# -- Component-level: vectorized P and F (after) ---------------------------------
cat("\n=== Component benchmarks — AFTER optimization (50 reps each) ===\n")

s1_vec <- rep(meshpts, times = n)
s0_vec <- rep(meshpts, each  = n)

mb_components_after <- microbenchmark(
  P_vectorized = {
    P_vals <- P_xEC(
      s1_vec, s0_vec,
      ctrl_1$delta_time, BA_comp_intra, BA_comp_inter, env$.MAT_scl, env$.MAP_scl,
      sp_pars[["growth"]], sp_pars[["mort"]], c(0, 0, 0)
    )
    h * matrix(P_vals, nrow = n, ncol = n)
  },
  F_vectorized_no_truncnorm = {
    mu  <- sp_pars[["sizeIngrowth"]]["size_int"] + sp_pars[["sizeIngrowth"]]["phi_time"] * ctrl_1$delta_time
    sig <- sp_pars[["sizeIngrowth"]]["sigma_size"]
    ingrowth_scalar <- ingrowth_f(
      sp_pars[["rec"]], ctrl_1$delta_time, 400, BAplot_intra, BAplot_total,
      env$.MAT_scl, env$.MAP_scl, 0
    )
    F_vals <- ingrowth_scalar * dnorm(s1_vec, mean = mu, sd = sig) / (1 - pnorm(127, mean = mu, sd = sig))
    h * matrix(F_vals, nrow = n, ncol = n)
  },
  times = 50L
)
print(mb_components_after)

# -- Baseline comparison: original outer() approach (10 reps — slow) -------------
cat("\n=== Original outer() approach — BEFORE optimization (10 reps) ===\n")
cat("(Using original outer() calls — this is the slow path)\n")

# Temporarily define original ingrowth_lk using truncnorm for comparison
ingrowth_lk_original <- function(
  size_t1, size_t0, delta_time, plot_size, BA_adult_sp, BA_adult, Temp, Prec,
  parsIngrowth, parsSizeIngrowth, plot_random
) {
  ingrowth_f(parsIngrowth, delta_time, plot_size, BA_adult_sp, BA_adult, Temp, Prec, plot_random) *
  truncnorm::dtruncnorm(
    size_t1, a = 127, b = Inf,
    mean = parsSizeIngrowth["size_int"] + parsSizeIngrowth["phi_time"] * delta_time,
    sd   = parsSizeIngrowth["sigma_size"]
  )
}

mb_components_before <- microbenchmark(
  P_outer_original = {
    h * outer(meshpts, meshpts, P_xEC,
      ctrl_1$delta_time, BA_comp_intra, BA_comp_inter, env$.MAT_scl, env$.MAP_scl,
      sp_pars[["growth"]], sp_pars[["mort"]], c(0, 0, 0))
  },
  F_outer_original = {
    h * outer(meshpts, meshpts, ingrowth_lk_original,
      ctrl_1$delta_time, 400, BAplot_intra, BAplot_total, env$.MAT_scl, env$.MAP_scl,
      sp_pars[["rec"]], sp_pars[["sizeIngrowth"]], 0)
  },
  times = 10L
)
print(mb_components_before)

# -- Summary statistics -----------------------------------------------------------
cat("\n=== Summary ===\n")

extract_median_ms <- function(mb) {
  median(mb$time) / 1e6  # nanoseconds to milliseconds
}

p_before_ms <- extract_median_ms(mb_components_before[mb_components_before$expr == "P_outer_original", ])
f_before_ms <- extract_median_ms(mb_components_before[mb_components_before$expr == "F_outer_original", ])
p_after_ms  <- extract_median_ms(mb_components_after[mb_components_after$expr  == "P_vectorized", ])
f_after_ms  <- extract_median_ms(mb_components_after[mb_components_after$expr  == "F_vectorized_no_truncnorm", ])
mkkernel_after_ms <- extract_median_ms(mb_optimized)

cat(sprintf("P matrix: before=%.1f ms, after=%.1f ms, speedup=%.1fx\n",
            p_before_ms, p_after_ms, p_before_ms / p_after_ms))
cat(sprintf("F matrix: before=%.1f ms, after=%.1f ms, speedup=%.1fx\n",
            f_before_ms, f_after_ms, f_before_ms / f_after_ms))
cat(sprintf("Full mkKernel(): baseline=~761ms (Phase 06-01 measurement), after=%.1f ms, speedup=%.1fx\n",
            mkkernel_after_ms, 761.0 / mkkernel_after_ms))
cat(sprintf("Estimated 100-yr run: baseline=~76s, after=~%.0fs\n",
            mkkernel_after_ms * 100 / 1000))
