# profile_benchmark.R
# Standalone profiling script — not shipped with the package.
# Run from the package root after devtools::load_all().
# Requires: profvis, microbenchmark (install if missing)
#
# Usage:
#   devtools::load_all()
#   source(".planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profile_benchmark.R")

# -- Install profiling packages if missing -----------------------------------------

if (!requireNamespace("profvis",        quietly = TRUE)) install.packages("profvis")
if (!requireNamespace("microbenchmark", quietly = TRUE)) install.packages("microbenchmark")
if (!requireNamespace("htmlwidgets",    quietly = TRUE)) install.packages("htmlwidgets")

devtools::load_all(".")
library(profvis)
library(microbenchmark)

phase_dir <- ".planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains"

# -- Benchmark inputs ---------------------------------------------------------------

# Single-species: QUERUB (Lmax = 2361mm — largest mesh, stresses outer() kernel build)
df_single <- data.frame(
  size_mm    = seq(130, 600, by = 20),   # 24 individuals across size range
  species_id = "QUERUB",
  plot_size  = 400
)
s_single    <- stand(df_single)
mod_single  <- species_model(s_single)
pars_single <- parameters(mod_single, draw = "mean")
env         <- env_condition(MAT = 10, MAP = 1100)
ctrl_5      <- control(years = 5,  compute_lambda = FALSE, progress = FALSE)
ctrl_100    <- control(years = 100, compute_lambda = FALSE, progress = FALSE)
ctrl_1      <- control(years = 1,  compute_lambda = FALSE)

# Multi-species: QUERUB + QUEPRI (exposes competition update and dual-species cost)
df_multi <- data.frame(
  size_mm    = c(seq(130, 600, by = 20), seq(130, 400, by = 30)),
  species_id = c(rep("QUERUB", length(seq(130, 600, by = 20))),
                 rep("QUEPRI", length(seq(130, 400, by = 30)))),
  plot_size  = 400
)
s_multi    <- stand(df_multi)
mod_multi  <- species_model(s_multi)
pars_multi <- parameters(mod_multi, draw = "mean")
ctrl_multi <- control(years = 3, compute_lambda = FALSE, progress = FALSE)

# Log mesh sizes
nvec_list <- forestIPM:::.stand_to_nvec(s_single, "QUERUB", pars_single, ctrl_1$bin_width)
mesh_n    <- length(nvec_list[["QUERUB"]]$N_con$meshpts)
message(sprintf("QUERUB mesh: %d points (Lmax=2361mm, bin_width=1mm, range 127-2362mm)", mesh_n))
message(sprintf("Kernel dimensions: %d x %d = %d cells per kernel build", mesh_n, mesh_n, mesh_n^2))

# -- profvis: project() single-species (5 timesteps) --------------------------------
cat("\nRunning profvis on 5-year single-species project()...\n")
p_project_single <- profvis({
  project(mod_single, pars_single, s_single, env, ctrl_5)
}, interval = 0.005, height = "600px")

# -- profvis: project() multi-species (3 timesteps) ---------------------------------
cat("Running profvis on 3-year multi-species project()...\n")
p_project_multi <- profvis({
  project(mod_multi, pars_multi, s_multi, env, ctrl_multi)
}, interval = 0.005, height = "600px")

# -- Save flamegraphs ---------------------------------------------------------------
htmlwidgets::saveWidget(
  p_project_single,
  file.path(phase_dir, "profvis_output.html"),
  selfcontained = TRUE
)
cat("Saved:", file.path(phase_dir, "profvis_output.html"), "\n")

htmlwidgets::saveWidget(
  p_project_multi,
  file.path(phase_dir, "profvis_multi_output.html"),
  selfcontained = TRUE
)
cat("Saved:", file.path(phase_dir, "profvis_multi_output.html"), "\n")

# -- microbenchmark: component isolation --------------------------------------------
cat("\nRunning microbenchmark on mkKernel() components (10 reps each)...\n")

sp_pars      <- pars_single$species_params[["QUERUB"]]$fixed
meshpts      <- nvec_list[["QUERUB"]]$N_con$meshpts
h            <- nvec_list[["QUERUB"]]$N_con$h

# Pre-compute BA metrics once (isolate kernel outer calls)
BA_comp_intra <- size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
BA_comp_inter <- size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con,
                                N_inter = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
BAplot_intra  <- size_to_BAplot(N = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
BAplot_inter  <- size_to_BAplot(N = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
BAplot_total  <- BAplot_intra + BAplot_inter

mb_components <- microbenchmark(
  mkKernel_full = {
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
  P_outer_only = {
    h * outer(meshpts, meshpts, P_xEC,
      ctrl_1$delta_time, BA_comp_intra, BA_comp_inter, env$.MAT_scl, env$.MAP_scl,
      sp_pars[["growth"]], sp_pars[["mort"]], c(0, 0, 0))
  },
  F_outer_only = {
    h * outer(meshpts, meshpts, ingrowth_lk,
      ctrl_1$delta_time, 400, BAplot_intra, BAplot_total, env$.MAT_scl, env$.MAP_scl,
      sp_pars[["rec"]], sp_pars[["sizeIngrowth"]], 0)
  },
  BA_metrics_only = {
    size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
    size_to_BAcomp(N_intra = nvec_list[["QUERUB"]]$N_con,
                   N_inter = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
    size_to_BAplot(N = nvec_list[["QUERUB"]]$N_con, plot_size = 400)
    size_to_BAplot(N = nvec_list[["QUERUB"]]$N_het, plot_size = 400)
  },
  times = 10L
)
print(mb_components)

# -- Print profvis hotspot summary ---------------------------------------------------
cat("\n=== profvis hotspot summary (single-species, 5 yr) ===\n")
prof_df        <- p_project_single$x$message$prof
total_samples  <- max(prof_df$time)
interval_ms    <- p_project_single$x$message$interval
labels_df      <- aggregate(rep(1, nrow(prof_df)), list(label = prof_df$label), sum)
labels_df$pct  <- round(100 * labels_df$x / total_samples, 1)
labels_df      <- labels_df[order(-labels_df$x), ]
cat(sprintf("Total elapsed: ~%d ms (%d samples x %d ms interval)\n",
    total_samples * interval_ms, total_samples, interval_ms))
print(head(labels_df, 20))
