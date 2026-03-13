# Phase 5 Context: Compare Old vs New Package Results

## Summary

Create a standalone comparison script in a new `simulations/` subfolder. The script runs a stratified subset of the `covariates_perturbation` simulation with the new package code and checks that all lambda values match `final_output.RDS` within 1e-10.

---

## Decisions

### A. Source of Truth

- **File:** `simulations/covariates_perturbation/final_output.RDS`
- **Structure:** A list where element `[[i]]` corresponds to row `i` of `simulation_pars.RDS`. Each element contains 100 replications with columns: `lambda_base`, `par.BA_con`, `par.BA_het`, `par.temp`, `par.prec`.
- **Reproducibility guarantee:** `set.seed(i)` in `ipm_i.R` means row `i` is exactly reproducible — this is why the file can be used as a deterministic reference.

### B. Comparison Scope

Compare all 5 metrics per replication:
- `lambda_base` — lambda at reference conditions
- `par.BA_con` — sensitivity to conspecific basal area perturbation
- `par.BA_het` — sensitivity to heterospecific basal area perturbation
- `par.temp` — sensitivity to temperature perturbation
- `par.prec` — sensitivity to precipitation perturbation

Perturbation size is 0.01 in scaled covariate space (same as original `ipm_i.R`).

### C. Tolerance

Maximum absolute difference per lambda value: **1e-10**. Apply `all(abs(new - old) < 1e-10)` across all 100 replications × 5 metrics for each selected simulation row.

### D. Script Location and Format

- **New folder:** `simulations/compare_versions/` (not git-tracked, add to `.gitignore` or leave untracked)
- **Script type:** Standalone R script (not testthat) — runs interactively
- **No output file required** — print a pass/fail table to console; fail entries show max absolute difference

### E. Simulation Row Selection

Use a **stratified sample** from `simulation_pars.RDS` to cover:
- Multiple species (e.g., 3–5 species)
- Rows with and without heterospecific competition (`het_dbh` NULL vs non-NULL)
- A range of climate conditions

~10–20 rows total is sufficient.

### F. New API Workflow

The script uses `devtools::load_all()` to load the package, then adapts the `ipm_i.R` loop to the new API. The primary target is testing the core math changes:
- `getEigenValues()` (RcppEigen, Phase 02-01) replaces `max(Re(eigen(...)$values))`
- Vectorized `mkKernel()` (Phase 06) replaces the original loop-based version
- `dnorm/pnorm` in `ingrowth_lk` (Phase 06-03) replaces `truncnorm::dtruncnorm`

---

## Parameter Draw Challenge (Researcher/Planner Must Resolve)

The old `ipm_i.R` loads parameters as:
```r
set.seed(array_id)
pop_pars <- readRDS('pop_pars.RDS') |> filter(species_id == Sp) |> slice_sample(n = 100)
plot_pars <- readRDS(paste0('plot_parameters/', Sp, '.RDS')) |> ...
# then for rep i:
pop_pars_i <- pars_to_list(pop_pars[i, ])
re_pars_i  <- plot_pars[i, ] |> unlist() |> unname()
```

The new `parameters()` function selects a single draw by draw index (1–2000):
```r
parameters(mod, draw = draw_id)   # draw_id ∈ 1:2000
```

These are **different draw selection mechanisms**. The planner must choose one of:

1. **Direct `mkKernel` + `getEigenValues` approach** *(recommended)*: Call `devtools::load_all()`, load `pop_pars.RDS` and `plot_parameters/` directly (same as old code), run the same parameter loop, but replace `max(Re(eigen(K)$values))` with `max(getEigenValues(K))`. This directly tests the math changes while guaranteeing identical parameter inputs.

2. **Full new API pipeline**: Use `stand()` → `species_model()` → `parameters(draw = i, seed = row_id)` → `lambda()`. This is only viable if the researcher confirms that `species_model()` can read from `simulations/covariates_perturbation/pop_pars.RDS` and that the draw index maps to the same posterior row as the old `slice_sample` with the same seed.

The planner should investigate which approach ensures **identical parameters** between old and new runs. Exact match at 1e-10 requires bit-identical parameter inputs.

---

## Code Context

### Old simulation entry points
- `simulations/covariates_perturbation/ipm_i.R` — per-row simulation logic; 100 replications, 5 lambda computations each
- `simulations/covariates_perturbation/simulation_pars.RDS` — rows with species_id, plot_id, year_measured, bio_01_mean, bio_12_mean, plot_size, con_dbh (list), het_dbh (list or NULL)
- `simulations/covariates_perturbation/pop_pars.RDS` — posterior draws, all species
- `simulations/covariates_perturbation/plot_parameters/{Sp}.RDS` — per-species plot random effects
- `simulations/covariates_perturbation/final_output.RDS` — reference output (list indexed by row)

### New package functions under test
- `R/kernel.R`: `mkKernel()` — kernel assembly; vectorized after Phase 06
- `src/eigen.cpp` + `R/RcppExports.R`: `getEigenValues()` — RcppEigen solver (Phase 02-01)
- `R/vital_rates.R`: `ingrowth_lk()` — uses `dnorm/pnorm` after Phase 06-03
- `R/BasalArea_competition.R`: `dbh_to_sizeDist()`, `size_to_BAplot()` — unchanged

### Scaling utilities (must be reproduced in comparison script)
The old code uses `scale_vars()` (local helper in `ipm_i.R`) and `dbh_range.RDS` for DBH scaling. The comparison script needs the same `vars_rg` object from `../data/climate_scaleRange.RDS` and `dbh_range.RDS`.

### Perturbation logic (from `ipm_i.R`)
- BA_con perturbation: scale each con_dbh value, add 0.01, unscale → new N_con_pertb
- BA_het perturbation: same for het_dbh (skip if NULL, lambda_bahet = lambda_base)
- Temp perturbation: `temp_scl + 0.01`
- Prec perturbation: `prec_scl + 0.01`
- Sensitivity: `abs(lambda_perturbed - lambda_base) / delta_covariate`
