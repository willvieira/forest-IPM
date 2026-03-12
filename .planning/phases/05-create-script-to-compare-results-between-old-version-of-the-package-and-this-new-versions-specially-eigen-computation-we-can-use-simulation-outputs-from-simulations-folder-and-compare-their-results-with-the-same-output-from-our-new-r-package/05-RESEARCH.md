# Phase 5: Compare Old vs New Package Results - Research

**Researched:** 2026-03-12
**Domain:** R package regression comparison — eigenvalue solver, vectorized kernel, ingrowth truncation replacement
**Confidence:** HIGH

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **Source of truth:** `simulations/covariates_perturbation/final_output.RDS` — list where `[[i]]` is the 100-rep tibble for simulation row `i`
- **Comparison scope:** 5 metrics per rep: `lambda_base`, `par.BA_con`, `par.BA_het`, `par.temp`, `par.prec`
- **Tolerance:** `all(abs(new - old) < 1e-10)` across all 100 reps × 5 metrics per selected row
- **Script location:** `simulations/compare_versions/` — new subfolder, not git-tracked
- **Script format:** Standalone R script (not testthat), runs interactively
- **Output:** Print pass/fail table to console; fail entries show max absolute difference
- **Row selection:** ~10–20 stratified rows covering multiple species, NULL vs non-NULL `het_dbh`, range of climate
- **Perturbation size:** `0.01` in scaled covariate space (same as original `ipm_i.R`)
- **New API entry:** `devtools::load_all()` — Direct `mkKernel` + `getEigenValues` approach (see Parameter Draw Resolution below)

### Claude's Discretion

- Specific stratified row indices to select (must cover multiple species, het/no-het, climate range)
- Whether to use `tryCatch` per-row or fail-fast on first mismatch

### Deferred Ideas (OUT OF SCOPE)

- Full new API pipeline (`stand()` → `species_model()` → `parameters()` → `lambda()`)
- Extending coverage beyond the 5 listed metrics
- Automated CI integration
</user_constraints>

---

## Summary

Phase 5 produces a single standalone R script (`simulations/compare_versions/compare_versions.R`) that reproduces the old `ipm_i.R` logic using the new package internals for a stratified subset of simulation rows, then checks all 5 lambda metrics against `final_output.RDS` within tolerance 1e-10.

The core verification target is narrow: `getEigenValues()` (RcppEigen, Phase 02-01) must reproduce `max(Re(eigen(K)$values))`, and the vectorized `mkKernel()` plus the `dnorm/pnorm` ingrowth replacement must produce numerically identical kernels. The `final_output.RDS` was generated with `truncnorm::dtruncnorm`; the math equivalence has been confirmed at 6.9e-18 maximum difference (machine epsilon, well inside 1e-10).

**Primary recommendation:** Use the Direct `mkKernel` + `getEigenValues` approach — load `pop_pars.RDS` and `plot_parameters/` directly (identical to old code), reproduce `set.seed(row_id)` + `slice_sample`, and call the new functions in place of the old ones. This guarantees bit-identical parameter inputs and isolates the math changes as the only variable.

---

## Parameter Draw Resolution (Critical Finding)

The CONTEXT.md flagged this as the key decision. Research confirms:

**Resolved: Use the Direct approach (Option 1 from CONTEXT.md).**

Reasons:
1. `pop_pars.RDS` has 4000 rows per species (verified). `slice_sample(n=100)` with `set.seed(array_id)` is deterministic and reproducible on the same machine.
2. `plot_parameters/{Sp}.RDS` stores a tibble with a `plot_pars` list-column; each element is a 100-row tibble with columns `growth`, `mort`, `recruit` — directly aligned with the 100 reps.
3. The new `parameters()` API uses draw indices 1–2000 from `mod$params`, a completely different data structure from `pop_pars.RDS`. There is no guaranteed row-to-row mapping between the two.
4. `pars_to_list()` is NOT in the current `R/` files — it existed in the old `params.R` (now in `.history/`). The comparison script must inline it or copy it locally.
5. `init_pop()` is NOT in the current `kernel.R` — it was removed during package cleanup. The comparison script needs it to construct `N_ref` for `dbh_to_sizeDist`. It must be inlined from `.history/R/kernel_20260224225955.R`.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| devtools | >=2.4 | `load_all()` to source current package | Standard R package dev workflow |
| tidyverse | >=2.0 | `filter`, `select`, `slice_sample` (already a dep) | Already in DESCRIPTION Imports |
| base R | - | `set.seed`, `readRDS`, `all`, `abs`, `max`, `rbind` | No new dependencies needed |

### No New Dependencies Required
The script uses only what the package already declares (tidyverse/purrr) plus base R. It does NOT need `truncnorm` (replaced) and does NOT need any new packages.

**To load the package:**
```r
devtools::load_all()
```

---

## Architecture Patterns

### Recommended Script Structure

```
simulations/compare_versions/
└── compare_versions.R    # standalone, sources devtools::load_all() from project root
```

The script must be run with working directory set to the project root (so `devtools::load_all()` works and relative paths to `simulations/covariates_perturbation/` are valid).

### Pattern: Direct mkKernel + getEigenValues Approach

**What:** Replicate `ipm_i.R` exactly, replacing only:
1. `max(Re(eigen(K)$values))` → `max(getEigenValues(K))`
2. Old `mkKernel` (outer-based, truncnorm) → new `mkKernel` (vectorized, dnorm/pnorm)

**When to use:** Always — this is the locked decision.

**Skeleton:**
```r
# Run from project root
devtools::load_all()
library(tidyverse)

SIM_DIR  <- "simulations/covariates_perturbation"
DATA_DIR <- readLines("_data.path")

# --- helpers carried from old code ---
pars_to_list <- function(pars) { ... }   # inlined from .history/R/params_20260224225907.R
init_pop     <- function(params, L, h, N, ...) { ... }  # inlined from .history/R/kernel_20260224225955.R
scale_vars   <- function(value, cov, direction, range_dt) { ... }  # inline from ipm_i.R

# --- load scaling ranges ---
vars_rg <- append(
  setNames(readRDS(file.path(DATA_DIR, "climate_scaleRange.RDS")), c("temp", "prec")),
  list("dbh" = readRDS(file.path(SIM_DIR, "dbh_range.RDS")))
)

# --- stratified rows (selected subset ~15 rows) ---
selected_rows <- c(...)

# --- main comparison loop ---
results <- lapply(selected_rows, function(array_id) {
  sim_pars <- readRDS(file.path(SIM_DIR, "simulation_pars.RDS"))[array_id, ]
  Sp       <- sim_pars$species_id
  set.seed(array_id)
  pop_pars <- readRDS(file.path(SIM_DIR, "pop_pars.RDS")) |>
    filter(species_id == Sp) |>
    select(-species_id) |>
    slice_sample(n = 100)
  plot_pars <- readRDS(file.path(SIM_DIR, "plot_parameters", paste0(Sp, ".RDS"))) |>
    filter(plot_id == sim_pars$plot_id, year_measured == sim_pars$year_measured) |>
    pull(plot_pars) |> .[[1]]
  ref_output <- readRDS(file.path(SIM_DIR, "final_output.RDS"))[[array_id]]

  # replicate 100-rep loop...
  # compare new vs ref_output[, c("lambda_base","par.BA_con","par.BA_het","par.temp","par.prec")]
})
```

### Pattern: pars_to_list (must be inlined)

The function converts a wide tibble row of `growth.*`, `mort.*`, `rec.*`, `sizeIngrowth.*` columns into a named list of named numeric vectors — the format `mkKernel(pars = ...)` expects.

```r
pars_to_list <- function(pars) {
  pars |>
    select(contains(".")) |>
    pivot_longer(cols = everything()) |>
    mutate(
      vr  = str_replace(name, "\\..*", ""),
      par = str_replace(name, paste0(vr, "."), ""),
      vr  = case_match(vr, "recruit" ~ "rec", .default = vr)
    ) |>
    select(-name) |>
    group_split(vr) |>
    set_names(map_chr(., ~.x$vr[1])) |>
    map(~.x |> select(-vr) |> pivot_wider(names_from = par) |> unlist(use.names = TRUE))
}
```

### Pattern: init_pop (must be inlined, simplified)

The old `init_pop()` with `N=0` only creates a mesh grid (no stochastic sampling needed):

```r
init_pop <- function(params, L = 127, h = 1, N = 0) {
  Lmax <- round(params[["growth"]]["Lmax"], 0)
  m    <- length(seq(L, Lmax, h))
  msh  <- L + ((seq_len(m)) - 0.5) * h
  if (N == 0) {
    return(list(meshpts = msh, Nvec = rep(0, m), h = h))
  }
  stop("N > 0 not needed for comparison script")
}
```

The comparison script only ever calls `init_pop(..., N = 0)` (to build the mesh grid that `dbh_to_sizeDist` bins into). The full stochastic path (`N > 0`) is never exercised.

### Anti-Patterns to Avoid

- **Calling `parameters()` + `species_model()`:** The new API uses a different parameter structure from `pop_pars.RDS` — no guaranteed row-to-row mapping with `set.seed` + `slice_sample`. Would not produce identical inputs.
- **Sourcing old `ipm_i.R` directly:** It uses `source("R/vital_rates.R")` etc., which references the old functions (with `truncnorm`, old `outer`-based kernel). Would undermine the test.
- **Using `pars_to_list` from `purrr::as_vector()`:** The old history version uses `as_vector()` which is deprecated in purrr >= 1.0.0. Use `unlist(use.names = TRUE)` instead (already fixed in BUG-03).

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Eigenvalue computation | Custom R loop | `getEigenValues(K)` | That's the whole point of the test |
| Kernel assembly | Re-implementing P/F matrix | `mkKernel()` after `load_all()` | The new vectorized version is what's being verified |
| Size distribution binning | Custom cut/table logic | `dbh_to_sizeDist()` | Already in package, unchanged |
| BA metric | Custom BA formula | `size_to_BAplot()`, `size_to_BAcomp()` | Already in package, unchanged |

**Key insight:** The comparison script is thin glue — it sets up inputs from the old data files, calls the new package functions, and diffs outputs. Almost nothing should be hand-rolled.

---

## Common Pitfalls

### Pitfall 1: Wrong working directory when sourcing helpers
**What goes wrong:** `devtools::load_all()` fails or data paths are wrong.
**Why it happens:** Script is run from `simulations/compare_versions/` instead of project root.
**How to avoid:** Document at top of script: "Run from project root with `source('simulations/compare_versions/compare_versions.R')`". Add `stopifnot(file.exists("DESCRIPTION"))` as guard.
**Warning signs:** "Cannot open file 'DESCRIPTION'" error on `load_all()`.

### Pitfall 2: set.seed must precede ALL random calls
**What goes wrong:** Sampling in `pop_pars` produces different rows than `final_output.RDS`.
**Why it happens:** Any random call before `set.seed(array_id)` shifts the RNG state.
**How to avoid:** Call `set.seed(array_id)` immediately before `filter(species_id == Sp) |> slice_sample(n = 100)`, with no intervening random calls.
**Warning signs:** Non-zero max diff for the very first rep.

### Pitfall 3: final_output index offset
**What goes wrong:** `final_output.RDS[[i]]` returns results for the wrong row.
**Why it happens:** Assuming 0-based indexing, or that element `i` = array_id `i`.
**How to avoid:** Verified — `final_output[[i]]$array_id[1] == i` for all tested elements. The file is 1-indexed and index-aligned. 6 rows are missing (IDs: 32794, 33256, 33257, 171478, 195665, 201283) — skip these in stratified selection.
**Warning signs:** Selecting a missing row ID causes `is.null(final_output[[i]])`.

### Pitfall 4: pars_to_list uses deprecated purrr::as_vector()
**What goes wrong:** `Error: as_vector() was deprecated in purrr 1.0.0`.
**Why it happens:** The `.history/` version uses `as_vector()`.
**How to avoid:** Replace with `unlist(use.names = TRUE)` in the inlined copy (BUG-03 already established this fix).

### Pitfall 5: het_dbh NULL handling must match original logic exactly
**What goes wrong:** Perturbation logic diverges for NULL `het_dbh` rows.
**Why it happens:** Missing the `if(is.null(unlist(sim_pars$het_dbh)))` guard around BA_het perturbation.
**How to avoid:** Copy the exact NULL guard from `ipm_i.R` lines 185–217. When `het_dbh` is NULL, `lambda_bahet = lambda_base` (no perturbation computed).

### Pitfall 6: Sensitivity formula uses absolute difference (abs)
**What goes wrong:** Sign error in `par.BA_con` etc.
**Why it happens:** The formula in `ipm_i.R` is `abs(lambda_perturbed - lambda_base) / delta_covariate`.
**How to avoid:** Use `abs()` in the numerator, matching `ipm_i.R` lines 283–292 exactly.

### Pitfall 7: climate_scaleRange.RDS is outside the sim directory
**What goes wrong:** File not found when constructing `vars_rg`.
**Why it happens:** `ipm_i.R` accesses it as `'../data/climate_scaleRange.RDS'` — relative to `simulations/covariates_perturbation/`, meaning it's at `simulations/data/`.
**How to avoid:** Construct the path from `_data.path`: `file.path(readLines("_data.path"), "climate_scaleRange.RDS")`. Verified: file exists at `/Users/wvieira/ownCloud/thesisData/climate_scaleRange.RDS` and returns `list(bio_01_mean = c(-5.60, 25.14), bio_12_mean = c(65, 2688))`.

---

## Code Examples

### Reproducing the exact parameter loading sequence
```r
# Source: ipm_i.R lines 24-48, adapted for non-cluster use
array_id <- selected_rows[1]
set.seed(array_id)   # MUST be before any random call

sim_pars <- readRDS(file.path(SIM_DIR, "simulation_pars.RDS"))[array_id, ]
Sp       <- sim_pars$species_id

pop_pars <- readRDS(file.path(SIM_DIR, "pop_pars.RDS")) |>
  filter(species_id == Sp) |>
  select(-species_id) |>
  slice_sample(n = 100)   # 100 draws per rep, same seed as original

plot_pars <- readRDS(file.path(SIM_DIR, "plot_parameters", paste0(Sp, ".RDS"))) |>
  filter(plot_id == sim_pars$plot_id, year_measured == sim_pars$year_measured) |>
  pull(plot_pars) |> .[[1]]   # 100-row tibble: growth, mort, recruit
```

### The new eigenvalue call replacing the old one
```r
# OLD (ipm_i.R):
lambda_base <- max(Re(eigen(mkKernel(...)$K)$values))

# NEW (comparison script):
lambda_base <- max(getEigenValues(mkKernel(...)$K))
```

### Stratified row selection
```r
# Source: analysis of simulation_pars.RDS structure
# 31 species, 216027 rows, 4947 NULL het_dbh, 211080 non-NULL het_dbh
# Missing IDs (must exclude): 32794, 33256, 33257, 171478, 195665, 201283
sim_pars <- readRDS(file.path(SIM_DIR, "simulation_pars.RDS"))
missing   <- c(32794L, 33256L, 33257L, 171478L, 195665L, 201283L)

# Example stratified selection (planner finalizes exact indices):
# - 4 species with het_dbh non-NULL, spanning temp/prec range
# - 2 rows where het_dbh IS NULL (species with NULL present in rows with IDs != missing)
# - Cover cold (temp < 5), temperate (5-15), warm (temp > 15) climate zones
```

### Tolerance check per row
```r
# Source: CONTEXT.md Decision C
check_row <- function(new_out, ref_out, metrics = c("lambda_base","par.BA_con","par.BA_het","par.temp","par.prec")) {
  diffs <- sapply(metrics, function(m) max(abs(new_out[[m]] - ref_out[[m]])))
  pass  <- all(diffs < 1e-10)
  list(pass = pass, max_diffs = diffs)
}
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `max(Re(eigen(K)$values))` | `max(getEigenValues(K))` | Phase 02-01 | RcppEigen EigenSolver, correct non-symmetric solver |
| `outer(meshpts, meshpts, FUN)` in mkKernel | `rep()/matrix()` vectorized | Phase 06-01 | 1.4x speedup, mathematically identical |
| `truncnorm::dtruncnorm(x, a=127, b=Inf, mu, sig)` | `dnorm(x,mu,sig)/(1-pnorm(127,mu,sig))` | Phase 06-03 | 2.4x F-matrix speedup; max diff from final_output.RDS: 6.9e-18 |
| `exists('fct')` global check in init_pop | `fct <- 1` initialized before loop | Phase 02-01 | BUG-02 fix; comparison script inlines simplified N=0-only version |
| `purrr::as_vector()` | `unlist(use.names=TRUE)` | Phase 02-01 | BUG-03 fix; must use new form in inlined pars_to_list |

**Key verified finding:** `dnorm/pnorm` formula vs `truncnorm::dtruncnorm` produces maximum absolute difference of 6.9e-18 (< 1e-10 tolerance). The 1e-10 threshold is achievable.

---

## Data Files: What Exists and Where

| File | Path | Structure | Size |
|------|------|-----------|------|
| `simulation_pars.RDS` | `simulations/covariates_perturbation/` | tibble, 216027 rows × 11 cols | ~15MB |
| `pop_pars.RDS` | `simulations/covariates_perturbation/` | tibble, 124000 rows (4000/species × 31 species) | ~20MB |
| `plot_parameters/{Sp}.RDS` | `simulations/covariates_perturbation/plot_parameters/` | grouped tibble: plot_id, year_measured, plot_pars (list-col of 100-row tibbles) | per species |
| `final_output.RDS` | `simulations/covariates_perturbation/` | list of 216021 elements; `[[i]]` is 100×15 tibble | 842MB |
| `dbh_range.RDS` | `simulations/covariates_perturbation/` | numeric(2): c(127, 1351.28) | tiny |
| `climate_scaleRange.RDS` | `{_data.path}/` | list: bio_01_mean=c(-5.60, 25.14), bio_12_mean=c(65, 2688) | tiny |

**Missing final_output rows (skip in stratified selection):** 32794, 33256, 33257, 171478, 195665, 201283

---

## Open Questions

1. **Which exact 15 row indices to hard-code in the script?**
   - What we know: 31 species, 4947 NULL het_dbh rows, 211080 non-NULL; temp range −4.64 to 24.74 °C; prec range 320–2232 mm
   - What's unclear: Planner needs to run a one-time selection query against `simulation_pars.RDS` to pick ~5 species × 3 climate zones × het/no-het = ~15 rows
   - Recommendation: Planner task includes a code snippet that selects and hard-codes the indices; they should be constant in the committed script so runs are reproducible.

2. **Should the script reload `simulation_pars.RDS` and `pop_pars.RDS` once (outside the loop) or per-row?**
   - What we know: `final_output.RDS` is 842MB and must NOT be loaded per-row; `pop_pars.RDS` is 20MB.
   - Recommendation: Load `simulation_pars.RDS`, `pop_pars.RDS`, and `final_output.RDS` once at the top, outside the row loop. Filter inside the loop.

---

## Validation Architecture

> `workflow.nyquist_validation` key is absent from `.planning/config.json` — treating as enabled.

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Standalone R script (not testthat) — per CONTEXT.md Decision D |
| Config file | none — interactive script |
| Quick run command | `source("simulations/compare_versions/compare_versions.R")` from project root |
| Full suite command | same — script is its own suite |

### Phase Requirements → Test Map
| Behavior | Test Type | Automated Command | File Exists? |
|----------|-----------|-------------------|-------------|
| `getEigenValues()` matches `max(Re(eigen()$values))` within 1e-10 | regression | `source("simulations/compare_versions/compare_versions.R")` | Wave 0 |
| Vectorized `mkKernel()` matches old `outer()`-based kernel | regression (implicit, via lambda) | same | Wave 0 |
| `dnorm/pnorm` ingrowth matches `truncnorm::dtruncnorm` | regression (implicit, via lambda) | same | Wave 0 |
| All 5 sensitivity derivatives match within 1e-10 | regression | same | Wave 0 |
| NULL `het_dbh` rows handled correctly (`lambda_bahet = lambda_base`) | regression | same | Wave 0 |

### Sampling Rate
- **Per task commit:** Run script manually against 1 selected row to confirm it executes
- **Per wave merge:** Run full script against all ~15 selected rows
- **Phase gate:** All rows PASS printed to console before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `simulations/compare_versions/compare_versions.R` — the entire deliverable of this phase
- [ ] `simulations/compare_versions/` directory — must be created

*(No existing test infrastructure covers this phase — the comparison script IS the test.)*

---

## Sources

### Primary (HIGH confidence)
- Direct inspection of `simulations/covariates_perturbation/ipm_i.R` — full logic verified
- Direct inspection of `R/kernel.R`, `R/vital_rates.R`, `R/BasalArea_competition.R`, `R/lambda.R`, `R/stand.R`, `R/parameters.R`, `src/eigen.cpp` — all current function signatures and logic confirmed
- `.history/R/kernel_20260224225955.R` — `init_pop()` and old `mkKernel()` (with truncnorm + outer) verified
- `.history/R/params_20260224225907.R` — `pars_to_list()` implementation verified
- R runtime verification: pop_pars structure (4000 draws/species × 31 species = 124000 rows), plot_parameters structure (100-row plot_pars tibble), final_output alignment (index = array_id, 216021 elements, 6 missing)
- R runtime verification: `truncnorm::dtruncnorm` vs `dnorm/pnorm` max diff = 6.9e-18
- R runtime verification: `climate_scaleRange.RDS` at `_data.path` — ranges confirmed

### Secondary (MEDIUM confidence)
- STATE.md: "regression tests confirmed 1e-10 tolerance" for Phase 06 optimizations — corroborates that the math changes are within the tolerance target

---

## Metadata

**Confidence breakdown:**
- Parameter draw approach: HIGH — confirmed by runtime inspection of both data structures; mapping impossibility between old and new APIs verified
- File structure: HIGH — all files read and structure confirmed with R
- Math equivalence (truncnorm vs dnorm/pnorm): HIGH — runtime verified at 6.9e-18 max diff
- init_pop/pars_to_list inlining need: HIGH — confirmed absent from current R/ files
- Tolerance achievability: HIGH — math equivalence verified, regression tests in Phase 06 already used 1e-10

**Research date:** 2026-03-12
**Valid until:** 2026-04-12 (stable codebase; only invalidated by changes to kernel.R, vital_rates.R, or eigen.cpp)
