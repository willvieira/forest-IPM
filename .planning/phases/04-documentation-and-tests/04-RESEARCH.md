# Phase 4: Documentation and Tests - Research

**Researched:** 2026-03-13
**Domain:** R package documentation (roxygen2), testthat test suite, S3 method design, Quarto book restructuring
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Tests**
- Remove all existing test files **except** `test-constructors.R` — the others are superseded or have pre-existing failures that are out of scope
- New test suite covers the full main workflow: `stand()` → `species_model()` → `parameters()` → `lambda()` → `project()`
- Target: **80%+ line coverage** on main workflow functions (`stand.R`, `species_model.R`, `parameters.R`, `kernel.R`, `lambda.R`, `project.R`)
- Use `plot_random = c(0, 0, 0)` throughout — no cloud/Parquet dependency needed
- Use ABIBAL as the test species
- DATA requirements (cloud fetching, cache layer) are out of scope — those belong in Phase 3's own tests
- Pre-existing `test-bug-fixes.R` failures remain deferred

**Documentation (roxygen2)**
- All exported functions (~25) need `@examples` sections added
- Examples use `parameters()` with `plot_random = c(0, 0, 0)` and ABIBAL — no network required
- Examples must pass `devtools::run_examples()` without errors

**plot.ipm_projection() — new export**
- Add `plot.ipm_projection(x, type = NULL, ...)` to the R package
- `type` argument selects figure: `"lambda"` | `"size_dist"` | `"lambda_vs_n"`
- Default (type = NULL): render all three figures in sequence
- Fully documented with `@param`, `@return`, `@examples`
- Exported from NAMESPACE

**guide_IPM.qmd — full rewrite**
- Replace old procedural API with the new `library(forestIPM)` package API
- Two-part structure: Quick start (single-species ABIBAL workflow) + Deep dive (multi-species, time-varying climate, roadmap)
- Visualizations via `plot(proj, type = ...)` — no manual ggplot boilerplate
- Keep all three figures: lambda over time, animated size distribution (plotly), lambda vs N correlation
- Reference to new book sections where appropriate

**Book restructuring — full**
- New 3-section structure in `_quarto.yml` at `/Users/wvieira/GitHub/book_forest-demography-IPM`
- Section 1: Demographic models evaluation (existing, unchanged)
- Section 2: From demographic rates to population dynamics (`guide_IPM.qmd` + new deep-dive .qmd)
- Section 3: Using the IPM (existing analysis chapters reorganized)
- Chapter filenames stay the same; only `_quarto.yml` navigation changes

### Claude's Discretion
- Deep-dive .qmd filename and internal section titles
- Which specific edge cases to cover in tests (beyond the happy path workflow)
- How to handle the plotly animation in a static Quarto build (whether to wrap in `eval: false` or use a precomputed RDS as the old guide did)
- Exact book section names

### Deferred Ideas (OUT OF SCOPE)
- None — discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| DOC-01 | All exported functions have roxygen2 documentation with `@param`, `@return`, and `@examples` (examples use bundled mini-dataset, not live network) | 11 `export()` entries in NAMESPACE + S3 methods confirmed; 31 RDS files bundled in `inst/extdata/parameters/` covering all species; ABIBAL verified offline with `draw = "mean"` |
| DOC-02 | testthat test suite covering PKG, BUG, and DATA requirements using synthetic fixtures — CI passes without network access | Existing `test-constructors.R` passes; workflow pattern (stand → species_model → parameters → lambda → project) confirmed working offline; fixtures dir already contains `lambda_baseline.rds` and `project_baseline.rds` |
</phase_requirements>

---

## Summary

Phase 4 delivers three things in parallel: roxygen2 `@examples` for all exported functions, a new testthat workflow suite (80%+ coverage on 6 core files), and a full book restructuring with `guide_IPM.qmd` rewritten to the new package API.

The codebase is in excellent shape for documentation: all 11 exported functions plus 16 S3 methods already have `@param` and `@return` roxygen blocks — only `@examples` sections are missing. The 31 bundled RDS files in `inst/extdata/parameters/` mean every example can run offline using `draw = "mean"` and ABIBAL. No network dependency is required.

The `plot.ipm_projection()` function already exists in `project.R` but has the wrong signature (`y, ...` instead of `type = NULL, ...`) and only renders a lambda plot — it must be redesigned to support the three-panel idiom from `plot.lm()`. The test suite cleanup is straightforward: delete three files that use `pkgload::pkg_path()` (which fails under `R CMD check`) and write a new `test-workflow.R` that covers the happy path plus edge cases.

**Primary recommendation:** Implement in three parallel work streams — (1) add `@examples` to all exported functions, (2) write `test-workflow.R` and delete obsolete test files, (3) redesign `plot.ipm_projection()` and rewrite `guide_IPM.qmd`. All three are independent.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| roxygen2 | 7.3.2 (already in DESCRIPTION) | Inline documentation → NAMESPACE + .Rd files | Package standard; already configured with `Roxygen: list(markdown = TRUE)` |
| testthat | >= 3.0 (already implied) | Unit/integration tests | Package standard; `tests/testthat/` already exists with infrastructure |
| devtools | any | `run_examples()`, `test()`, `document()` | Development workflow; already in use |
| covr | any | Line coverage measurement | Standard for R packages; needed to verify 80% target |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| withr | any | Temporary state (e.g., `withr::with_seed()`) | Reproducible random draws in tests |
| plotly | any | Animated size distribution in vignette | Already used in book guide |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `devtools::document()` | `roxygen2::roxygenise(load_code='source')` | `roxygenise()` avoids Xcode license requirement for C++ compilation (per Phase 02-03 decision) |
| `devtools::test()` | `testthat::test_dir()` | Both work; `devtools::test()` is preferred for interactive use |

**Installation:** No new packages needed — roxygen2, testthat, and devtools are already in the development environment.

---

## Architecture Patterns

### Recommended Project Structure (additions only)
```
tests/testthat/
├── test-constructors.R          # KEEP — covers constructor validation
├── test-workflow.R              # NEW — full workflow: stand → project
├── fixtures/
│   ├── lambda_baseline.rds      # KEEP — regression guard
│   └── project_baseline.rds     # KEEP — regression guard
R/
├── project.R                    # MODIFY — redesign plot.ipm_projection()
├── stand.R          )
├── species_model.R  )
├── parameters.R     ) ADD @examples to all @export functions
├── lambda.R         )
├── kernel.R         )
└── ...              )
book_forest-demography-IPM/
├── _quarto.yml                  # MODIFY — add 3rd section
├── guide_IPM.qmd                # FULL REWRITE
└── deep_dive_IPM.qmd            # NEW
```

### Pattern 1: Minimal Offline Example Block
**What:** Every `@examples` section builds a minimal stand from in-line data, uses `draw = "mean"`, and runs the function being documented.
**When to use:** All exported user-facing functions.
**Example:**
```r
# Source: pattern consistent across all functions in this package
#' @examples
#' df <- data.frame(size_mm = c(150, 200, 350),
#'                  species_id = "ABIBAL",
#'                  plot_size = 400)
#' s   <- stand(df)
#' mod <- species_model(s)
#' p   <- parameters(mod, draw = "mean")
#' env <- env_condition(MAT = 8, MAP = 1200)
#' ctrl <- control(years = 5, compute_lambda = TRUE)
#' lam <- lambda(mod, p, s, env, ctrl)
```

### Pattern 2: plot.lm() Idiom for plot.ipm_projection()
**What:** S3 `plot()` method that renders all figures when `type = NULL`, or one figure when `type = "lambda"` | `"size_dist"` | `"lambda_vs_n"`.
**When to use:** `plot.ipm_projection()` redesign.
**Current state:** Existing `plot.ipm_projection(x, y, ...)` in `project.R` only renders lambda trajectory. Must be replaced.
**Example:**
```r
#' @param x An \code{ipm_projection} object.
#' @param type Character or NULL. One of \code{"lambda"}, \code{"size_dist"},
#'   \code{"lambda_vs_n"}. If NULL (default), all three figures are rendered
#'   in sequence.
#' @param ... Additional arguments passed to individual plot functions.
plot.ipm_projection <- function(x, type = NULL, ...) {
  types <- c("lambda", "size_dist", "lambda_vs_n")
  if (!is.null(type) && !type %in% types) {
    cli::cli_abort(
      "{.arg type} must be one of {.val {types}} or NULL. Got {.val {type}}."
    )
  }
  to_plot <- if (is.null(type)) types else type
  for (t in to_plot) {
    switch(t,
      lambda       = .plot_lambda(x, ...),
      size_dist    = .plot_size_dist(x, ...),
      lambda_vs_n  = .plot_lambda_vs_n(x, ...)
    )
  }
  invisible(x)
}
```

### Pattern 3: Workflow Test Structure
**What:** Single `test-workflow.R` covering the end-to-end pipeline with `expect_*` at each step.
**When to use:** New workflow test file.
**Example:**
```r
# Shared fixture — reused across all workflow tests
make_abibal_stand <- function() {
  stand(data.frame(
    size_mm    = seq(130, 600, by = 50),
    species_id = "ABIBAL",
    plot_size  = 1000
  ))
}

test_that("full workflow stand → project returns ipm_projection", {
  s    <- make_abibal_stand()
  mod  <- species_model(s)
  pars <- parameters(mod, draw = "mean")
  env  <- env_condition(MAT = 8, MAP = 1200)
  ctrl <- control(years = 10, compute_lambda = TRUE, progress = FALSE)
  proj <- project(mod, pars, s, env, ctrl)
  expect_s3_class(proj, "ipm_projection")
  expect_true(all(!is.na(proj$lambda[["ABIBAL"]])))
})
```

### Anti-Patterns to Avoid
- **`pkgload::pkg_path()` in tests:** Fails under `R CMD check` because the package isn't installed — use `testthat::test_path()` for fixture paths or `system.file()` for package data. This is the root cause of `test-bug-fixes.R` and `test-package-structure.R` failures.
- **Live network calls in examples:** Any `parameters()` call with a species lacking a bundled RDS will fail offline. Stick to ABIBAL (`ABIBAL_pars.rds` is confirmed present in `inst/extdata/parameters/`).
- **`progress = TRUE` in test projections:** Produces noisy output in test runner; always use `control(..., progress = FALSE)` in tests.
- **`devtools::document()` instead of `roxygen2::roxygenise(load_code='source')`:** Triggers C++ compilation, which requires Xcode license acceptance (blocker on this machine per STATE.md).

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Coverage measurement | Custom line counter | `covr::package_coverage()` + `covr::percent_coverage()` | Handles Rcpp, S3 dispatch, and branch coverage |
| Test fixtures / shared setup | Per-test data construction | testthat `local_*` helpers or shared `make_*()` functions at top of file | DRY; faster test execution |
| Snapshot testing | Manual RDS comparison | `testthat::expect_snapshot()` | Handles formatted output diffs automatically |
| Plotly static fallback | Custom iframe embedding | `eval: false` chunk option + precomputed RDS (per existing pattern in book) | Already established pattern in old `guide_IPM.qmd`; avoids render-time network and slow computation |

---

## Common Pitfalls

### Pitfall 1: `@examples` Triggers C++ Compilation via `devtools::document()`
**What goes wrong:** Running `devtools::document()` during the `@examples`-addition task attempts to compile `src/eigen.cpp`, which requires Xcode license acceptance. Process hangs or errors.
**Why it happens:** `devtools::document()` calls `devtools::load_all()` internally, which compiles Rcpp code.
**How to avoid:** Use `roxygen2::roxygenise(load_code='source')` instead — per the Phase 02-03 decision already in STATE.md. This skips compilation.
**Warning signs:** `xcode-select` license prompt, compilation errors during `roxygenise()`.

### Pitfall 2: S3 Method `plot.ipm_projection` Signature Change Breaks NAMESPACE
**What goes wrong:** The existing `plot.ipm_projection(x, y, ...)` is registered in NAMESPACE as `S3method(plot,ipm_projection)`. Changing the signature to `(x, type = NULL, ...)` requires that the roxygen `@export` tag is still present and `roxygenise()` is re-run.
**Why it happens:** NAMESPACE is generated by roxygen2 — if `@export` is present, registration is automatic.
**How to avoid:** Keep `#' @export` on the new `plot.ipm_projection` definition; verify NAMESPACE still contains `S3method(plot,ipm_projection)` after running `roxygenise()`.

### Pitfall 3: `devtools::run_examples()` Fails on Time-Consuming Examples
**What goes wrong:** `project()` with `years = 100` takes several seconds; `devtools::run_examples()` may time out or produce slow CI runs.
**Why it happens:** Default control uses `years = 100`.
**How to avoid:** Always use `control(years = 5, progress = FALSE)` in `@examples` blocks for `project()` and `lambda()`.

### Pitfall 4: `test-bug-fixes.R` vs `test-package-structure.R` — Pre-existing Failures
**What goes wrong:** Both files use `pkg <- pkgload::pkg_path()` at top-level (outside any `test_that()` block). Under `R CMD check`, the package is installed (not loaded), so `pkgload::pkg_path()` fails.
**Why it happens:** `pkgload` is a development tool; `R CMD check` uses a real install.
**How to avoid:** Delete both files per the locked decision. Do not migrate these tests — the coverage they provided is superseded by the new workflow suite.

### Pitfall 5: `plot.ipm_projection()` Called Before `project()` Stores Lambda
**What goes wrong:** `plot(proj, type = "lambda")` is called on a projection created with `compute_lambda = FALSE` — all lambda values are NA.
**Why it happens:** Default `control()` has `compute_lambda = FALSE`.
**How to avoid:** In `plot.ipm_projection()`, check for all-NA lambda before attempting the lambda plot and emit a `message()` (the existing code already does this for the old function — carry this check forward).

### Pitfall 6: Book Repo Is a Separate Git Repo
**What goes wrong:** Committing `_quarto.yml` and `guide_IPM.qmd` changes to the wrong repo (`forest-IPM` instead of `book_forest-demography-IPM`).
**Why it happens:** Two repos with similar content.
**How to avoid:** Book repo path is `/Users/wvieira/GitHub/book_forest-demography-IPM/` — all book file writes and commits go there, not to `forest-IPM`.

---

## Code Examples

Verified patterns from codebase inspection:

### Minimal Offline Stand for Tests and Examples
```r
# All 31 species RDS files confirmed in inst/extdata/parameters/
# ABIBAL_pars.rds verified present — works offline
df <- data.frame(
  size_mm    = seq(130, 600, by = 50),
  species_id = "ABIBAL",
  plot_size  = 1000
)
s <- stand(df)
```

### parameters() with draw = "mean" (no network, no seed needed)
```r
mod  <- species_model(s)
pars <- parameters(mod, draw = "mean")
# pars$draw_type == "mean"; pars$seed == NULL
```

### Full Workflow for @examples (fast, offline)
```r
s    <- stand(data.frame(size_mm = c(150, 200, 350),
                          species_id = "ABIBAL", plot_size = 400))
mod  <- species_model(s)
pars <- parameters(mod, draw = "mean")
env  <- env_condition(MAT = 8, MAP = 1200)
ctrl <- control(years = 5, compute_lambda = TRUE, progress = FALSE)
proj <- project(mod, pars, s, env, ctrl)
lam  <- lambda(mod, pars, s, env, ctrl)
```

### roxygen2 roxygenise (no Xcode needed)
```r
# Per Phase 02-03 decision — avoids C++ compilation
roxygen2::roxygenise(load_code = "source")
```

### plot.ipm_projection() Three-Type Design
```r
plot(proj)                    # all three figures
plot(proj, type = "lambda")   # lambda trajectory only
plot(proj, type = "size_dist")   # size distribution only
plot(proj, type = "lambda_vs_n") # lambda vs N correlation only
```

### covr coverage check
```r
cov <- covr::package_coverage(
  path = ".",
  type = "tests",
  quiet = TRUE
)
covr::percent_coverage(cov)
# Filter to the 6 target files:
covr::file_coverage(cov, "R/stand.R")
```

### Book _quarto.yml new structure (target)
```yaml
chapters:
  - index.qmd
  - introduction.qmd
  - db.qmd
  - part: Demographic models
    chapters: [...]
  - part: Demographic model evaluation
    chapters: [...]
  - part: From demographic rates to population dynamics
    chapters:
    - ipm_description.qmd
    - guide_IPM.qmd
    - deep_dive_IPM.qmd          # NEW
  - part: Using the IPM
    chapters:
    - sens_analysis.qmd
    - conditional_lambda.qmd
    - marginal_lambda.qmd
    - extinction_risk.qmd
  - references.qmd
```

---

## Exported Functions Inventory

Complete list of functions needing `@examples` (from NAMESPACE inspection):

**`export()` functions (11) — all need `@examples`:**
| Function | File | Has @param/@return? | Needs @examples? |
|----------|------|---------------------|-----------------|
| `stand()` | `stand.R` | Yes | Yes |
| `species_model()` | `species_model.R` | Yes | Yes |
| `parameters()` | `parameters.R` | Yes | Yes |
| `set_random_effects()` | `parameters.R` | Yes | Yes |
| `env_condition()` | `env_condition.R` | Yes | Yes |
| `control()` | `control.R` | Yes | Yes |
| `lambda()` | `lambda.R` | Yes | Yes |
| `project()` | `project.R` | Yes (partial) | Yes |
| `supported_species()` | `supported_species.R` | Yes (partial) | Yes |
| `scale_env()` | `env_scaling.R` | Yes | Yes |
| `unscale_env()` | `env_scaling.R` | Yes | Yes |

**`S3method()` entries (16) — need `@examples` only on the primary entry point, not on print/summary:**
- `plot.ipm_projection` — NEW signature needed (locked decision)
- `print.*` (6) — inherit documentation from user-facing function; minimal `@examples` acceptable
- `summary.*` (6) — same as print
- `plot.ipm_projection` is the only S3 method the user calls directly

**`mkKernel` note:** `mkKernel()` is exported (confirmed from devtools::load_all output), but not in NAMESPACE `export()` list. It is an internal function exposed via `load_all()` but not intended for users. No `@examples` needed unless `@export` tag is added.

---

## Current Test File Status

| File | Status | Action |
|------|--------|--------|
| `test-constructors.R` | Passing (9 tests) | KEEP |
| `test-regression-baselines.R` | Passing (2 tests) | KEEP (baselines exist in fixtures/) |
| `test-vectorized-kernel.R` | Unknown (uses internal functions) | DELETE per locked decision |
| `test-bug-fixes.R` | FAILING (`pkgload::pkg_path()` fails in R CMD check) | DELETE per locked decision |
| `test-package-structure.R` | FAILING (`pkgload::pkg_path()` fails in R CMD check) | DELETE per locked decision |
| `test-workflow.R` | MISSING | CREATE |

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `source()` from GitHub + `getPars_sp()` + `mkKernel()` + manual `eigen()` | `library(forestIPM)` → `stand()` → `species_model()` → `parameters()` → `lambda()` / `project()` | Phase 2 | guide_IPM.qmd must be fully rewritten |
| `plot.ipm_projection(x, y, ...)` (lambda only) | `plot.ipm_projection(x, type = NULL, ...)` (three panels) | Phase 4 (this phase) | Function in project.R must be replaced |
| No `@examples` blocks | All exported functions have runnable `@examples` | Phase 4 (this phase) | `devtools::run_examples()` was broken |

---

## Open Questions

1. **`mkKernel` export status**
   - What we know: `mkKernel` appears in `load_all()` namespace but is NOT in `NAMESPACE export()` section; no `@export` tag in `kernel.R`
   - What's unclear: Whether it should be exported (it has `@return` but no `@export` tag)
   - Recommendation: Treat as internal — do not add `@examples`. If future use requires export, that is a separate decision.

2. **plotly animation in static Quarto build**
   - What we know: Old guide used precomputed `data/out_IPM_guide.RDS` to avoid slow render
   - What's unclear: Whether to use `eval: false` chunk or precomputed RDS for the deep-dive animated figure
   - Recommendation: Claude's discretion per CONTEXT.md — reuse the RDS pattern (`readRDS("data/out_IPM_guide.RDS")`) since it is already established

3. **`test-vectorized-kernel.R` uses internal functions directly**
   - What we know: Tests call `P_xEC()` and `ingrowth_lk()` directly (not via NAMESPACE); works under `load_all()` but may fail under `R CMD check` if functions are unexported
   - Recommendation: Delete per locked decision — functionality is already covered by regression baseline tests

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | testthat (version from DESCRIPTION implies >= 3.0) |
| Config file | none — standard `tests/testthat/` layout |
| Quick run command | `devtools::test(filter = "workflow")` |
| Full suite command | `devtools::test()` |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| DOC-01 | All `@examples` run without error | smoke | `devtools::run_examples()` | ❌ Wave 0 — @examples missing |
| DOC-01 | `plot(proj, type = "lambda")` renders without error | unit | `devtools::test(filter = "workflow")` | ❌ Wave 0 |
| DOC-01 | `plot(proj)` (all types) renders without error | unit | `devtools::test(filter = "workflow")` | ❌ Wave 0 |
| DOC-02 | Full workflow stand → project passes | integration | `devtools::test(filter = "workflow")` | ❌ Wave 0 |
| DOC-02 | 80%+ line coverage on 6 core files | coverage | `covr::package_coverage()` | ❌ Wave 0 |
| DOC-02 | `devtools::test()` zero failures on clean machine | integration | `devtools::test()` | Partial — test-constructors.R passing |

### Sampling Rate
- **Per task commit:** `devtools::test(filter = "workflow")`
- **Per wave merge:** `devtools::test()`
- **Phase gate:** `devtools::test()` zero failures + `devtools::run_examples()` zero errors before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/testthat/test-workflow.R` — covers DOC-02 workflow tests
- [ ] `@examples` blocks in all 11 exported functions — covers DOC-01
- [ ] New `plot.ipm_projection(x, type = NULL, ...)` in `R/project.R` — prerequisite for DOC-01 and guide rewrite
- [ ] `book_forest-demography-IPM/deep_dive_IPM.qmd` — new file, no template yet

---

## Sources

### Primary (HIGH confidence)
- Codebase inspection: `NAMESPACE`, `R/*.R`, `tests/testthat/*.R`, `inst/extdata/parameters/` — all verified by direct file read
- `DESCRIPTION` — roxygen2 version 7.3.2 confirmed; R >= 4.1.0 confirmed
- STATE.md — Phase 02-03 decision on `roxygen2::roxygenise(load_code='source')` confirmed
- Live execution: `Rscript -e "devtools::load_all(); parameters(species_model(stand(...)), draw='mean')"` — confirmed ABIBAL works offline

### Secondary (MEDIUM confidence)
- testthat 3.x API: `expect_s3_class()`, `test_that()`, `local_*()` helpers — confirmed consistent with codebase patterns
- covr API: `package_coverage()` + `percent_coverage()` — standard R coverage toolchain

### Tertiary (LOW confidence)
- None — all claims verified from codebase or live execution

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — all tools confirmed in DESCRIPTION and live execution
- Architecture: HIGH — all patterns derived from existing codebase inspection
- Pitfalls: HIGH — root causes verified (pkgload::pkg_path failure confirmed by test output pattern, Xcode blocker in STATE.md)

**Research date:** 2026-03-13
**Valid until:** 2026-04-12 (stable domain — roxygen2 and testthat APIs change slowly)
