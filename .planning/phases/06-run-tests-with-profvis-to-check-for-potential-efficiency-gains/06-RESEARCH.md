# Phase 6: Run tests with profvis to check for potential efficiency gains - Research

**Researched:** 2026-03-11
**Domain:** R performance profiling, tidyverse style, Rcpp C++ migration
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**Sub-step 1: Tidyverse-style revamp**
- Goal: Consistency and readability — make all R/ code consistently idiomatic (pipes, purrr::map, dplyr verbs). Performance improvement is a secondary benefit, not the primary goal.
- Scope: All R/ files — params.R, project.R, lambda.R, kernel.R, stand.R, BasalArea_competition.R, vital_rates.R, and all other files in R/. Full pass for consistency.
- Boundary: R/ source files only — do not touch test files or roxygen2 @examples (those belong to Phase 4).
- Skill to use: `r-skills:tidyverse-patterns`
- Verification: `devtools::check()` (R CMD check) + existing testthat suite must pass after revamp.

**Sub-step 2: profvis profiling**
- Benchmark inputs:
  - Single-species run: 1 species (QUERUB — large max size stresses big kernel), 100 timesteps via `project()`
  - Multi-species run: 2-3 species including QUERUB and QUEPRI, 100 timesteps — exposes competition update costs
  - Rationale: single vs multi-species expose different code paths
- Skill to use: `r-skills:r-performance`
- Output deliverables:
  - A written markdown summary (`PROFILING.md`) in the phase directory documenting: top hotspots, call stack breakdown, time percentages, and recommended optimization targets for sub-step 3
  - An HTML profvis flamegraph saved to the phase directory as an artifact
- Script location: `.planning/phases/06-.../` (not shipped with the package)

**Sub-step 3: Targeted optimizations**
- Decision framework: Apply if measurable speedup (>10-15% on benchmark) AND implementation is low-risk (no change to math/statistical outputs). Document anything evaluated but skipped and why.
- C++ bias: C++ migrations (via Rcpp) should be weighted favorably — the package already depends on RcppEigen, so adding more C++ functions is low dependency overhead.
- Strategies in scope (all conditional on profvis evidence):
  - Vectorize `outer()` calls in `mkKernel()` — replace `outer(meshpts, meshpts, FUN)` with pre-allocated matrix + vectorized operations to avoid per-cell R function call overhead
  - Cache competition metrics — `BA_comp_intra/inter` are recomputed inside `mkKernel()` every call; cache when Nvec hasn't changed between timesteps
  - C++ migration of inner loops — move `P_xEC` or `ingrowth_lk` to C++ via Rcpp if they appear as hotspots in profvis flamegraph
- Final check: Run benchmark comparison (before/after) to quantify gains from each applied optimization.

### Claude's Discretion
- Exact profvis flamegraph presentation format (HTML styling, chart options)
- Which specific dplyr/purrr patterns to use in each file (follow tidyverse-patterns skill guidance)
- How to structure the PROFILING.md report
- Whether to combine optimizations into one commit or separate commits per optimization

### Deferred Ideas (OUT OF SCOPE)
- None — discussion stayed within phase scope
</user_constraints>

---

## Summary

Phase 6 has three sequential sub-steps: a tidyverse-style revamp of all `R/` source files, profvis-based profiling of the IPM computation pipeline, and targeted code optimizations based on profiling evidence. The dominant computational hotspot is `mkKernel()`, which calls `outer(meshpts, meshpts, FUN)` twice per species per timestep — this pattern invokes an R function (either `P_xEC` or `ingrowth_lk`) once for every cell of an n×n matrix, creating n² R-level function calls per kernel build. For QUERUB with its large maximum DBH, the mesh is correspondingly wide, amplifying this cost. The `project()` loop rebuilds the kernel every timestep × every species, so a 100-timestep, 2-species run calls `mkKernel()` 200 times.

The profiling script itself is a standalone `.R` file in the phase directory, not shipped with the package. It uses `profvis()` to wrap both `lambda()` and `project()` calls, then exports the flamegraph as HTML using `htmlwidgets::saveWidget()`. The profvis output drives a written `PROFILING.md` summary. Optimizations are conditional on the profvis evidence: vectorizing `outer()` calls, caching competition metrics across unchanged timesteps, and optionally migrating `P_xEC` or `ingrowth_lk` to C++ via Rcpp are the three candidate strategies. A before/after benchmark (using `microbenchmark` or `bench`) quantifies each optimization's gain.

**Primary recommendation:** Run profvis on the 100-timestep QUERUB benchmark first; the flamegraph will confirm whether `outer()` call overhead or eigenvalue computation dominates. Apply vectorization of `outer()` calls as the first optimization (high expected gain, deterministic), then evaluate the others based on actual profile data.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| profvis | >= 0.3.8 | Interactive flame graph profiler; wraps Rprof() with visualization | The canonical R profiling tool; used by Advanced R, RStudio, and Tidyverse teams |
| htmlwidgets | >= 1.6.0 | Save profvis HTML output to file | Required dependency for `htmlwidgets::saveWidget()` — transitive dep of profvis |
| microbenchmark | >= 1.4-7 | Micro-benchmark before/after comparisons; 100 reps by default | Standard benchmarking pair alongside profvis |
| bench | >= 1.1.3 | Alternative benchmark with memory tracking and ggplot output | Better than microbenchmark for memory-sensitive comparisons |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Rcpp | >= 1.0.10 | Already a package dependency — C++ function registration | When migrating `P_xEC` or `ingrowth_lk` inner loops to C++ |
| RcppEigen | already linked | Eigen matrix operations from C++ | Existing `src/eigen.cpp` is the template; use for any new matrix C++ functions |
| truncnorm | >= 1.0-9 | Already a package dependency — truncated normal density | Used in `ingrowth_lk()` via `dtruncnorm()`; if migrating to C++, this call must be replicated (use Eigen + custom truncated normal, or call back to R) |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| profvis | Rprof() directly | profvis wraps Rprof() with interactive HTML output — always prefer profvis for visual inspection |
| microbenchmark | system.time() | system.time() is single-run; microbenchmark provides distribution across 100 runs, revealing variance |
| bench | microbenchmark | bench tracks memory allocation; prefer bench when GC overhead is suspected |

**Installation (for profiling script only — not added to DESCRIPTION):**
```r
install.packages(c("profvis", "microbenchmark", "bench"))
```

---

## Architecture Patterns

### Sub-step 1: Tidyverse-style revamp

#### Files in scope (all `R/` source files)
```
R/
├── BasalArea_competition.R   # sapply() -> purrr::map_dbl() or vectorized
├── control.R                 # Minimal; check for base pipe opportunities
├── env_condition.R           # Check for tidyverse idioms
├── env_scaling.R             # Check for tidyverse idioms
├── globals.R                 # @importFrom declarations — add new symbols here
├── kernel.R                  # outer() calls, mixed style
├── lambda.R                  # vapply() -> purrr::map_dbl()
├── parameters.R              # Already partially tidyverse; full pass
├── project.R                 # for loops -> functional where applicable
├── RcppExports.R             # Auto-generated; do NOT touch
├── species_model.R           # Review for style consistency
├── stand.R                   # Review for style consistency
├── supported_species.R       # Review for style consistency
└── vital_rates.R             # Pure math functions; style pass for readability
```

#### Pattern 1: sapply/vapply -> purrr::map family
**What:** Replace base R iteration with typed purrr equivalents
**When to use:** Any loop producing a vector or list from a list input
**Example:**
```r
# Before (BasalArea_competition.R)
sapply(N_intra$meshpts, function(x) sum(ba_N_plot[N_focal$meshpts > x]))

# After
purrr::map_dbl(N_intra$meshpts, \(x) sum(ba_N_plot[N_focal$meshpts > x]))
```

#### Pattern 2: Native pipe |> over magrittr %>%
**What:** Replace `%>%` with `|>` (R >= 4.1.0 — this package requires R >= 4.1.0)
**When to use:** All pipeline expressions in R/ files
**Example:**
```r
# Before
sp_data$draws %>% dplyr::filter(.data$draw == draw_id) %>% ...

# After
sp_data$draws |> dplyr::filter(.data$draw == draw_id) |> ...
```

#### Pattern 3: Bare .data$ pronoun consistency
**What:** Ensure all dplyr column references use `.data$col` to avoid R CMD check NOTEs
**When to use:** Any `dplyr::filter()`, `dplyr::select()`, `dplyr::mutate()` call referencing bare column names

#### Anti-Patterns to Avoid
- **Touching RcppExports.R:** Auto-generated by `Rcpp::compileAttributes()` — never edit manually
- **Changing function signatures:** Public API is locked from Phase 1 — style revamp changes internals only
- **Over-purring tight numeric loops:** `purrr::map_dbl()` has some overhead versus `vapply()`; in `size_to_BAcomp()` the `sapply()` is called once per kernel build (not per matrix cell), so the gain is readability, not speed
- **Adding new @importFrom without updating globals.R:** Any new purrr/dplyr symbols introduced must be declared in `globals.R`

### Sub-step 2: profvis profiling

#### Profiling script structure
```
.planning/phases/06-.../
├── profile_benchmark.R     # standalone script (not shipped with package)
├── PROFILING.md            # written summary of hotspots (deliverable)
└── profvis_output.html     # saved flamegraph artifact (deliverable)
```

#### Pattern: profvis() wrapping
**What:** Wrap the benchmark code in `profvis({...})`, capture the result, save to HTML
**Example:**
```r
# Source: https://profvis.r-lib.org/articles/profvis.html
library(forestIPM)
library(profvis)

# Build benchmark inputs (synthetic QUERUB stand)
df_single <- data.frame(
  size_mm    = seq(200, 800, by = 50),
  species_id = "QUERUB",
  plot_size  = 1000
)
s    <- stand(df_single)
mod  <- species_model(s)
pars <- parameters(mod, draw = "mean")
env  <- env_condition(MAT = 10, MAP = 1100)
ctrl <- control(years = 100, progress = FALSE)

# Profile project()
p_project <- profvis({
  project(mod, pars, s, env, ctrl)
})
htmlwidgets::saveWidget(p_project, "profvis_output.html", selfcontained = TRUE)

# Profile lambda() separately (single kernel build path)
p_lambda <- profvis({
  for (i in 1:200) lambda(mod, pars, s, env, ctrl)
})
```

#### Pattern: benchmark comparison (before/after)
```r
library(microbenchmark)

mb <- microbenchmark(
  original  = project(mod, pars, s, env, ctrl),
  optimized = project_optimized(mod, pars, s, env, ctrl),
  times = 20
)
summary(mb)
```

#### Key profvis parameters
| Parameter | Recommended Value | Reason |
|-----------|------------------|--------|
| `interval` | 0.005 (5ms) | Default 10ms may undersample fast functions; 5ms is minimum reliable value |
| `height` | "600px" | Wider view for deep call stacks from `mkKernel()` -> `outer()` -> `P_xEC` |

### Sub-step 3: Targeted optimizations

#### Optimization 1: Vectorize outer() calls in mkKernel()
**What:** The two `outer(meshpts, meshpts, FUN)` calls in `mkKernel()` invoke an R function once per matrix cell. For n mesh points, this is n² R function calls per kernel build. Replacing with vectorized computation avoids per-cell call overhead.

**How:** Pre-expand `meshpts` into all (size_t0, size_t1) pairs using `rep()` and `rep()` with `each=`, evaluate `P_xEC` / `ingrowth_lk` on the full vectors at once (since all vital rate functions are already vectorized over scalars), then `matrix()` to reshape.

**Example (conceptual):**
```r
# Current (n^2 R function calls)
P <- h * outer(meshpts, meshpts, P_xEC, ...)

# Vectorized alternative (1 R function call on length-n^2 vectors)
n <- length(meshpts)
s0 <- rep(meshpts, each = n)   # size_t0 varies slowly
s1 <- rep(meshpts, times = n)  # size_t1 varies fast
P_vals <- P_xEC(s1, s0, ...)   # all vital rate functions are vectorized
P <- h * matrix(P_vals, nrow = n, ncol = n)
```

**Prerequisite verification:** Confirm `vonBertalanffy_lk()`, `survival_f()`, `vonBertalanffy_f()` are all vectorized over `size_t0` and `size_t1`. Current code uses `pars['r']`, `exp()`, arithmetic — all base R vectorized. `dnorm()` is vectorized. **This approach is safe.**

**IMPORTANT: `outer()` argument order.** In R, `outer(X, Y, FUN)` calls `FUN(X_expanded, Y_expanded)`. In `mkKernel()`, the call is `outer(meshpts, meshpts, P_xEC, ...)` — `P_xEC`'s first argument is `size_t1` and second is `size_t0`. In the vectorized replacement: `s1 <- rep(meshpts, each = n)` (outer's X maps to FUN's first arg) and `s0 <- rep(meshpts, times = n)`. Verify this mapping matches the existing `outer()` semantics before replacing.

#### Optimization 2: Cache competition metrics across timesteps
**What:** `BA_comp_intra`, `BA_comp_inter`, `BAplot_intra`, `BAplot_inter` are recomputed inside `mkKernel()` every call. Each involves a `sapply()` over meshpts (an O(n) operation called n times inside `size_to_BAcomp()`). These values change only when `Nvec` changes.

**How:** In `project()`, cache the BA metrics after each timestep's `mkKernel()` call and pass them in as arguments rather than recomputing. This requires adding optional `BA_comp_intra`, `BA_comp_inter` arguments to `mkKernel()` — only compute them if NULL.

**Risk:** Moderate — requires changing `mkKernel()` signature (internal function, not exported). Must verify caching logic is correct for multi-species case where `Nvec_inter` changes every timestep.

#### Optimization 3: C++ migration of P_xEC / ingrowth_lk
**What:** If `outer()` vectorization is not pursued or still leaves hotspots, migrating the per-cell evaluation function to C++ eliminates R function call overhead entirely. The existing `src/eigen.cpp` is the template.

**Integration pattern (using existing template):**
```cpp
// src/kernel_ops.cpp
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd buildPkernel(
  Eigen::VectorXd meshpts, double delta_time,
  double BA_comp_intra, double BA_comp_inter,
  double Temp, double Prec,
  Eigen::VectorXd parsGrowth, Eigen::VectorXd parsMort,
  double plot_random_growth, double plot_random_mort
) {
  // ... vectorized C++ implementation
}
```

After adding, run `Rcpp::compileAttributes()` to regenerate `R/RcppExports.R`.

**Note on truncnorm:** `ingrowth_lk()` calls `truncnorm::dtruncnorm()`. Migrating `ingrowth_lk` to C++ requires either: (a) calling R from C++ via `Rcpp::Function`, which partially negates the benefit, or (b) implementing truncated normal density in C++ (straightforward: `dnorm(x)/pnorm_complement`). Either approach is viable.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| R code profiling | Custom timer with `proc.time()` around sections | `profvis()` | profvis samples the full call stack at intervals — custom timers miss nested call costs and GC |
| Saving HTML output | Custom HTML templating | `htmlwidgets::saveWidget(profvis_result, "file.html", selfcontained = TRUE)` | profvis returns an htmlwidget; saveWidget handles all serialization correctly |
| Benchmarking before/after | `system.time()` single run | `microbenchmark()` or `bench::mark()` | Single runs are dominated by noise; 20-100 reps give reliable median and IQR |
| Vectorized outer | Custom nested for-loop | `rep()` + direct vectorized evaluation + `matrix()` reshape | The rep/matrix pattern is idiomatic R; nested for-loops are much slower than even `outer()` |
| C++ compilation | Manual gcc calls | `devtools::load_all()` after editing `src/` + `Rcpp::compileAttributes()` | Rcpp build system handles include paths, linking, attribute registration automatically |

---

## Common Pitfalls

### Pitfall 1: outer() argument order mismatch during vectorization
**What goes wrong:** When replacing `outer(meshpts, meshpts, P_xEC, ...)` with the `rep()` + vectorized pattern, the expansion order of `rep(meshpts, each = n)` vs `rep(meshpts, times = n)` must match `outer()`'s own expansion: outer's X is the first arg to FUN, outer's Y is the second.
**Why it happens:** The distinction between `each=` and `times=` is easily reversed, silently producing transposed kernel matrices with wrong math.
**How to avoid:** Add a unit test that builds a tiny (3×3) kernel with both approaches and checks element-wise equality against the original `outer()` result before replacing.
**Warning signs:** Lambda values that change significantly after vectorization refactor.

### Pitfall 2: R CMD check fails after adding new tidyverse symbols
**What goes wrong:** Using a new purrr or dplyr function (e.g., `purrr::list_rbind()`, `dplyr::pick()`) without adding an `@importFrom` declaration causes "no visible global function definition" NOTE in R CMD check.
**Why it happens:** `globals.R` only covers symbols currently in use; new symbols from the tidyverse revamp need explicit declaration.
**How to avoid:** After the tidyverse revamp, run `devtools::check()` and scan for NOTE lines about undeclared symbols. Add them to `@importFrom` in `globals.R`.
**Warning signs:** `devtools::check()` output with "no visible global function definition for 'X'" NOTEs.

### Pitfall 3: profvis undersampling fast functions
**What goes wrong:** With default `interval = 0.01` (10ms), functions that complete in < 10ms may appear as zero-time in the flamegraph, even if they are called millions of times and collectively dominate runtime.
**Why it happens:** The sampling profiler only captures the call stack at each interval tick; fast functions that start and finish between ticks are invisible.
**How to avoid:** Use `interval = 0.005` (5ms). Also run the benchmark for enough iterations (100+ timesteps) to give each function enough cumulative time to be sampled reliably.
**Warning signs:** Flat flamegraph that doesn't match intuition about where time should be spent.

### Pitfall 4: Caching competition metrics breaks multi-species correctness
**What goes wrong:** If BA metrics are cached across timesteps in `project()`, but the cache is not invalidated when `Nvec_inter` (the competitor's size distribution) changes, the focal species will use stale competition values.
**Why it happens:** In the multi-species loop, each species' `Nvec_inter` is updated each timestep via `.update_N_het()` — the cache must be per-species and per-timestep for the inter-species component.
**How to avoid:** Only cache `BA_comp_intra` (which depends only on focal species `Nvec`, updated after the kernel is built). Recompute `BA_comp_inter` every timestep unless doing single-species runs.
**Warning signs:** Multi-species lambda values diverge from single-species runs in ways not explained by competition.

### Pitfall 5: Rcpp::compileAttributes() not run after adding C++ files
**What goes wrong:** New C++ functions with `// [[Rcpp::export]]` are not available in R until `Rcpp::compileAttributes()` regenerates `R/RcppExports.R` and `src/RcppExports.cpp`.
**Why it happens:** `devtools::load_all()` compiles C++ but does NOT automatically regenerate the attribute registration files.
**How to avoid:** Run `Rcpp::compileAttributes()` before `devtools::load_all()` whenever adding or modifying `// [[Rcpp::export]]` function signatures.
**Warning signs:** `devtools::load_all()` succeeds but the new function is not found in R.

### Pitfall 6: profvis script accidentally shipped in package build
**What goes wrong:** The profiling script in `.planning/phases/06-.../` gets included in the package tarball, bloating it.
**Why it happens:** `.planning/` is not listed in `.Rbuildignore`.
**How to avoid:** Verify `.Rbuildignore` contains `^\.planning$` (it was added in Phase 2 plan 02-03). If not, add it.
**Warning signs:** `devtools::build()` produces a tarball larger than expected.

---

## Code Examples

### profvis() basic profiling pattern
```r
# Source: https://profvis.r-lib.org/articles/profvis.html
library(profvis)
p <- profvis({
  # code to profile
}, interval = 0.005, height = "600px")
htmlwidgets::saveWidget(p, "profvis_output.html", selfcontained = TRUE)
```

### Outer-to-vectorized replacement pattern (conceptual)
```r
# Original: n^2 R function calls
P <- h * outer(meshpts, meshpts, P_xEC,
               delta_time, BA_comp_intra, BA_comp_inter,
               Temp, Prec, parsGrowth, parsMort, plot_random)

# Vectorized: 1 call on length-n^2 vectors
n   <- length(meshpts)
s1  <- rep(meshpts, each = n)   # outer's X -> FUN's first arg (size_t1)
s0  <- rep(meshpts, times = n)  # outer's Y -> FUN's second arg (size_t0)
P_vals <- P_xEC(s1, s0, delta_time, BA_comp_intra, BA_comp_inter,
                Temp, Prec, parsGrowth, parsMort, plot_random)
P   <- h * matrix(P_vals, nrow = n, ncol = n)
```

### microbenchmark before/after comparison
```r
# Source: https://adv-r.hadley.nz/perf-measure.html pattern
library(microbenchmark)
mb <- microbenchmark(
  original   = { mkKernel(Nvec_intra, Nvec_inter, dt, ps, T, P, pars, pr) },
  vectorized = { mkKernel_v2(Nvec_intra, Nvec_inter, dt, ps, T, P, pars, pr) },
  times = 50
)
print(mb)
boxplot(mb)
```

### New Rcpp function registration flow
```r
# After editing src/kernel_ops.cpp:
Rcpp::compileAttributes()   # regenerates R/RcppExports.R
devtools::load_all()        # compiles and loads
devtools::check()           # verify 0 errors, 0 warnings
```

### sapply -> purrr::map_dbl (BasalArea_competition.R)
```r
# Before
sapply(N_intra$meshpts, function(x) sum(ba_N_plot[N_focal$meshpts > x]))

# After (tidyverse revamp)
purrr::map_dbl(N_intra$meshpts, \(x) sum(ba_N_plot[N_focal$meshpts > x]))
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `%>%` magrittr pipe | `|>` native pipe | R 4.1.0 (2021) | No magrittr dependency for piping; slightly faster due to no function wrapping |
| `purrr::map_dfr()` | `purrr::map() |> list_rbind()` | purrr 1.0.0 (2022) | map_dfr() is deprecated; list_rbind() is the current idiom |
| `purrr::as_vector()` | `unlist(use.names = TRUE)` | purrr 1.0.0 (2022) | as_vector() removed in purrr >= 1.0.0 (already fixed in Phase 2) |
| `group_by() + summarise() + ungroup()` | `summarise(.by = x)` | dplyr 1.1.0 (2023) | Transient grouping avoids the ungroup() footgun |
| `Vectorize()` wrapper | Direct vectorized arithmetic or `purrr::map_dbl()` | Style evolution | Vectorize() has higher overhead than direct vectorization |

**Deprecated/outdated:**
- `purrr::as_vector()`: Removed in purrr >= 1.0.0; replaced with `unlist()`. Already fixed in Phase 2.
- `purrr::map_dfr()` / `map_dfc()`: Deprecated in purrr 1.0.0; use `map() |> list_rbind()` / `list_cbind()`.
- `magrittr::%>%`: Not deprecated but `|>` is now preferred for new R code.

---

## Open Questions

1. **Mesh size for QUERUB**
   - What we know: QUERUB (Quercus rubra) has large maximum DBH; `bin_width = 1` (default in `control()`) determines number of mesh points as `(max_size - 127) / bin_width`
   - What's unclear: The actual maximum size parameter for QUERUB from the RDS data — determines the n in the n×n kernel matrix. A 500mm max DBH at bin_width=1 would be ~373 mesh points, giving a 373×373 matrix with 139,000+ cells.
   - Recommendation: Check `max(params$sizeIngrowth['Lmax'])` for QUERUB in the benchmark script setup and log the mesh size in PROFILING.md.

2. **truncnorm in C++ migration**
   - What we know: `ingrowth_lk()` calls `truncnorm::dtruncnorm()` — this is a compiled R package function, not a base R function. Calling it from C++ requires `Rcpp::Function("dtruncnorm")` or reimplementing the density.
   - What's unclear: Whether migrating `ingrowth_lk` to C++ is worth the added complexity vs just vectorizing the R function.
   - Recommendation: Only consider C++ migration for `P_xEC` (pure arithmetic: dnorm + logistic) before `ingrowth_lk`. Profile results will determine if `ingrowth_lk` is even a bottleneck.

3. **Does profvis reliably capture C++ time?**
   - What we know: profvis uses Rprof() which samples R's call stack. C++ code called via Rcpp shows up as the R-level wrapper function (e.g., `getEigenValues`) — the internal C++ time is attributed to the wrapper.
   - What's unclear: Whether time spent inside `getEigenValues()` → Eigen::EigenSolver is accurately captured or appears as a single block.
   - Recommendation: If eigenvalue computation is suspected to dominate, benchmark `getEigenValues()` separately with `microbenchmark` on representative kernel sizes.

---

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | testthat (version via DESCRIPTION Suggests) |
| Config file | `tests/testthat.R` (exists) |
| Quick run command | `devtools::test(filter = "bug-fixes")` |
| Full suite command | `devtools::check()` |

### Phase Requirements -> Test Map

Phase 6 has no formal requirement IDs (phase requirement IDs: null). The verification criteria are behavioral:

| Behavior | Test Type | Automated Command | File Exists? |
|----------|-----------|-------------------|-------------|
| R CMD check passes after tidyverse revamp | integration | `devtools::check()` | Yes — existing suite |
| testthat suite passes after revamp | unit/integration | `devtools::test()` | Yes — existing suite |
| `lambda()` returns same value before and after optimization | regression | custom assertion in benchmark script | No — Wave 0 gap |
| `project()` returns same lambda trajectory before and after optimization | regression | custom assertion in benchmark script | No — Wave 0 gap |

### Sampling Rate
- **Per task commit:** `devtools::test()` — confirms existing tests pass
- **Per wave merge:** `devtools::check()` — full R CMD check, 0 errors, 0 warnings
- **Phase gate:** Full suite green before verification

### Wave 0 Gaps
- [ ] Regression assertions in `profile_benchmark.R` — verify `lambda(QUERUB)` value matches pre-optimization baseline before declaring any optimization correct. Use `testthat::expect_equal(lambda_new, lambda_baseline, tolerance = 1e-10)`.
- [ ] Baseline capture step: run `lambda()` and `project()` before any optimization and save numeric output to a `.rds` fixture file for comparison after each optimization step.

---

## Sources

### Primary (HIGH confidence)
- https://profvis.r-lib.org/articles/profvis.html — profvis() API, interval parameter, height parameter, basic usage pattern
- https://adv-r.hadley.nz/perf-improve.html — vectorization patterns; replacing per-element loops with vectorized C operations
- https://adv-r.hadley.nz/perf-measure.html — profvis + microbenchmark workflow; sampling profiler interpretation
- https://profvis.r-lib.org/reference/profvis.html — profvis() function reference; htmlwidgets::saveWidget() for HTML export

### Secondary (MEDIUM confidence)
- https://purrr.tidyverse.org/news/index.html — purrr 1.0.0 changelog; map_dfr() deprecation, as_vector() removal
- https://dplyr.tidyverse.org/articles/programming.html — .data$ pronoun for R CMD check compliance
- https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf — RcppEigen matrix operations; existing `src/eigen.cpp` confirmed as correct template

### Tertiary (LOW confidence)
- https://privefl.github.io/blog/Tip-Optimize-your-Rcpp-loops/ — Rcpp loop optimization patterns; flagged LOW because specific claim about 2-3x speedup is not independently verified for this codebase

---

## Metadata

**Confidence breakdown:**
- Standard stack (profvis, microbenchmark): HIGH — verified against official r-lib.org docs
- Architecture (sub-step sequence, outer() vectorization pattern): HIGH — verified against Advanced R official docs and direct code inspection
- Pitfalls (outer() argument order, Rcpp::compileAttributes()): HIGH — derived from direct code inspection + official Rcpp docs
- Pitfalls (caching correctness in multi-species case): MEDIUM — derived from code inspection of `project.R`, not confirmed with a live test
- C++ migration patterns: MEDIUM — template from existing `src/eigen.cpp` is confirmed; truncnorm replacement is inferred

**Research date:** 2026-03-11
**Valid until:** 2026-06-11 (stable ecosystem; profvis, purrr, dplyr APIs are stable)
