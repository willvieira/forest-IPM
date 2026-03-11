# Phase 6: Run tests with profvis to check for potential efficiency gains - Context

**Gathered:** 2026-03-11
**Status:** Ready for planning

<domain>
## Phase Boundary

Profile the core IPM computation pipeline with profvis, revamp R/ code toward consistent tidyverse style, and apply targeted optimizations based on profiling evidence and feasibility. Three sequential sub-steps: (1) tidyverse revamp, (2) profvis profiling, (3) code optimization.

</domain>

<decisions>
## Implementation Decisions

### Sub-step 1: Tidyverse-style revamp

- **Goal:** Consistency and readability — make all R/ code consistently idiomatic (pipes, purrr::map, dplyr verbs). Performance improvement is a secondary benefit, not the primary goal.
- **Scope:** All R/ files — params.R, project.R, lambda.R, kernel.R, stand.R, BasalArea_competition.R, vital_rates.R, and all other files in R/. Full pass for consistency.
- **Boundary:** R/ source files only — do not touch test files or roxygen2 @examples (those belong to Phase 4).
- **Skill to use:** `r-skills:tidyverse-patterns`
- **Verification:** `devtools::check()` (R CMD check) + existing testthat suite must pass after revamp.

### Sub-step 2: profvis profiling

- **Benchmark inputs:**
  - Single-species run: 1 species (QUERUB — large max size stresses big kernel), 100 timesteps via `project()`
  - Multi-species run: 2–3 species including QUERUB, 100 timesteps — exposes competition update costs
  - Rationale: single vs multi-species expose different code paths
- **Skill to use:** `r-skills:r-performance`
- **Output deliverables:**
  - A written markdown summary (`PROFILING.md`) in the phase directory documenting: top hotspots, call stack breakdown, time percentages, and recommended optimization targets for sub-step 3
  - An HTML profvis flamegraph saved to the phase directory as an artifact
- **Script location:** `.planning/phases/06-.../` (not shipped with the package)

### Sub-step 3: Targeted optimizations

- **Decision framework:** Apply if measurable speedup (>10–15% on benchmark) AND implementation is low-risk (no change to math/statistical outputs). Document anything evaluated but skipped and why.
- **C++ bias:** C++ migrations (via Rcpp) should be weighted favorably — the package already depends on RcppEigen, so adding more C++ functions is low dependency overhead.
- **Strategies in scope** (all conditional on profvis evidence):
  - Vectorize `outer()` calls in `mkKernel()` — replace `outer(meshpts, meshpts, FUN)` with pre-allocated matrix + vectorized operations to avoid per-cell R function call overhead
  - Cache competition metrics — `BA_comp_intra/inter` are recomputed inside `mkKernel()` every call; cache when Nvec hasn't changed between timesteps
  - C++ migration of inner loops — move `P_xEC` or `ingrowth_lk` to C++ via Rcpp if they appear as hotspots in profvis flamegraph
- **Final check:** Run benchmark comparison (before/after) to quantify gains from each applied optimization.

### Claude's Discretion

- Exact profvis flamegraph presentation format (HTML styling, chart options)
- Which specific dplyr/purrr patterns to use in each file (follow tidyverse-patterns skill guidance)
- How to structure the PROFILING.md report
- Whether to combine optimizations into one commit or separate commits per optimization

</decisions>

<code_context>
## Existing Code Insights

### Reusable Assets
- `src/eigen.cpp`: Existing C++ file with RcppEigen — template for adding new Rcpp functions if inner loops get migrated
- `R/RcppExports.R`: Auto-generated Rcpp glue — any new C++ functions will appear here after `Rcpp::compileAttributes()`

### Established Patterns
- `mkKernel()`: Two `outer(meshpts, meshpts, FUN)` calls (P matrix and F matrix) — likely the dominant compute cost; each builds an n×n matrix calling R functions per cell
- `project()` loop: calls `mkKernel()` + matrix-vector multiply (`K %*% Nvec`) per timestep per species — repeated kernel reconstruction is the main cost driver
- `BasalArea_competition.R`: `size_to_BAcomp()` and `size_to_BAplot()` called inside `mkKernel()` every iteration — potential cache target
- `vital_rates.R`: `vonBertalanffy_f()`, `survival_f()`, `ingrowth_f()` — called inside `P_xEC` and `ingrowth_lk` which are called once per matrix cell by `outer()`

### Integration Points
- Any C++ additions go in `src/` and must be registered via `Rcpp::compileAttributes()`
- Style revamp must not change function signatures or return types (public API is locked from Phase 1)
- `NAMESPACE` and `globals.R` may need updating if new purrr/dplyr functions are introduced

</code_context>

<specifics>
## Specific Ideas

- Use QUERUB (Quercus rubra) as the benchmark species — large maximum DBH means larger mesh (more size classes), which stresses the `outer()` kernel build more than small-statured species
- The profvis script should benchmark both `lambda()` and `project()` separately, since they exercise slightly different code paths (lambda does one kernel build; project does N builds in a loop)

</specifics>

<deferred>
## Deferred Ideas

- None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains*
*Context gathered: 2026-03-11*
