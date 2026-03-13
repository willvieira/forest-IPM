# Phase 4: Documentation and Tests - Context

**Gathered:** 2026-03-13
**Status:** Ready for planning

<domain>
## Phase Boundary

Deliver three things: (1) a testthat suite with 80%+ coverage of the main workflow, (2) roxygen2 docs with runnable examples for all exported functions, and (3) a complete book restructuring — rewrite `guide_IPM.qmd` to the new package API, add a new deep-dive vignette, and reorganize the book sections in `_quarto.yml`. A new `plot.ipm_projection()` S3 method is also in scope as it is required by the vignette.

Phase 3 (Remote Data Layer / cloud Parquet) is NOT a dependency — `parameters()` works with `plot_random = c(0, 0, 0)` for all examples and tests.

</domain>

<decisions>
## Implementation Decisions

### Tests

- Remove all existing test files **except** `test-constructors.R` — the others are superseded or have pre-existing failures that are out of scope
- New test suite covers the full main workflow: `stand()` → `species_model()` → `parameters()` → `lambda()` → `project()`
- Target: **80%+ line coverage** on main workflow functions (`stand.R`, `species_model.R`, `parameters.R`, `kernel.R`, `lambda.R`, `project.R`)
- Use `plot_random = c(0, 0, 0)` throughout — no cloud/Parquet dependency needed
- Use ABIBAL as the test species
- DATA requirements (cloud fetching, cache layer) are out of scope — those belong in Phase 3's own tests
- Pre-existing `test-bug-fixes.R` failures remain deferred

### Documentation (roxygen2)

- All exported functions (~25) need `@examples` sections added
- Examples use `parameters()` with `plot_random = c(0, 0, 0)` and ABIBAL — no network required
- Examples must pass `devtools::run_examples()` without errors

### plot.ipm_projection() — new export

- Add `plot.ipm_projection(x, type = NULL, ...)` to the R package
- `type` argument selects figure: `"lambda"` | `"size_dist"` | `"lambda_vs_n"`
- Default (type = NULL): render all three figures in sequence
- Fully documented with `@param`, `@return`, `@examples`
- Exported from NAMESPACE

### guide_IPM.qmd — full rewrite

- Replace old procedural API (source from GitHub + `getPars_sp()` + `mkKernel()` + manual `eigen()`) with the new `library(forestIPM)` package API
- **Two-part structure:**
  1. **Quick start** — simple single-species workflow: `stand()` → `species_model()` → `parameters()` → `env_condition()` → `control()` → `lambda()` → `project()`. Species: ABIBAL, no competition (`plot_random = c(0,0,0)`), static climate
  2. **Deep dive** — multi-species example, time-varying climate (environmental change), and a roadmap section for future developments (ML decoder, coexistence IPM, CRAN)
- Visualizations rendered via `plot(proj, type = ...)` — no manual ggplot boilerplate in vignette
- Keep all three figures: lambda over time, animated size distribution (plotly), lambda vs N correlation
- Reference to new book sections where appropriate

### Book restructuring — full

- Current structure has 2 sections: "Demographic models evaluation" + "Using demographic models to infer population level performance"
- **New 3-section structure:**
  1. Section 1: **Demographic models evaluation** — existing chapters (unchanged)
  2. Section 2: **From demographic rates to population dynamics** — `guide_IPM.qmd` (rewritten) + new deep-dive .qmd
  3. Section 3: **Using the IPM** — existing analysis chapters reorganized here (conditional lambda, coexistence, extinction risk, suitable probability, sensitivity analysis)
- Update `_quarto.yml` in `/Users/wvieira/GitHub/book_forest-demography-IPM` to reflect new section structure
- Chapter filenames stay the same; only `_quarto.yml` navigation and section labels change

### Claude's Discretion

- Deep-dive .qmd filename and internal section titles
- Which specific edge cases to cover in tests (beyond the happy path workflow)
- How to handle the plotly animation in a static Quarto build (whether to wrap in `eval: false` or use a precomputed RDS as the old guide did)
- Exact book section names (user gave examples: "Using the IPM to infer population performance", "evaluating IPM", "Using the IPM" — any of these works)

</decisions>

<specifics>
## Specific Ideas

- `plot.ipm_projection()` should mirror `plot.lm()` idiom: `plot(proj)` renders all, `plot(proj, type = "lambda")` renders one
- Quick start example: ABIBAL, single species, no competition, static climate (MAT = 8, MAP = 1200 or similar temperate values)
- Roadmap section in deep-dive vignette: ML decoder for unobserved locations, multi-species coexistence IPM, CRAN submission
- Book section name examples (user's words): "From demographic rates to population dynamics", "Using the IPM to infer population performance"

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets

- `tests/testthat/test-constructors.R`: Keep as-is — covers constructor validation
- `tests/testthat/fixtures/`: `lambda_baseline.rds`, `project_baseline.rds` — regression baselines may be reusable
- All `R/*.R` files already have `@param` and `@return` — only `@examples` sections are missing
- `/Users/wvieira/GitHub/book_forest-demography-IPM/guide_IPM.qmd`: Source for rewrite; plots and narrative structure can be adapted

### Established Patterns

- `plot_random = c(0, 0, 0)` is the zero-effect baseline already used in Phase 2 and compare_versions.R
- ABIBAL is the reference species used throughout the book and existing simulation scripts
- Three-layer S3 pattern: `new_<class>()` → `validate_<class>()` → user-facing helper — `plot.ipm_projection()` follows this
- Old guide uses a precomputed RDS (`data/out_IPM_guide.RDS`) to avoid slow plotly animation at render time — same pattern appropriate for new guide

### Integration Points

- `plot.ipm_projection()` integrates with the `ipm_projection` S3 class returned by `project()`
- Book repo at `/Users/wvieira/GitHub/book_forest-demography-IPM` — separate git repo, changes committed there independently
- `_quarto.yml` in book repo controls section structure and chapter order

</code_context>

<deferred>
## Deferred Ideas

- None — discussion stayed within phase scope

</deferred>

---

*Phase: 04-documentation-and-tests*
*Context gathered: 2026-03-13*
