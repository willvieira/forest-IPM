# Roadmap: forestIPM

## Overview

Convert a working simulation codebase into a distributable R package that any researcher can install from GitHub and use to run Bayesian hierarchical IPMs for 31 tree species — no local 15 GB dataset required. Three phases: first make the package installable and bug-free, then wire in cloud parameter access via Arrow + Parquet, then finish with documentation and tests that allow CI to pass without network access.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: API Design** - Define the public API surface of the package before implementation
- [x] **Phase 2: Package Skeleton and Bug Fixes** - Installable R package with zero `R CMD check` errors and all known code defects resolved (completed 2026-03-04)
- [ ] **Phase 3: Remote Data Layer** - Cloud Parquet parameter access via Arrow with local cache — the core value proposition
- [ ] **Phase 4: Documentation and Tests** - roxygen2 docs on all exported functions and a testthat suite that passes on CI without network access

## Phase Details

### Phase 1: API Design
**Goal**: Define the public API surface of the package — exported function signatures, argument names, return types, and user-facing conventions — before any implementation begins
**Depends on**: Nothing (first phase)
**Requirements**: API-01, API-02, API-03, API-04, API-05
**Success Criteria** (what must be TRUE):
  1. All exported functions listed with finalized signatures and documented intent
  2. Naming conventions agreed upon and applied consistently
  3. API design document produced that guides Phase 2 implementation
**Plans**: 1 plan

Plans:
- [x] 01-01-PLAN.md — Write API_DESIGN.md: complete interface contract (all 8 exported functions, 6 S3 classes, naming conventions, validation rules)

### Phase 2: Package Skeleton and Bug Fixes
**Goal**: Researchers can install the package from GitHub and call core functions without errors — `devtools::install_github("wvieira/forest-IPM")` succeeds, `devtools::load_all()` works, and `R CMD check` produces zero errors and zero warnings
**Depends on**: Phase 1
**Requirements**: PKG-01, PKG-02, PKG-03, PKG-04, PKG-05, BUG-01, BUG-02, BUG-03
**Success Criteria** (what must be TRUE):
  1. `devtools::install_github("wvieira/forest-IPM")` completes without errors on a fresh R environment
  2. `R CMD check` produces zero errors and zero warnings (NOTEs acceptable temporarily)
  3. `devtools::load_all()` compiles the Rcpp C++ eigenvalue solver without errors
  4. `mkKernel()` and `getEigenValues()` return correct results — eigenvalue computed with general (non-symmetric) solver matches base R `eigen()` output on a test kernel
  5. Calling `init_pop()` with a function argument does not depend on or inspect any global variable named `fct`
**Plans**: 3 plans

Plans:
- [ ] 02-01-PLAN.md — Fix three code bugs: replace SelfAdjointEigenSolver with EigenSolver (BUG-01), remove exists('fct') global check (BUG-02), replace as_vector() with unlist() (BUG-03)
- [ ] 02-02-PLAN.md — Implement the Phase 1 API specification: 5 constructors (stand, species_model, parameters, env_condition, control), 2 engines (lambda, project), 1 utility (supported_species), all S3 classes and methods
- [ ] 02-03-PLAN.md — Create package infrastructure: DESCRIPTION, NAMESPACE via roxygen2/devtools::document(), R/globals.R with all @importFrom declarations; iterate R CMD check to zero errors and warnings

### Phase 3: Remote Data Layer
**Goal**: A researcher can fetch parameters for any of the 31 supported species from the cloud with a single function call — no local RDS files, no manual downloads, with results cached locally for subsequent calls
**Depends on**: Phase 2
**Requirements**: DATA-01, DATA-02, DATA-03, DATA-04, DATA-05, DATA-06
**Success Criteria** (what must be TRUE):
  1. `get_params("ACRU", model = "intcpt_plot_comp_clim", method = "mean")` returns a parameter object without downloading the full 15 GB dataset (verified by network traffic or Arrow predicate pushdown logs)
  2. A second call to `get_params()` for the same species reads from local disk cache — no HTTP request made
  3. The cache directory survives package reinstallation (located under `tools::R_user_dir("forestIPM", "cache")`)
  4. Examples and tests that use `run_ipm()` work without any internet connection using the bundled `sysdata.rda` mini-dataset
  5. Attempting to fetch parameters when the cloud host is unreachable produces an informative error message, not a cryptic timeout
**Plans**: TBD

### Phase 4: Documentation and Tests
**Goal**: Every exported function is documented, examples run without network access, and the testthat suite passes on CI without requiring the 15 GB parameter dataset
**Depends on**: Phase 3
**Requirements**: DOC-01, DOC-02
**Success Criteria** (what must be TRUE):
  1. `?get_params` and `?mkKernel` (and all other exported functions) open a help page with `@param`, `@return`, and `@examples` sections
  2. `devtools::run_examples()` completes without errors or warnings — all examples use the bundled mini-dataset, not live network calls
  3. `devtools::test()` passes with zero failures on a machine with no internet access
**Plans:** 3 plans

Plans:
- [ ] 04-01-PLAN.md — Redesign plot.ipm_projection() (three-panel type dispatch) and establish test infrastructure (delete obsolete tests, write test-workflow.R)
- [ ] 04-02-PLAN.md — Add @examples to all 11 exported functions and verify offline coverage
- [ ] 04-03-PLAN.md — Rewrite guide_IPM.qmd, create deep_dive_IPM.qmd, restructure _quarto.yml with 3-section layout

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3 → 4

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. API Design | 1/1 | Complete | 2026-03-03 |
| 2. Package Skeleton and Bug Fixes | 3/3 | Complete | 2026-03-04 |
| 3. Remote Data Layer | 0/TBD | Not started | - |
| 4. Documentation and Tests | 0/TBD | Not started | - |
| 5. Compare Old vs New Package Results | 1/1 | Complete | 2026-03-13 |

### Phase 5: Create script to compare results between old version of the package and this new versions (specially eigen computation). We can use simulation outputs from simulations folder and compare their results with the same output from our new R package

**Goal:** Verify that all Phase 02-01 (EigenSolver fix), Phase 06-01 (vectorized kernel), and Phase 06-03 (truncnorm elimination) changes produce lambda values numerically identical (within 1e-10) to the original cluster simulation outputs in final_output.RDS
**Requirements**: COMPARE-01, COMPARE-02, COMPARE-03
**Depends on:** Phase 4
**Plans:** 1 plan

Plans:
- [x] 05-01-PLAN.md — Write compare_versions.R: standalone script that loads new package, reproduces ipm_i.R parameter draw sequence, and diffs all 5 lambda metrics against final_output.RDS for ~15 stratified rows

### Phase 6: Run tests with profvis to check for potential efficiency gains

**Goal:** Profile the core IPM computation pipeline with profvis, revamp R/ code to consistent tidyverse style, and apply targeted optimizations based on profiling evidence
**Requirements**: TBD
**Depends on:** Phase 5
**Plans:** 3/3 plans complete

Plans:
- [ ] 06-01-PLAN.md — Capture pre-optimization baselines and apply tidyverse-style revamp to all R/ source files
- [ ] 06-02-PLAN.md — Write and run profvis benchmarking script, produce PROFILING.md with hotspot analysis and flamegraph
- [ ] 06-03-PLAN.md — Apply approved optimizations to kernel.R/project.R, run before/after benchmark, update PROFILING.md
