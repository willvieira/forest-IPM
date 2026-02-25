# Roadmap: forestIPM

## Overview

Convert a working simulation codebase into a distributable R package that any researcher can install from GitHub and use to run Bayesian hierarchical IPMs for 31 tree species — no local 15 GB dataset required. Three phases: first make the package installable and bug-free, then wire in cloud parameter access via Arrow + Parquet, then finish with documentation and tests that allow CI to pass without network access.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [ ] **Phase 1: Package Skeleton and Bug Fixes** - Installable R package with zero `R CMD check` errors and all known code defects resolved
- [ ] **Phase 2: Remote Data Layer** - Cloud Parquet parameter access via Arrow with local cache — the core value proposition
- [ ] **Phase 3: Documentation and Tests** - roxygen2 docs on all exported functions and a testthat suite that passes on CI without network access

## Phase Details

### Phase 1: Package Skeleton and Bug Fixes
**Goal**: Researchers can install the package from GitHub and call core functions without errors — `devtools::install_github("wvieira/forest-IPM")` succeeds, `devtools::load_all()` works, and `R CMD check` produces zero errors and zero warnings
**Depends on**: Nothing (first phase)
**Requirements**: PKG-01, PKG-02, PKG-03, PKG-04, PKG-05, BUG-01, BUG-02, BUG-03
**Success Criteria** (what must be TRUE):
  1. `devtools::install_github("wvieira/forest-IPM")` completes without errors on a fresh R environment
  2. `R CMD check` produces zero errors and zero warnings (NOTEs acceptable temporarily)
  3. `devtools::load_all()` compiles the Rcpp C++ eigenvalue solver without errors
  4. `mkKernel()` and `getEigenValues()` return correct results — eigenvalue computed with general (non-symmetric) solver matches base R `eigen()` output on a test kernel
  5. Calling `init_pop()` with a function argument does not depend on or inspect any global variable named `fct`
**Plans**: TBD

### Phase 2: Remote Data Layer
**Goal**: A researcher can fetch parameters for any of the 31 supported species from the cloud with a single function call — no local RDS files, no manual downloads, with results cached locally for subsequent calls
**Depends on**: Phase 1
**Requirements**: DATA-01, DATA-02, DATA-03, DATA-04, DATA-05, DATA-06
**Success Criteria** (what must be TRUE):
  1. `get_params("ACRU", model = "intcpt_plot_comp_clim", method = "mean")` returns a parameter object without downloading the full 15 GB dataset (verified by network traffic or Arrow predicate pushdown logs)
  2. A second call to `get_params()` for the same species reads from local disk cache — no HTTP request made
  3. The cache directory survives package reinstallation (located under `tools::R_user_dir("forestIPM", "cache")`)
  4. Examples and tests that use `run_ipm()` work without any internet connection using the bundled `sysdata.rda` mini-dataset
  5. Attempting to fetch parameters when the cloud host is unreachable produces an informative error message, not a cryptic timeout
**Plans**: TBD

### Phase 3: Documentation and Tests
**Goal**: Every exported function is documented, examples run without network access, and the testthat suite passes on CI without requiring the 15 GB parameter dataset
**Depends on**: Phase 2
**Requirements**: DOC-01, DOC-02
**Success Criteria** (what must be TRUE):
  1. `?get_params` and `?mkKernel` (and all other exported functions) open a help page with `@param`, `@return`, and `@examples` sections
  2. `devtools::run_examples()` completes without errors or warnings — all examples use the bundled mini-dataset, not live network calls
  3. `devtools::test()` passes with zero failures on a machine with no internet access
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Package Skeleton and Bug Fixes | 0/TBD | Not started | - |
| 2. Remote Data Layer | 0/TBD | Not started | - |
| 3. Documentation and Tests | 0/TBD | Not started | - |
