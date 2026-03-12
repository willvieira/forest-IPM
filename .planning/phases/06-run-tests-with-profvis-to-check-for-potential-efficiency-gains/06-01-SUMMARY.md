---
phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains
plan: "01"
subsystem: testing
tags: [purrr, tidyverse, regression-testing, ipm, r-package]

# Dependency graph
requires:
  - phase: 02-package-skeleton-and-bug-fixes
    provides: Working R package with lambda(), project(), stand(), species_model(), parameters(), env_condition() public API
provides:
  - Pre-optimization regression fixtures for QUERUB lambda and project() lambda trajectory (1e-10 tolerance)
  - Tidyverse-style revamp of all R/ source files (native pipe, purrr typed maps, .data$ pronouns)
affects: [07-optimization, any future profvis profiling phase]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Typed purrr maps (map_dbl/map_lgl/map_chr) instead of sapply/vapply throughout R/"
    - "Native pipe |> replaces magrittr %>% everywhere"
    - "Regression fixtures at tests/testthat/fixtures/*.rds as pre-optimization anchors"
    - ".gitignore !exception for test fixtures to override *.rds blanket exclusion"

key-files:
  created:
    - tests/testthat/test-regression-baselines.R
    - tests/testthat/fixtures/lambda_baseline.rds
    - tests/testthat/fixtures/project_baseline.rds
  modified:
    - R/BasalArea_competition.R
    - R/globals.R
    - R/lambda.R
    - R/project.R
    - R/species_model.R
    - NAMESPACE
    - .gitignore

key-decisions:
  - "Force-add test fixtures with git add -f to override *.rds gitignore rule; add !exception pattern for future tracking"
  - "species_model(s) without species= argument (function signature is species_model(x, on_missing)); plan template had wrong call"
  - "Pre-existing test failures in test-bug-fixes.R and test-package-structure.R are out of scope and deferred to deferred-items.md"

patterns-established:
  - "All new R code uses purrr typed maps (map_dbl/map_lgl/map_chr) — no sapply/vapply"
  - "All future code uses native |> pipe — no magrittr %>%"
  - "Regression fixtures at tests/testthat/fixtures/ lock pre-optimization numeric output"

requirements-completed: []

# Metrics
duration: 85min
completed: 2026-03-12
---

# Phase 06 Plan 01: Pre-optimization Baselines and Tidyverse Revamp Summary

**QUERUB regression fixtures at tolerance 1e-10 saved to tests/testthat/fixtures/ and all R/ source files converted from sapply/vapply/magrittr to typed purrr maps and native pipe**

## Performance

- **Duration:** ~85 min
- **Started:** 2026-03-12T00:07:57Z
- **Completed:** 2026-03-12T01:33:00Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments

- Captured QUERUB lambda (1.201558) and 10-year project() lambda trajectory as RDS fixtures in tests/testthat/fixtures/
- Wrote two regression test assertions (tolerance 1e-10) that guard any subsequent optimization from numeric drift
- Replaced all sapply/vapply calls in R/ with typed purrr equivalents (map_dbl, map_lgl, map_chr) across 4 files
- Removed @importFrom magrittr %>% from globals.R; all pipe usage already used native |>
- Regenerated NAMESPACE with correct importFrom declarations for new purrr symbols

## Task Commits

Each task was committed atomically:

1. **Task 1: Capture pre-optimization baselines and write regression test** - `c81658d` (test)
2. **Task 2: Tidyverse-style revamp of all R/ source files** - `758b05a` (feat)

_Note: Task 1 used TDD flow: RED (tests skip with missing fixtures) -> fixtures generated -> GREEN (both tests pass)_

## Files Created/Modified

- `tests/testthat/test-regression-baselines.R` - Two regression tests asserting lambda() and project() output matches saved baselines at 1e-10 tolerance
- `tests/testthat/fixtures/lambda_baseline.rds` - QUERUB lambda=1.201558 (ipm_lambda class, named vector)
- `tests/testthat/fixtures/project_baseline.rds` - QUERUB 10-year project() lambda trajectory (10 numeric values)
- `R/BasalArea_competition.R` - sapply() -> purrr::map_dbl() in size_to_BAcomp()
- `R/lambda.R` - vapply(focal_species, ..., numeric(1)) -> purrr::map_dbl() with setNames
- `R/project.R` - 3x vapply(..., logical(1)) -> purrr::map_lgl()
- `R/species_model.R` - vapply(bad_ids, ..., character(1)) -> purrr::map_chr(); vapply(params, is.null, logical(1)) -> purrr::map_lgl()
- `R/globals.R` - Added map_dbl/map_lgl/map_chr to @importFrom purrr; removed @importFrom magrittr %>%
- `NAMESPACE` - Regenerated via roxygen2::roxygenise(load_code='source')
- `.gitignore` - Added !tests/testthat/fixtures/*.rds exception

## Decisions Made

- Force-added test fixture RDS files with `git add -f` since *.rds is gitignored package-wide; added `!tests/testthat/fixtures/*.rds` exception to .gitignore for future tracking
- The plan template used `species_model(s, species = "QUERUB")` which is not the actual function signature; corrected to `species_model(s)` (the function only accepts `x` and `on_missing` arguments)
- Pre-existing test failures in `test-bug-fixes.R` (BUG-02/BUG-03 referencing wrong file `R/params.R` instead of `R/parameters.R`, and `pkgload::pkg_path()` failing in R CMD check context) were documented in `deferred-items.md` and excluded from this plan's scope

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] .gitignore blocked committing test fixture RDS files**
- **Found during:** Task 1 (commit step)
- **Issue:** `*.rds` blanket rule in .gitignore prevented `git add tests/testthat/fixtures/*.rds`
- **Fix:** Added `!tests/testthat/fixtures/*.rds` exception line to .gitignore; force-added the files with `git add -f` to override the ignore rule for already-tracked exceptions
- **Files modified:** .gitignore
- **Verification:** Both fixture files staged and committed successfully
- **Committed in:** c81658d (Task 1 commit)

**2. [Rule 1 - Bug] Plan template used wrong species_model() argument**
- **Found during:** Task 1 (test writing)
- **Issue:** Plan template had `species_model(s, species = "QUERUB")` but the function signature is `species_model(x, on_missing = "error")` — no `species` argument exists
- **Fix:** Corrected to `species_model(s)` in both the fixture generation script and the regression test file
- **Files modified:** tests/testthat/test-regression-baselines.R
- **Verification:** Tests pass; lambda output matches baseline exactly (diff = 0)
- **Committed in:** c81658d (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered

- **Pre-existing test failures:** `test-bug-fixes.R` has 3 pre-existing failures (BUG-02 test expects `fct <- 1` in kernel.R which isn't there, BUG-03 tests reference non-existent `R/params.R`). These are out of scope and documented in `deferred-items.md`.
- **R CMD check errors:** `pkgload::pkg_path()` called at top-level outside `test_that()` in test-bug-fixes.R and test-package-structure.R causes errors in the R CMD check environment. Pre-existing issue, documented in `deferred-items.md`.
- **`devtools::check()` produces 1 error, 3 warnings, 5 notes** — all pre-existing; regression and non-bug-fixes tests pass cleanly.

## Next Phase Readiness

- Regression baselines established — any optimization can be validated against them
- All R/ source files use idiomatic tidyverse style (native pipe, typed purrr maps, .data$ pronouns)
- Pre-existing test infrastructure issues documented in deferred-items.md for Phase 06-02 or maintenance work

## Self-Check: PASSED

- FOUND: tests/testthat/test-regression-baselines.R
- FOUND: tests/testthat/fixtures/lambda_baseline.rds
- FOUND: tests/testthat/fixtures/project_baseline.rds
- FOUND: 06-01-SUMMARY.md
- FOUND commit c81658d (Task 1)
- FOUND commit 758b05a (Task 2)
- No %>% in R/ code
- No sapply/vapply in R/ code

---
*Phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains*
*Completed: 2026-03-12*
