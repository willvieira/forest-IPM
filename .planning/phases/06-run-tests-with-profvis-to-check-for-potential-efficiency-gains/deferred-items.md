# Deferred Items — Phase 06

## Pre-existing test failures discovered during 06-01 execution

These failures existed before any changes in plan 06-01 and are out of scope for this plan.

### 1. BUG-02 test: `fct <- 1` expected in R/kernel.R

- **File:** tests/testthat/test-bug-fixes.R:35
- **Issue:** Test asserts `any(grepl("fct <- 1", src_lines))` in R/kernel.R but this line does not exist.
- **Root cause:** The `fct <- 1` fix from Phase 02-01 may have been applied to a different location (e.g., a separate `init_pop.R` that was later merged or renamed), or the test was written for a file structure that no longer matches.
- **Action needed:** Investigate where `fct <- 1` was placed during Phase 02-01 and either update the test to point to the correct file or add the initializer to kernel.R.

### 2. BUG-03 tests: `R/params.R` file not found

- **Files:** tests/testthat/test-bug-fixes.R:40, :45
- **Issue:** Tests read `R/params.R` but the file is named `R/parameters.R`.
- **Root cause:** The file was renamed during Phase 02-02 or 02-03 but the tests were not updated.
- **Action needed:** Update test-bug-fixes.R to reference `R/parameters.R` instead of `R/params.R`.

### 3. `pkgload::pkg_path()` fails in R CMD check context

- **Files:** test-bug-fixes.R:3, test-package-structure.R:3
- **Issue:** `pkg <- pkgload::pkg_path()` at the top-level outside `test_that()` blocks fails when tests run from the installed package temp directory during `devtools::check()`.
- **Root cause:** `pkgload::pkg_path()` requires a DESCRIPTION file in the working directory, which is not present in the check temp environment.
- **Action needed:** Wrap `pkg_path()` calls inside `test_that()` blocks or use `skip_if_not(...)` guards.
