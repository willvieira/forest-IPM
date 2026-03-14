---
phase: 07-add-github-actions-ci-with-code-coverage-and-pkgdown-site-with-readme-pointing-to-book-vignettes
plan: "01"
subsystem: ci
tags: [github-actions, ci, code-coverage, codecov, r-cmd-check]
dependency_graph:
  requires: []
  provides: [CI_R_CMD_CHECK, CI_COVERAGE]
  affects: []
tech_stack:
  added:
    - r-lib/actions/check-r-package@v2
    - r-lib/actions/setup-r-dependencies@v2
    - codecov/codecov-action@v5
    - covr (R package for coverage measurement)
  patterns:
    - GitHub Actions YAML workflow
    - Cobertura XML format for codecov upload
key_files:
  created:
    - .github/workflows/R-CMD-check.yaml
    - .github/workflows/test-coverage.yaml
  modified: []
decisions:
  - Linux-only (ubuntu-latest) with R release — no matrix, no Windows, no macOS rows to keep CI simple and fast
  - Two separate workflow files — coverage failures do not block check result and vice versa
  - codecov/codecov-action@v5 with fail_ci_if_error conditional on push (not PR) to avoid noise
  - CODECOV_TOKEN is optional for public repos; note added for user to add it after activating on codecov.io
  - C++ coverage (src/eigen.cpp) excluded from covr by design; 70% target applies to R code only
metrics:
  duration_minutes: 10
  completed_date: "2026-03-14"
  tasks_completed: 2
  files_changed: 2
requirements_fulfilled:
  - TST-01
  - TST-02
---

# Phase 07 Plan 01: GitHub Actions CI Workflows Summary

**One-liner:** R CMD check and covr/codecov CI via r-lib/actions on ubuntu-latest with separate workflow files for independent failure modes.

## What Was Built

Two GitHub Actions workflow files that provide automated quality gates on every push and pull request to main/master:

1. `.github/workflows/R-CMD-check.yaml` — runs R CMD check on ubuntu-latest with R release using r-lib/actions/check-r-package@v2; setup-r-dependencies automatically resolves system dependencies for Rcpp/RcppEigen via pak.

2. `.github/workflows/test-coverage.yaml` — measures R code coverage via `covr::package_coverage()`, converts to Cobertura XML via `covr::to_cobertura()`, and uploads to codecov.io using codecov/codecov-action@v5; covr and xml2 are installed as extra-packages alongside the standard dependency resolution.

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Create R-CMD-check.yaml workflow | f8b734b | .github/workflows/R-CMD-check.yaml |
| 2 | Create test-coverage.yaml workflow | 56886b3 | .github/workflows/test-coverage.yaml |

## Decisions Made

- **Linux-only CI:** ubuntu-latest only, R release only — no matrix expansion. Keeps CI simple, fast, and directly relevant to the primary deployment target.
- **Separate workflow files:** R-CMD-check and test-coverage are independent. A coverage upload failure (e.g., missing CODECOV_TOKEN) does not block the check result badge.
- **pak-based system deps:** `setup-r-dependencies` detects and installs build-essential for Rcpp/RcppEigen automatically. Manual `apt-get install` only needed if `g++: not found` appears in first CI run.
- **codecov-action@v5 conditional failure:** `fail_ci_if_error` is true on push (non-PR) to catch persistent upload issues, but false on PRs to avoid blocking contributor PRs from missing CODECOV_TOKEN.
- **C++ coverage excluded:** `src/eigen.cpp` is not measured by covr — this is expected; covr instruments R code only. The 70% coverage target applies to R code.

## Deviations from Plan

None - plan executed exactly as written.

## Pre-existing Issues (Out of Scope)

The test suite has one pre-existing failure in `test-regression-baselines.R:23` — `lambda()` output attribute mismatch for QUERUB species. This failure pre-dates this plan and is unrelated to the workflow files created here. Logged to deferred-items for future investigation.

## Verification

- `.github/workflows/R-CMD-check.yaml` — EXISTS with ubuntu-latest runner, R release, check-r-package@v2
- `.github/workflows/test-coverage.yaml` — EXISTS with covr::package_coverage(), covr::to_cobertura(), codecov-action@v5
- No macOS or Windows runners in either workflow
- Test suite: 70 PASS, 1 FAIL (pre-existing regression baseline issue, unrelated to this plan)

## Self-Check: PASSED

- /Users/wvieira/GitHub/forest-IPM/.github/workflows/R-CMD-check.yaml — FOUND
- /Users/wvieira/GitHub/forest-IPM/.github/workflows/test-coverage.yaml — FOUND
- Commit f8b734b — FOUND (git log confirms)
- Commit 56886b3 — FOUND (git log confirms)
