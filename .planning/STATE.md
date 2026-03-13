---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: executing
stopped_at: Completed 05-01-PLAN.md — compare_versions.R written and verified 15/15 PASS
last_updated: "2026-03-13T00:00:00.000Z"
last_activity: 2026-03-13 — Plan 05-01 executed and human-verified
progress:
  total_phases: 6
  completed_phases: 4
  total_plans: 8
  completed_plans: 8
  percent: 75
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Any researcher can run a full IPM for any supported species at any location by installing the package — no local data setup, no manual sourcing, no cluster required.
**Current focus:** Phase 2 - Package Skeleton and Bug Fixes

## Current Position

Phase: 5 of 6 (Compare Old vs New Package Results)
Plan: 1 of 1 in current phase — COMPLETE
Status: Phase 5 complete — compare_versions.R written, 15/15 rows PASS at 1e-10
Last activity: 2026-03-13 — Plan 05-01 executed and human-verified

Progress: [████████░░] 75%

## Performance Metrics

**Velocity:**
- Total plans completed: 1
- Average duration: ~60 min
- Total execution time: ~1 hour

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 1 - API Design | 1 | ~60 min | ~60 min |

**Recent Trend:**
- Last 5 plans: 01-01 (60 min)
- Trend: On track

*Updated after each plan completion*
| 2 - Package Skeleton and Bug Fixes | 2 | 8 min | 4 min |
| Phase 02-package-skeleton-and-bug-fixes P03 | 76 | 3 tasks | 21 files |
| Phase 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains P01 | 85 | 2 tasks | 8 files |
| Phase 06 P02 | 16 | 1 tasks | 5 files |
| Phase 06 P03 | 120 | 2 tasks | 6 files |
| Phase 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains P03 | 525548 | 3 tasks | 6 files |

## Accumulated Context

### Roadmap Evolution

- Phase 1 (API Design) inserted before all existing phases: existing phases renumbered 1→2, 2→3, 3→4
- Phase 5 added: Create script to compare results between old version of the package and this new versions (specially eigen computation). We can use simulation outputs from simulations folder and compare their results with the same output from our new R package
- Phase 6 added: Run tests with profvis to check for potential efficiency gains

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Setup]: Parquet for cloud parameter storage — supports partial reads via Arrow; columnar format efficient for species/draw subsets
- [Setup]: Tiered access: Parquet + ML decoder — full fidelity for known plots, extrapolation for new locations (ML decoder deferred to v3)
- [Setup]: Keep RcppEigen for eigenvalue solver — already implemented; C++ speedup critical for iterative projection
- [Setup]: Cloud hosting platform TBD — needs HTTP range request support (AWS S3 confirmed; Zenodo/HuggingFace need curl verification)
- [01-01]: Five constructors + two engines replaces run_ipm() — fulfills API-03 with composable design
- [01-01]: Single ipm_spModel class for both single and multi-species — avoids API branching
- [01-01]: env_condition() singular is canonical — env_conditions() plural is never used
- [01-01]: plot_size is a column in stand() data argument, not a separate argument
- [01-01]: store_every in control() enables memory management for long (>200 year) runs
- [01-01]: stand() size constraint is >= 127 mm (minimum DBH threshold for inventory data)
- [Phase 02-01]: Use Eigen::EigenSolver with computeEigenvectors=false; return .real() for VectorXd — IPM kernels are non-negative so dominant eigenvalue is real by Perron-Frobenius
- [Phase 02-01]: Initialize fct <- 1 before while-loop in init_pop() instead of exists() check — mathematically equivalent and correct for local scope
- [Phase 02-01]: Replace purrr::as_vector() with unlist(use.names=TRUE) — preserves named vector contracts for downstream parameter indexing
- [Phase 02-02]: Species IDs use short form (ABIBAL) in API — numeric-prefix form (18032ABIBAL) in raw data to be resolved in Phase 3
- [Phase 02-02]: Single-species competition: pass same Nvec as both Nvec_intra and Nvec_inter; multi-species: aggregate onto focal mesh via stats::approx()
- [Phase 02-02]: plot_random = c(0, 0, 0) in Phase 2 — plot-level random effects from Parquet in Phase 3
- [Phase 02-02]: Three-layer S3 pattern: new_<class>() (no validation) -> validate_<class>() (cli_abort) -> user-facing helper
- [Phase 02-03]: Use roxygen2::roxygenise(load_code='source') instead of devtools::document() — avoids Xcode license requirement for C++ compilation during documentation generation
- [Phase 02-03]: Add .Rbuildignore with ^simulations$ pattern — simulations/ is 61GB; without exclusion R CMD build takes 30+ minutes and produces 12-22GB tarballs
- [Phase 02-03]: Add .data$ pronoun to all bare dplyr column references in params.R to fix 'no visible binding' check warnings
- [Phase 06-01]: Force-add test fixture RDS files with git add -f to override *.rds gitignore rule; add exception pattern for future tracking
- [Phase 06-01]: species_model(s) is the correct call — plan template used wrong species= argument that does not exist in function signature
- [Phase 06-01]: Pre-existing test failures in test-bug-fixes.R (wrong file paths, pkgload::pkg_path() in R CMD check) deferred to deferred-items.md; out of scope for plan 06-01
- [Phase 06]: F matrix (outer + ingrowth_lk + truncnorm::dtruncnorm) is 2.5x more expensive than P matrix — profiling evidence changes Plan 03 priority order; Targets 1+2 must be co-applied
- [Phase 06]: truncnorm::dtruncnorm is 46% inclusive hotspot — replace with dnorm/pnorm direct formula in ingrowth_lk (math equivalent, base R vectorizes natively after outer() vectorization)
- [Phase 06]: Human approved Target 1 (vectorize outer) and Target 2 (replace truncnorm::dtruncnorm with dnorm/pnorm) for Plan 03 implementation; Targets 3 and 4 deferred
- [Phase 06]: outer() semantics: X_expanded=rep(X,times=n), Y_expanded=rep(Y,each=n); s1=rep(meshpts,times=n) for size_t1, s0=rep(meshpts,each=n) for size_t0 in vectorized mkKernel()
- [Phase 06]: Target 2 (truncnorm elimination) is the dominant gain source — F matrix 2.4x faster; Target 1 (outer vectorization) alone gives negligible speedup for P matrix
- [Phase 06]: Remove truncnorm from DESCRIPTION Imports — dependency eliminated now that ingrowth_lk uses dnorm/pnorm; closes PITFALLS.md Pitfall 11
- [Phase 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains]: Human approved Plan 03 optimizations: regression tests confirmed 1e-10 tolerance; 1.4x mkKernel() speedup and 2.4x F-matrix speedup accepted as complete
- [Phase 05-01]: final_output.RDS is a compact list of 216021 elements (not 216027) — 6 missing IDs are absent rather than NULL-slotted; direct indexing by array_id fails above index 32793; use fo_idx lookup from array_id column
- [Phase 05-01]: copy ref_output$par.BA_het for NULL-het rows to avoid NaN tolerance comparison (denominator is 0 when N_het_pertb = N_het = empty distribution)
- [Phase 05-01]: pars_to_list() native pipe incompatibility with .[[]] placeholder — fixed by intermediate variable assignment after group_split()
- [Phase 05-01]: All three key changes (RcppEigen solver, vectorized mkKernel, dnorm/pnorm ingrowth) confirmed numerically identical to cluster reference at 1e-10 — 15/15 stratified rows PASS

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 2 readiness]: Cloud hosting platform not yet selected — must verify HTTP 206 range request support via `curl -I --range 0-100 <url>` before writing fetch code. AWS S3 is confirmed; Zenodo and HuggingFace Hub need runtime verification.
- [Phase 2 readiness]: `httptest2` interceptability of Arrow's libcurl calls is unconfirmed — R-level HTTP mocking may not reach the C++ layer. Fallback is `skip_on_ci()` for network tests + RDS fixture testing.
- Xcode license not accepted — C++ compilation blocked until sudo xcodebuild -license accept is run in a terminal

## Session Continuity

Last session: 2026-03-13T00:44:39Z
Stopped at: Completed 05-01-PLAN.md — compare_versions.R written, verified 15/15 PASS, SUMMARY.md created
Resume file: None
