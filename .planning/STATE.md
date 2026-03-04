---
gsd_state_version: 1.0
milestone: v1.0
milestone_name: milestone
status: executing
stopped_at: "Phase 02 plan 02-02 complete. Eight API functions implemented: 5 constructors, 2 engines, 1 utility. Ready for plan 02-03 (DESCRIPTION/NAMESPACE)."
last_updated: "2026-03-04T15:25:11Z"
last_activity: 2026-03-04 — Plan 02-02 executed and verified
progress:
  total_phases: 4
  completed_phases: 1
  total_plans: 4
  completed_plans: 3
  percent: 62
---

# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Any researcher can run a full IPM for any supported species at any location by installing the package — no local data setup, no manual sourcing, no cluster required.
**Current focus:** Phase 2 - Package Skeleton and Bug Fixes

## Current Position

Phase: 2 of 4 (Package Skeleton and Bug Fixes)
Plan: 2 of 3 in current phase — COMPLETE
Status: Phase 2 in progress — plan 02-02 complete, ready for 02-03
Last activity: 2026-03-04 — Plan 02-02 executed and verified

Progress: [██████░░░░] 62%

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

## Accumulated Context

### Roadmap Evolution

- Phase 1 (API Design) inserted before all existing phases: existing phases renumbered 1→2, 2→3, 3→4

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

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 2 readiness]: Cloud hosting platform not yet selected — must verify HTTP 206 range request support via `curl -I --range 0-100 <url>` before writing fetch code. AWS S3 is confirmed; Zenodo and HuggingFace Hub need runtime verification.
- [Phase 2 readiness]: `httptest2` interceptability of Arrow's libcurl calls is unconfirmed — R-level HTTP mocking may not reach the C++ layer. Fallback is `skip_on_ci()` for network tests + RDS fixture testing.
- Xcode license not accepted — C++ compilation blocked until sudo xcodebuild -license accept is run in a terminal

## Session Continuity

Last session: 2026-03-04T15:25:11Z
Stopped at: Phase 02 plan 02-02 complete. Eight API functions implemented: 5 constructors, 2 engines, 1 utility. Ready for plan 02-03 (DESCRIPTION/NAMESPACE).
Resume file: None
