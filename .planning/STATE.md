# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Any researcher can run a full IPM for any supported species at any location by installing the package — no local data setup, no manual sourcing, no cluster required.
**Current focus:** Phase 1 - API Design

## Current Position

Phase: 1 of 4 (API Design)
Plan: 1 of 1 in current phase — COMPLETE
Status: Phase 1 complete — ready for Phase 2
Last activity: 2026-03-03 — Plan 01-01 executed and verified

Progress: [██░░░░░░░░] 25%

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

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 2 readiness]: Cloud hosting platform not yet selected — must verify HTTP 206 range request support via `curl -I --range 0-100 <url>` before writing fetch code. AWS S3 is confirmed; Zenodo and HuggingFace Hub need runtime verification.
- [Phase 2 readiness]: `httptest2` interceptability of Arrow's libcurl calls is unconfirmed — R-level HTTP mocking may not reach the C++ layer. Fallback is `skip_on_ci()` for network tests + RDS fixture testing.

## Session Continuity

Last session: 2026-03-03
Stopped at: Phase 1 plan 01-01 complete. API_DESIGN.md written, verified, and committed. Phase 2 ready to plan.
Resume file: None
