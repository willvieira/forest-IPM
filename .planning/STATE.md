# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Any researcher can run a full IPM for any supported species at any location by installing the package — no local data setup, no manual sourcing, no cluster required.
**Current focus:** Phase 1 - API Design

## Current Position

Phase: 1 of 4 (API Design)
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-02-25 — Roadmap created

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**
- Total plans completed: 0
- Average duration: -
- Total execution time: 0 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| - | - | - | - |

**Recent Trend:**
- Last 5 plans: -
- Trend: -

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

### Pending Todos

None yet.

### Blockers/Concerns

- [Phase 2 readiness]: Cloud hosting platform not yet selected — must verify HTTP 206 range request support via `curl -I --range 0-100 <url>` before writing fetch code. AWS S3 is confirmed; Zenodo and HuggingFace Hub need runtime verification.
- [Phase 2 readiness]: `httptest2` interceptability of Arrow's libcurl calls is unconfirmed — R-level HTTP mocking may not reach the C++ layer. Fallback is `skip_on_ci()` for network tests + RDS fixture testing.

## Session Continuity

Last session: 2026-02-25
Stopped at: Roadmap created, files written. Phase 1 ready to plan.
Resume file: None
