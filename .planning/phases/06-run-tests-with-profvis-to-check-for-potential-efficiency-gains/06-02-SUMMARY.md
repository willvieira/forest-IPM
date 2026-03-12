---
phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains
plan: "02"
subsystem: testing
tags: [profvis, profiling, performance, ipm, r-package, microbenchmark, flamegraph]

# Dependency graph
requires:
  - phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains
    plan: "01"
    provides: Pre-optimization regression fixtures + tidyverse-revamped R/ source files
provides:
  - profvis flamegraph HTML artifact (QUERUB single-species 5-year project())
  - profvis multi-species flamegraph HTML artifact (QUERUB+QUEPRI 3-year project())
  - PROFILING.md with numeric hotspot analysis and three optimization targets
  - profile_benchmark.R standalone script for reproducible profiling
affects: [06-03-optimization, any future performance work]

# Tech tracking
tech-stack:
  added:
    - profvis (profiling — not in DESCRIPTION, dev-only)
    - microbenchmark (benchmarking — not in DESCRIPTION, dev-only)
    - htmlwidgets (saveWidget for flamegraph HTML — not in DESCRIPTION, dev-only)
  patterns:
    - "Profiling scripts live in .planning/phases/ — not shipped with package"
    - "htmlwidgets::saveWidget(selfcontained=TRUE) for portable flamegraph artifacts"
    - "Component microbenchmark: isolate outer(P), outer(F), BA metrics separately to quantify each cost"

key-files:
  created:
    - .planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profile_benchmark.R
    - .planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/PROFILING.md
    - .planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profvis_output.html
    - .planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profvis_multi_output.html
  modified:
    - .gitignore

key-decisions:
  - "F matrix (outer + ingrowth_lk + truncnorm::dtruncnorm) is 2.5x more expensive than P matrix — priority order changed from plan: Target 2 (replace truncnorm) must be done alongside Target 1 (vectorize outer)"
  - "truncnorm::dtruncnorm is largest leaf hotspot at 46% inclusive time — replace with dnorm/pnorm direct formula (math equivalent, base R vectorizes natively)"
  - "Force-add profvis HTML artifacts with .gitignore !exception pattern — same approach as RDS fixtures in 06-01"
  - "Profile at 5 timesteps (not 100) to get profvis data in under 10s — scales linearly so 100-year projection estimate (~76s) derived from microbenchmark"

patterns-established:
  - "Profile at representative but short run (5 timesteps) for profvis; use microbenchmark for absolute timing at 10+ reps"
  - "Decompose mkKernel() cost into P-outer, F-outer, BA-metrics to identify which component to optimize first"

requirements-completed: []

# Metrics
duration: 16min
completed: 2026-03-12
---

# Phase 06 Plan 02: Profvis Profiling Summary

**profvis flamegraphs and component microbenchmarks reveal F matrix (truncnorm::dtruncnorm) as dominant hotspot at 46% of runtime — Priority 1 for Plan 03 optimization alongside outer() vectorization**

## Performance

- **Duration:** ~16 min
- **Started:** 2026-03-12T01:45:15Z
- **Completed:** 2026-03-12T02:01:33Z
- **Tasks:** 2 of 2 (Task 2 checkpoint completed after human flamegraph review)
- **Files modified:** 5

## Accomplishments

- Wrote `profile_benchmark.R` — standalone profiling script covering QUERUB single-species and QUERUB+QUEPRI 2-species benchmarks via profvis and microbenchmark
- Saved `profvis_output.html` (603KB self-contained flamegraph, 5-year QUERUB project()) and `profvis_multi_output.html` (3-year 2-species)
- Produced `PROFILING.md` with numeric hotspot breakdown: mkKernel() 86% of time, F matrix 64%, P matrix 26%, BA metrics 10%, truncnorm::dtruncnorm 46% inclusive
- Key unexpected finding: F matrix is 2.5× more expensive than P matrix — changes optimization priority order vs. plan expectations

## Task Commits

1. **Task 1: Write and run profiling script, produce PROFILING.md** - `2874adc` (feat)
2. **Task 2: Checkpoint — record approved optimization targets in PROFILING.md** - `7d615d0` (docs)

## Files Created/Modified

- `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profile_benchmark.R` - Standalone profiling script; benchmarks mkKernel components, runs profvis, saves flamegraphs
- `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/PROFILING.md` - Hotspot analysis with numeric data, call stack breakdown, 4 optimization targets with evidence-based assessments
- `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profvis_output.html` - 603KB self-contained flamegraph (QUERUB 5-year)
- `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profvis_multi_output.html` - 2-species flamegraph (QUERUB+QUEPRI 3-year)
- `.gitignore` - Added `!.planning/phases/*/profvis_output.html` and `!.planning/phases/*/profvis_multi_output.html` exceptions

## Decisions Made

- The plan anticipated P matrix as dominant; profiling revealed F matrix (via `truncnorm::dtruncnorm`) is 2.5× more expensive. Priority updated: Target 1 (vectorize outer) and Target 2 (replace truncnorm with dnorm/pnorm) must be co-applied in Plan 03 since Target 2 benefits from vectorized input.
- Profiled at 5 timesteps rather than 100 to get profvis data quickly (~3s vs ~76s); extrapolated 100-year costs from microbenchmark median.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] .gitignore blocked committing *.html profiling artifacts**
- **Found during:** Task 1 (commit step)
- **Issue:** `*.html` blanket rule in .gitignore prevented `git add profvis_output.html`
- **Fix:** Added `!.planning/phases/*/profvis_output.html` and `!.planning/phases/*/profvis_multi_output.html` exception lines to .gitignore
- **Files modified:** .gitignore
- **Verification:** Both HTML files staged and committed successfully
- **Committed in:** 2874adc (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Fix necessary to commit profiling artifacts. No scope creep.

## Issues Encountered

- **htmlwidgets::saveWidget creates companion `_files/` directory** despite `selfcontained=TRUE`. The HTML embeds all JS/CSS (verified: no external references, 603KB). The `_files/` directory is a side effect that does not need to be committed.

## Next Phase Readiness

- Profiling complete and human-reviewed
- PROFILING.md has numeric evidence and explicit approval for Plan 03 optimization decisions
- Plan 03 is cleared to apply: (1) vectorize outer() in mkKernel(), (2) replace truncnorm::dtruncnorm with dnorm/pnorm in ingrowth_lk()
- Target 3 (BA caching) and Target 4 (C++) are deferred — evaluate after Targets 1+2

## Self-Check: PASSED

- FOUND: `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profvis_output.html` (603.4 KB)
- FOUND: `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/PROFILING.md` (9.3 KB)
- FOUND: `.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/profile_benchmark.R`
- FOUND commit 2874adc (Task 1)
- FOUND commit 7d615d0 (Task 2 — Approved Optimizations section in PROFILING.md)
- No unfilled {placeholders} in PROFILING.md
- PROFILING.md contains ## Approved Optimizations section with explicit go/no-go for each target

---
*Phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains*
*Completed: 2026-03-12*
