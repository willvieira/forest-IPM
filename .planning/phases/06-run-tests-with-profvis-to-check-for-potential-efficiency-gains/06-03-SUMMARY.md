---
phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains
plan: "03"
subsystem: performance
tags: [IPM, kernel, optimization, vectorization, outer, truncnorm, dnorm, pnorm, microbenchmark]

requires:
  - phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains
    plan: "02"
    provides: "PROFILING.md with approved optimizations (Target 1 + Target 2), regression baselines in tests/testthat/fixtures/"

provides:
  - "Optimized R/kernel.R: outer() vectorization (Target 1) + truncnorm replacement (Target 2)"
  - "PROFILING.md Before/After Benchmark Results section with real measured speedups"
  - "benchmark_before_after.R: reproducible before/after microbenchmark script"
  - "tests/testthat/test-vectorized-kernel.R: 3x3 correctness guard for vectorized P and F matrices"
  - "truncnorm removed from DESCRIPTION Imports (dependency eliminated)"

affects:
  - "Future phases using mkKernel() or project() — same API, faster execution"
  - "Phase 5 comparison script (simulations/) — 100-yr runs now ~22s faster"

tech-stack:
  added: []
  patterns:
    - "Vectorized kernel: rep(meshpts,times=n)/rep(meshpts,each=n) + direct FUN call + matrix() replaces outer()"
    - "Truncated normal density: dnorm(x,mu,sig)/(1-pnorm(a,mu,sig)) is base-R equivalent of truncnorm::dtruncnorm(x,a,Inf,mu,sig)"
    - "outer() semantics: X_expanded=rep(X,times=n), Y_expanded=rep(Y,each=n) — key for correct each/times ordering"

key-files:
  created:
    - ".planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/benchmark_before_after.R"
    - "tests/testthat/test-vectorized-kernel.R"
  modified:
    - "R/kernel.R"
    - "DESCRIPTION"
    - "tests/testthat/test-package-structure.R"
    - ".planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/PROFILING.md"

key-decisions:
  - "outer() semantics: X_expanded=rep(X,times=n), Y_expanded=rep(Y,each=n); s1=rep(meshpts,times=n), s0=rep(meshpts,each=n)"
  - "Target 2 (truncnorm elimination) is the dominant gain source; Target 1 (outer vectorization) alone gives negligible speedup for P matrix"
  - "Remove truncnorm from DESCRIPTION Imports — dependency eliminated now that dnorm/pnorm are used directly"
  - "Regression baseline fixtures remain valid — output is numerically identical to pre-optimization (tested to 1e-10)"

patterns-established:
  - "Vectorized kernel pattern: rep(x,times=n)/rep(x,each=n) + FUN(s1,s0,...) + matrix(result,n,n)"
  - "Truncated normal: dnorm/pnorm replacement formula for truncnorm::dtruncnorm with fixed lower bound"

requirements-completed: []

duration: ~120min (including 3x 5-min regression test runs ~300s each)
completed: 2026-03-12
---

# Phase 06 Plan 03: Apply Approved Kernel Optimizations Summary

**F matrix 2.4x faster by replacing truncnorm::dtruncnorm with dnorm/pnorm; combined mkKernel() 1.4x faster (761ms -> 542ms), 100-year run reduced from ~76s to ~54s**

## Performance

- **Duration:** ~120 min (dominated by 3 regression test runs at ~5 min each)
- **Started:** 2026-03-12
- **Completed:** 2026-03-12
- **Tasks:** 2 of 3 (Task 3 is checkpoint:human-verify — awaiting human approval)
- **Files modified:** 5

## Accomplishments

- Applied Target 1 (vectorize outer()) and Target 2 (replace truncnorm::dtruncnorm with dnorm/pnorm) to R/kernel.R
- Regression tests confirm lambda() and project() output numerically identical to pre-optimization baselines (tolerance 1e-10)
- Measured 2.4x speedup for F matrix (427.5ms -> 174.6ms) and 1.4x for full mkKernel() (761ms -> 542ms)
- Eliminated truncnorm package dependency from DESCRIPTION (closes PITFALLS.md Pitfall 11)
- Added 3x3 micro-kernel tests guarding against transposition bugs in vectorized implementation
- Documented before/after benchmark results in PROFILING.md with no placeholder strings

## Task Commits

Each task was committed atomically:

1. **TDD RED — 3x3 vectorized kernel correctness tests** - `fd3a0e5` (test)
2. **Task 1: Apply Target 1 + Target 2 optimizations** - `b275819` (feat)
3. **[Auto-fix Rule 1] Fix outer() each/times ordering** — discovered during regression test failure; corrected in same commit as `b275819` (fix within task)
4. **Task 2: Before/after benchmark + PROFILING.md update** - `05df204` (feat)
5. **[Auto-fix Rule 2] Remove truncnorm from DESCRIPTION** - `c3c3680` (fix)

## Files Created/Modified

- `/Users/wvieira/GitHub/forest-IPM/R/kernel.R` — Target 1 (outer->vectorized) + Target 2 (truncnorm->dnorm/pnorm) applied
- `/Users/wvieira/GitHub/forest-IPM/tests/testthat/test-vectorized-kernel.R` — 3x3 micro-kernel correctness tests (P and F matrices)
- `/Users/wvieira/GitHub/forest-IPM/.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/benchmark_before_after.R` — reproducible before/after benchmark script
- `/Users/wvieira/GitHub/forest-IPM/.planning/phases/06-run-tests-with-profvis-to-check-for-potential-efficiency-gains/PROFILING.md` — Before/After Benchmark Results section added
- `/Users/wvieira/GitHub/forest-IPM/DESCRIPTION` — truncnorm removed from Imports
- `/Users/wvieira/GitHub/forest-IPM/tests/testthat/test-package-structure.R` — truncnorm removed from expected imports list

## Decisions Made

- **outer() argument ordering**: `outer(X, Y, FUN)` passes `X_expanded=rep(X,times=n)` as FUN's first arg and `Y_expanded=rep(Y,each=n)` as second. In mkKernel: `s1=rep(meshpts,times=n)` (size_t1), `s0=rep(meshpts,each=n)` (size_t0). The initial implementation had this swapped — caught by regression test failure (>1e-3 difference vs baseline).
- **Target 2 is the dominant gain**: Replacing truncnorm::dtruncnorm with dnorm/pnorm gave 2.4x speedup for F matrix. Target 1 (vectorize outer) alone gave ~1.0x for P matrix — the R-to-C dispatch overhead was minimal for P_xEC.
- **Remove truncnorm dependency**: Since ingrowth_lk now uses dnorm/pnorm, truncnorm has no remaining usages in the package. Removed from DESCRIPTION per Rule 2 (correctness of declared dependencies).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed outer() each/times argument ordering in mkKernel() and test file**
- **Found during:** Task 1 — regression tests after first implementation attempt
- **Issue:** Initial vectorized implementation had `s1=rep(meshpts,each=n)` and `s0=rep(meshpts,times=n)` — this is transposed relative to outer()'s X/Y expansion semantics. Caused mean relative lambda difference of ~6.6e-4 vs baseline.
- **Fix:** Corrected to `s1=rep(meshpts,times=n)` and `s0=rep(meshpts,each=n)` in both R/kernel.R and test-vectorized-kernel.R
- **Files modified:** R/kernel.R, tests/testthat/test-vectorized-kernel.R
- **Verification:** Regression tests pass (FAIL 0, PASS 2), 3x3 tests confirm exact match with outer()
- **Committed in:** `b275819` (integrated into task commit)

**2. [Rule 2 - Missing Critical] Removed truncnorm from DESCRIPTION Imports**
- **Found during:** Task 1 post-implementation — noticed truncnorm still in DESCRIPTION after removing all usage
- **Issue:** After replacing `truncnorm::dtruncnorm` with `dnorm/pnorm`, `truncnorm` had zero remaining runtime usage but remained in DESCRIPTION Imports and package-structure test
- **Fix:** Removed `truncnorm (>= 1.0-9)` from DESCRIPTION Imports; updated required packages list in test-package-structure.R
- **Files modified:** DESCRIPTION, tests/testthat/test-package-structure.R
- **Verification:** test-package-structure tests pass (PASS 18, FAIL 0)
- **Committed in:** `c3c3680`

---

**Total deviations:** 2 auto-fixed (1 bug, 1 missing critical)
**Impact on plan:** Both auto-fixes were essential — the ordering bug would have silently broken the math, and the stale dependency declaration is a package correctness issue.

## Issues Encountered

- Initial outer() vectorization had swapped each/times — discovered via regression test failure (mean relative diff 6.6e-4, far outside 1e-10 tolerance). Careful reading of R's outer() documentation confirmed correct semantics.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- All regression tests pass (FAIL 0, WARN 0, SKIP 0, PASS 2 for regression; PASS 18 for package-structure)
- Phase 06 is complete pending human verification of this checkpoint
- 100-year runs now ~54s (down from ~76s), a ~22s saving per run
- If further speedup desired: Target 3 (BA caching, ~10%) and Target 4 (C++ migration, pending Xcode license) remain available

---
*Phase: 06-run-tests-with-profvis-to-check-for-potential-efficiency-gains*
*Completed: 2026-03-12*
