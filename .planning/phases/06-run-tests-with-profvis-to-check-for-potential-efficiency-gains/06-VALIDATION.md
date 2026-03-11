---
phase: 6
slug: run-tests-with-profvis-to-check-for-potential-efficiency-gains
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-11
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat (via DESCRIPTION Suggests) |
| **Config file** | `tests/testthat.R` (exists) |
| **Quick run command** | `devtools::test()` |
| **Full suite command** | `devtools::check()` |
| **Estimated runtime** | ~60 seconds |

---

## Sampling Rate

- **After every task commit:** Run `devtools::test()`
- **After every plan wave:** Run `devtools::check()`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** ~60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 6-01-01 | 01 | 0 | baseline | regression | `source("profile_benchmark.R")` — baseline capture | ❌ W0 | ⬜ pending |
| 6-01-02 | 01 | 1 | profiling | manual | `profvis::profvis({ lambda(sp) })` — inspect HTML | ✅ | ⬜ pending |
| 6-02-01 | 02 | 2 | vectorize outer() | regression | `devtools::test()` + lambda equality assert | ❌ W0 | ⬜ pending |
| 6-02-02 | 02 | 2 | BA cache | regression | `devtools::test()` + lambda equality assert | ❌ W0 | ⬜ pending |
| 6-03-01 | 03 | 3 | suite passes | integration | `devtools::check()` | ✅ | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/testthat/test-regression-baselines.R` — regression assertions comparing lambda() and project() output to pre-optimization baseline
- [ ] `tests/testthat/fixtures/lambda_baseline.rds` — saved numeric baseline for QUERUB lambda before any optimization
- [ ] `tests/testthat/fixtures/project_baseline.rds` — saved numeric baseline for project() trajectory before any optimization

*Baseline capture must happen BEFORE any optimization changes are applied.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| profvis flame graph identifies outer() as hotspot | Profiling | Visual inspection of HTML profvis output | Run `profvis::profvis({ lambda(sp) }, interval=0.005)`, open HTML, confirm mkKernel() dominates |
| Speedup is measurable | Performance | Requires human judgment on tradeoff | Run microbenchmark before/after, confirm ≥10% speedup or document why no gain found |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
