---
phase: 4
slug: documentation-and-tests
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-13
---

# Phase 4 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat >= 3.0 |
| **Config file** | none — standard `tests/testthat/` layout |
| **Quick run command** | `devtools::test(filter = "workflow")` |
| **Full suite command** | `devtools::test()` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `devtools::test(filter = "workflow")`
- **After every plan wave:** Run `devtools::test()`
- **Before `/gsd:verify-work`:** `devtools::test()` zero failures + `devtools::run_examples()` zero errors
- **Max feedback latency:** ~30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 4-01-01 | 01 | 0 | DOC-02 | unit | `devtools::test(filter = "workflow")` | ❌ W0 | ⬜ pending |
| 4-01-02 | 01 | 0 | DOC-01 | unit | `devtools::test(filter = "workflow")` | ❌ W0 | ⬜ pending |
| 4-02-01 | 02 | 1 | DOC-02 | integration | `devtools::test()` | ❌ W0 | ⬜ pending |
| 4-02-02 | 02 | 1 | DOC-02 | coverage | `covr::package_coverage()` | ❌ W0 | ⬜ pending |
| 4-03-01 | 03 | 1 | DOC-01 | smoke | `devtools::run_examples()` | ❌ W0 | ⬜ pending |
| 4-04-01 | 04 | 2 | DOC-01 | smoke | `devtools::run_examples()` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/testthat/test-workflow.R` — stubs for DOC-02 workflow tests (stand → project)
- [ ] `R/project.R` — new `plot.ipm_projection(x, type = NULL, ...)` replacing old signature (prerequisite for DOC-01 plot examples)
- [ ] Delete obsolete test files: `test-bug-fixes.R`, `test-package-structure.R`, `test-vectorized-kernel.R`

*Wave 0 must be complete before wave 1 tasks can be verified.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `guide_IPM.qmd` renders without errors in the book repo | DOC-01 | Quarto render not part of R `devtools::test()` pipeline | `cd /Users/wvieira/GitHub/book_forest-demography-IPM && quarto render guide_IPM.qmd` |
| `devtools::run_examples()` passes zero-error on clean machine | DOC-01 | Requires full package load; CI gate | Run `devtools::run_examples()` in fresh R session after all @examples added |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
