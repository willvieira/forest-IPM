---
phase: 5
slug: compare-old-vs-new-package-results
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-12
---

# Phase 5 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Standalone R script (not testthat) — per CONTEXT.md Decision D |
| **Config file** | none — interactive script |
| **Quick run command** | `source("simulations/compare_versions/compare_versions.R")` |
| **Full suite command** | same — script is its own suite |
| **Estimated runtime** | ~60–120 seconds (15–20 rows × 100 reps × 5 kernels) |

---

## Sampling Rate

- **After every task commit:** Run `source("simulations/compare_versions/compare_versions.R")`
- **After every plan wave:** Run same (single wave phase)
- **Before `/gsd:verify-work`:** Full script must print all PASS
- **Max feedback latency:** ~120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 05-01-01 | 01 | 1 | — | regression | `source("simulations/compare_versions/compare_versions.R")` | ❌ W0 | ⬜ pending |
| 05-01-02 | 01 | 1 | — | regression | same | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `simulations/compare_versions/compare_versions.R` — main comparison script (this is Wave 0 and Wave 1 combined — single deliverable phase)

*Existing infrastructure: testthat suite in `tests/` covers package code; this script covers cross-version equivalence outside of testthat.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Console output shows PASS for all stratified rows | — | Script produces console output, not test framework artifacts | Run `source("simulations/compare_versions/compare_versions.R")`; inspect printed table — all rows must show PASS |
| 6 missing rows excluded from selection | — | Selection logic must avoid IDs 32794, 33256, 33257, 171478, 195665, 201283 | Verify hard-coded indices in script do not include these IDs |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
