---
phase: 7
slug: add-github-actions-ci-with-code-coverage-and-pkgdown-site-with-readme-pointing-to-book-vignettes
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-13
---

# Phase 7 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat (existing) + devtools |
| **Config file** | `tests/testthat.R` |
| **Quick run command** | `Rscript -e "devtools::test()"` |
| **Full suite command** | `Rscript -e "devtools::test()"; Rscript -e "pkgdown::build_site()"` |
| **Estimated runtime** | ~60 seconds |

---

## Sampling Rate

- **After every task commit:** Run `Rscript -e "devtools::test()"` (verify existing tests still pass)
- **After every plan wave:** Run `Rscript -e "pkgdown::build_site()"` + `Rscript -e "covr::package_coverage()"`
- **Before `/gsd:verify-work`:** Full suite must be green + pkgdown builds locally
- **Max feedback latency:** 60 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 7-01-01 | 01 | 0 | DOC-04 | smoke | `Rscript -e "pkgdown::build_site()"` | ❌ W0 | ⬜ pending |
| 7-01-02 | 01 | 0 | TST-01 | smoke | push to main, observe workflow | ❌ W0 | ⬜ pending |
| 7-01-03 | 01 | 1 | TST-02 | measured | `Rscript -e "covr::package_coverage()"` | ✅ existing | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `_pkgdown.yml` — required before pkgdown smoke test can run
- [ ] `.github/workflows/R-CMD-check.yaml` — CI execution requires this file
- [ ] `.github/workflows/pkgdown.yaml` — pkgdown CI deploy requires this file

*Wave 0 creates all CI/infra files before remaining tasks can be validated.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| GitHub Actions CI runs green on push | TST-01 | Requires actual push + GitHub CI execution | Push to main, verify green check on GitHub |
| pkgdown site live at github.io URL | DOC-04 | Requires GitHub Pages enabled manually in repo settings | Enable Pages → gh-pages branch; verify site at https://willvieira.github.io/forest-IPM/ |
| Codecov badge visible in README | TST-02 | Requires manual repo activation on codecov.io | Activate repo on codecov.io, push coverage upload; verify badge renders |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 60s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
