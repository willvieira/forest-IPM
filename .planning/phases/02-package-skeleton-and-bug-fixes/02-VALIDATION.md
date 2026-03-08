---
phase: 02
slug: package-skeleton-and-bug-fixes
status: partial
nyquist_compliant: false
wave_0_complete: true
created: 2026-03-07
---

# Phase 02 — Validation Strategy

> Per-phase validation contract. Reconstructed from PLAN and SUMMARY artifacts (State B).

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | testthat 3.x |
| **Config file** | `tests/testthat.R` |
| **Quick run command** | `devtools::test()` |
| **Full suite command** | `devtools::test()` + `devtools::check()` |
| **Estimated runtime** | ~30 seconds |
| **Test files** | `tests/testthat/test-bug-fixes.R`, `tests/testthat/test-constructors.R`, `tests/testthat/test-package-structure.R` |

---

## Sampling Rate

- **After every task commit:** Run `devtools::test()`
- **After every plan wave:** Run `devtools::test()` + `devtools::check()`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** ~30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File | Status |
|---------|------|------|-------------|-----------|-------------------|------|--------|
| 02-01-01 | 01 | 1 | BUG-01 | unit | `devtools::test(filter="bug-fixes")` | `test-bug-fixes.R` | ✅ green |
| 02-01-01 | 01 | 1 | BUG-01 (static) | static | `devtools::test(filter="bug-fixes")` | `test-bug-fixes.R` | ✅ green |
| 02-01-02 | 01 | 1 | BUG-02 | static | `devtools::test(filter="bug-fixes")` | `test-bug-fixes.R` | ✅ green |
| 02-01-03 | 01 | 1 | BUG-03 | static | `devtools::test(filter="bug-fixes")` | `test-bug-fixes.R` | ✅ green |
| 02-02-01 | 02 | 2 | stand() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-02-01 | 02 | 2 | species_model() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-02-01 | 02 | 2 | parameters() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-02-01 | 02 | 2 | env_condition() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-02-01 | 02 | 2 | control() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-02-01 | 02 | 2 | supported_species() | unit | `devtools::test(filter="constructors")` | `test-constructors.R` | ✅ green |
| 02-03-01 | 03 | 3 | PKG-02 (DESCRIPTION) | static | `devtools::test(filter="package-structure")` | `test-package-structure.R` | ✅ green |
| 02-03-02 | 03 | 3 | PKG-03 (NAMESPACE) | static | `devtools::test(filter="package-structure")` | `test-package-structure.R` | ✅ green |
| 02-03-02 | 03 | 3 | PKG-03 (globals.R) | static | `devtools::test(filter="package-structure")` | `test-package-structure.R` | ✅ green |
| 02-03-03 | 03 | 3 | PKG-04 | compilation | `devtools::load_all()` | manual | ⚠️ manual |
| 02-03-03 | 03 | 3 | PKG-01 | install | `devtools::install_github()` | manual | ⚠️ manual |
| 02-03-03 | 03 | 3 | PKG-05 (full check) | integration | `devtools::check()` | manual | ⚠️ manual |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ manual*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| `devtools::load_all()` compiles `src/eigen.cpp` without errors | PKG-04 | Requires C++ build tools (Xcode) and cannot be run inside an automated testthat runner | Run `devtools::load_all()` after accepting Xcode license; expect no compilation errors |
| `max(getEigenValues(K))` equals `Re(eigen(K)$values[1])` within 1e-10 at runtime | BUG-01 (runtime) | Requires compiled shared library | `K <- matrix(c(0.8,0.1,0.05,0.9),nrow=2); isTRUE(all.equal(max(getEigenValues(K)), Re(eigen(K)$values[1]), tolerance=1e-10))` |
| `devtools::check()` exits with 0 errors, 0 warnings | PKG-05 | Full check requires C++ compilation; R-level --no-install check was 0/0/0 (documented in 02-03-SUMMARY) | `devtools::check()` after Xcode license acceptance |
| `devtools::install_github("wvieira/forest-IPM")` succeeds on fresh R | PKG-01 | Requires C++ compilation + network + clean R environment | Run from fresh R session in a temp directory |
| `lambda()` end-to-end with real parameters | — | Requires local RDS parameter files (Phase 3 scope) | `lambda(species_model(stand(...)), parameters(...), ..., env_condition(...))` once Phase 3 data is available |

---

## Notes

### C++ Compilation Fix (2026-03-07)

During validation, a compilation error was discovered: `RcppExports.cpp` used `Map<MatrixXd>` and `VectorXd` without `Eigen::` namespace qualifiers, because the `using Eigen::Map` declarations in `eigen.cpp` don't carry over to the generated file.

**Fix applied:**
- `src/eigen.cpp` updated to use fully-qualified `Eigen::VectorXd`, `Eigen::Map<Eigen::MatrixXd>`, `Eigen::EigenSolver<Eigen::MatrixXd>` in the function signature
- `Rcpp::compileAttributes()` regenerated `src/RcppExports.cpp` with correct qualified types
- `devtools::load_all()` now compiles successfully: `Match: TRUE` confirmed

---

## Validation Audit 2026-03-07

| Metric | Count |
|--------|-------|
| Gaps found | 15 |
| Automated tests written | 59 (passing) |
| Escalated to manual-only | 6 |
| C++ compilation fix applied | 1 |

---

_Validated: 2026-03-07_
_Validator: Claude (gsd-validate-phase)_
