# Project Research Summary

**Project:** forest-IPM
**Domain:** Scientific R package — Bayesian hierarchical Integral Projection Model for tree population dynamics
**Researched:** 2026-02-25
**Confidence:** HIGH (infrastructure, stack, architecture) / MEDIUM (cloud data layer specifics)

## Executive Summary

The forest-IPM project is converting a working simulation codebase into a distributable R package that allows researchers to run Bayesian hierarchical Integral Projection Models for 31 tree species in eastern North America without requiring a local 15 GB parameter dataset. The recommended approach is a four-layer architecture: (1) a standard R package infrastructure scaffold (DESCRIPTION, NAMESPACE, roxygen2, testthat), (2) an S3 domain object layer wrapping the existing math engine, (3) a remote data layer using Arrow + Parquet with HTTP predicate pushdown to fetch only the needed posterior draws on demand, and (4) a high-level user API (`run_ipm()`) that orchestrates the layers. The existing C++ eigenvalue solver (Rcpp + RcppEigen) and the core vital rate and kernel functions are already correct and require no rewriting — the packaging work is additive infrastructure.

The primary technical risk is the cloud data layer. The promise of the package — "no local data setup" — depends entirely on Parquet files being partitioned by species (hive partitioning) and hosted on a server that supports HTTP 206 range requests. If either condition fails, users download 15 GB per call, making the package unusable. AWS S3 public buckets are the highest-confidence hosting option; Zenodo and HuggingFace Hub are acceptable fallbacks pending verification with `curl -I --range`. A secondary risk is the existing codebase's script-style patterns (hardcoded relative paths, undeclared tidyverse imports, `%>%` / `|>` mixing, `exists('fct')` global state, deprecated `as_vector()`, `SelfAdjointEigenSolver` used on a non-symmetric matrix) — these will all break `R CMD check` and must be resolved in Phase 1 before any other work.

The roadmap should proceed in strict dependency order: package skeleton first (nothing else works until `devtools::load_all()` succeeds and `R CMD check` produces zero errors), then S3 classes (define the contracts all other layers return), then the remote data layer (the core value proposition), then the user-facing API and documentation. Testing infrastructure (synthetic fixtures, not real data) must be established in Phase 1 — tests written against real parameter files will fail on every CI run.

---

## Key Findings

### Recommended Stack

The toolchain is unambiguous: `devtools` + `usethis` + `roxygen2` + `testthat` (edition 3) is the universal R package standard. `usethis::create_package()` + `usethis::use_rcpp_eigen()` will scaffold the correct DESCRIPTION/NAMESPACE/CI files in minutes. The existing `src/eigen.cpp` already uses correct `[[Rcpp::export]]` and `[[Rcpp::depends(RcppEigen)]]` attributes — it will compile without modification once DESCRIPTION is in place.

The most consequential stack decision is `arrow` (>= 14.0.0) for remote Parquet access. It is the only mature R library that supports lazy Parquet reads with predicate pushdown over HTTP. The `open_dataset()` + `filter()` + `collect()` pattern is the correct implementation: filter before materializing. Parquet files must be written with hive partitioning by species (`write_dataset(df, path, partitioning = "species_id")`) to make pushdown effective — this is a write-time decision that cannot be changed after cloud upload without re-partitioning.

**Core technologies:**
- `devtools` / `usethis` / `roxygen2`: R package scaffold, documentation, NAMESPACE — unambiguous standard since 2016
- `testthat` (edition 3): unit testing — with `httptest2` or RDS fixtures for network-free CI
- `Rcpp` + `RcppEigen`: C++ eigenvalue solver — already integrated, keep as-is
- `arrow` (>= 14.0.0): remote Parquet reads with HTTP predicate pushdown — the only viable option for partial 15 GB dataset access
- Base R S3 classes: `ipm_params`, `ipm_kernel`, `ipm_population`, `ipm_projection` — idiomatic, serializable, no external dependencies
- `tools::R_user_dir()`: XDG-compliant local cache — CRAN-blessed location for persistent parameter cache
- GitHub Actions (`r-lib/actions`): CI/CD for `R CMD check` on Linux + macOS

### Expected Features

**Must have (table stakes for v1):**
- `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `NEWS.md` — no package without these; `R CMD check` fails
- Rcpp registration + zero-error `R CMD check` — compilation must succeed on install
- roxygen2 docs on all exported functions — `?function_name` must work
- Bundled `sysdata.rda` with posterior means for 2 species — enables tests and examples without network access
- `testthat` suite covering kernel assembly, lambda computation, `init_pop` edge cases — regression safety net
- Cloud Parquet + Arrow lazy reads replacing local RDS path in `getPars_sp()` — the core value proposition
- `run_ipm()` high-level orchestration function — users should not need to know the internal call order
- S3 classes for `ipm_kernel`, `ipm_params`, `ipm_population` with `print()` / `summary()` methods
- One getting-started vignette — minimum documentation for a scientific package
- `README.md` with install instructions and a quick example
- Local parameter cache via `tools::R_user_dir()` — fetching on every call is a non-starter
- Informative error messages on network unavailability — HPC clusters have no internet

**Should have (differentiators):**
- `pkgdown` site (GitHub Pages) — searchable function reference and vignettes
- Full posterior draw access (`n_draws` argument) — posterior mean alone is scientifically incomplete
- `plot.ipm_kernel()` method delegating to `matrix_image.R` — researchers need to visually inspect kernels
- `supported_species()` helper — users need to know which of the 31 species are available
- GitHub Actions CI (Linux + macOS) — signals package health; catches Rcpp regression silently
- `CITATION` file + `citation("forestipm")` — scientific packages must be citable
- Parameter version pinning in cloud paths — reproducibility across years
- Input validation with climate training-range warnings

**Defer to v2+:**
- ML decoder for unobserved sites (lat/lon to parameters) — novel research contribution, very high complexity
- Multiple vignettes (uncertainty propagation, sensitivity analysis)
- CRAN submission — premature until cloud data access pattern is proven stable
- Community/multi-species IPM — single-species API must stabilize first
- `lifecycle` badges — premature until API is stable
- Windows binary support — RcppEigen CI complexity deferred per PROJECT.md

### Architecture Approach

The architecture adds three production layers on top of an existing, working simulation engine. The existing math (vital rates, kernel assembly, competition, C++ eigenvalues) is correct and should not be rewritten. The work is: (1) adding a remote data layer that replaces local RDS reads with Arrow Parquet fetches partitioned by species, (2) adding an S3 domain object layer that makes IPM objects first-class R citizens with `print()`, `plot()`, `summary()` dispatch, and (3) adding package infrastructure that makes it installable via `devtools::install_github()`. Data flows strictly top-down: User API → S3 domain objects → {Remote Data, IPM Engine, C++ Acceleration}.

**Major components:**
1. **Package infrastructure** (`DESCRIPTION`, `NAMESPACE`, `src/` Rcpp compilation) — foundation; nothing else works without it
2. **Remote data layer** (`R/params.R` + `R/cache.R`) — Arrow `open_dataset()` with HTTP predicate pushdown; XDG cache via `tools::R_user_dir()`
3. **S3 domain objects** (`R/classes.R`) — thin wrappers over existing list structures; `ipm_params`, `ipm_kernel`, `ipm_population`, `ipm_projection`
4. **IPM engine** (existing `R/vital_rates.R`, `R/kernel.R`, `R/BasalArea_competition.R`, `src/eigen.cpp`) — no rewriting needed; wire into S3 layer
5. **User API** (`R/ipm.R`) — `run_ipm()`, `project_population()`, `lambda()` orchestrate the above
6. **Visualization** (`R/plot.R`) — S3 `plot()` methods delegating to existing `matrix_image.R`

Parquet layout must use hive partitioning by `species_id` so Arrow uses partition columns as zero-overhead filters. Within each species file, sort by `iter` and set `chunk_size = 100` at write time for efficient draw-level pushdown.

### Critical Pitfalls

1. **NAMESPACE pollution from undeclared tidyverse imports** — `params.R` uses `map`, `filter`, `pivot_wider`, `mutate`, `str_replace`, etc. without declaring them in DESCRIPTION `Imports:`. Run `devtools::check()` in the first commit and fix every NOTE before proceeding. Add `dplyr`, `tidyr`, `purrr`, `stringr` to `Imports:` and use `@importFrom` tags.

2. **`SelfAdjointEigenSolver` on a non-symmetric matrix** — The IPM kernel K = P + F is NOT symmetric; the self-adjoint solver returns silently wrong eigenvalues. Replace with `EigenSolver<MatrixXd>` for the full K matrix, or expose a high-level `compute_lambda(K)` that uses base R `eigen()` internally. Validate before wiring into the public API.

3. **Parquet row group misalignment defeats predicate pushdown** — If all 31 species share a single Parquet file with large mixed row groups, `filter(species == "ABBA")` downloads nearly the full 15 GB. Must use `write_dataset(df, path, partitioning = "species_id")` before cloud upload — the layout cannot be changed cheaply afterward.

4. **Cloud host without HTTP range request support** — Arrow's partial reads require HTTP 206. Verify with `curl -I --range 0-100 <url>` before writing any fetch code. AWS S3 public buckets are confirmed; Zenodo and HuggingFace need verification.

5. **Tests that call `getPars_sp()` against real files fail on CI** — Create `tests/testthat/fixtures/` with minimal synthetic parameter sets in Phase 1. Gate any integration tests behind `testthat::skip_if_not(has_remote_data())`. Never let a test assume the 15 GB data exists.

---

## Implications for Roadmap

Based on combined research, the build order is strict — each phase has hard dependencies on the previous one. There is no parallelism between the core phases.

### Phase 1: Package Skeleton and Code Cleanup

**Rationale:** Nothing else compiles or installs until `devtools::load_all()` succeeds and `R CMD check` produces zero errors/warnings/NOTEs. All five of the critical Phase 1 pitfalls (NAMESPACE pollution, `%>%` / `|>` mixing, `exists('fct')` global state, `as_vector()` deprecation, `matrix.image` naming, `truncnorm` hard dependency, hardcoded paths) will cause immediate CI failure and must be resolved before any feature work begins.

**Delivers:** An installable R package skeleton with the existing math engine working through `devtools::load_all()`. `R CMD check` passes with zero errors/warnings. Rcpp compilation succeeds. Synthetic test fixtures established.

**Addresses:** All table-stakes package infrastructure features (DESCRIPTION, NAMESPACE, LICENSE, Rcpp registration, usethis scaffold).

**Avoids:** Pitfalls 1, 7, 8, 9, 10, 11, 12, 13, 14 — all Phase 1 code cleanup issues.

**Research flag:** Standard, well-documented patterns. No additional research needed.

### Phase 2: S3 Domain Object Layer

**Rationale:** S3 classes define the contracts that the data layer (Phase 3) and user API (Phase 4) must return. Building classes before the remote data layer means the data layer can return properly typed objects from day one, and tests can validate object contracts rather than raw list structure.

**Delivers:** `ipm_params`, `ipm_kernel`, `ipm_population`, `ipm_projection` S3 classes with constructors, validators, and `print()` / `summary()` / `plot()` methods. Unit tests validating class dispatch.

**Addresses:** S3 class features for all core objects; `plot.ipm_kernel()` delegating to matrix visualization.

**Avoids:** Pitfall 2 (S3 method export without NAMESPACE registration — inspect generated NAMESPACE after every `devtools::document()` call).

**Research flag:** Standard patterns. No additional research needed.

### Phase 3: Remote Data Layer (Cloud Parquet + Local Cache)

**Rationale:** This is the core value proposition of the package and the highest-risk phase. The Parquet layout decision (partitioning strategy, row group sizing, cloud host selection) must be made and verified before any fetch code is written — repartitioning after cloud upload is expensive. Verify HTTP range request support first with `curl` before touching Arrow code.

**Delivers:** `getPars_sp()` / `get_params()` backed by Arrow `open_dataset()` with HTTP predicate pushdown. Local cache via `tools::R_user_dir()`. Parquet files converted from existing RDS and uploaded with hive partitioning by `species_id`. Bundled `sysdata.rda` with posterior means for 2 species for offline use and tests.

**Addresses:** Cloud parameter fetching without 15 GB local dataset; local caching; informative error messages on network unavailability; `supported_species()` helper.

**Avoids:** Pitfalls 3, 4, 5 (Parquet chunking alignment, cloud host range request support, correct Arrow filter-before-collect pattern). Pitfall 6 (tests against real files) — handled by sysdata.rda fixtures from Phase 1.

**Research flag:** NEEDS RESEARCH. The specific Arrow + HTTP URL behavior for non-S3 hosts (HuggingFace Hub dataset API, Zenodo direct file URLs) needs verification. Cloud host selection and cost model need decision. The `httptest2` / arrow interaction for CI mocking needs a proof-of-concept.

### Phase 4: High-Level User API and Documentation

**Rationale:** Only buildable once all three prior layers are functional. `run_ipm()` orchestrates: parameter fetch (Phase 3) → kernel assembly (Phase 1 engine) → eigenvalue computation (Phase 1 C++) → projection loop. Documentation (vignettes, pkgdown) depends on having a stable API.

**Delivers:** `run_ipm()`, `project_population()`, `lambda()` public API. Eigenvalue computation fixed for non-symmetric matrices (Pitfall 3). Full roxygen2 documentation on all exported functions. Getting-started vignette. `README.md`. `CITATION` file. `pkgdown` site on GitHub Pages.

**Addresses:** `run_ipm()` high-level function; pipe-friendly API; `supported_species()` helper; session info in vignettes; DOI citation.

**Avoids:** Pitfall 15 (eigenvalue discrepancy — expose `compute_lambda()` wrapper, not raw `getEigenValues()`).

**Research flag:** Standard patterns for API design and pkgdown. No additional research needed.

### Phase 5: CI/CD, Testing Hardening, and Quality Signals

**Rationale:** Testing infrastructure should be built throughout all phases (synthetic fixtures in Phase 1, class tests in Phase 2, data layer tests in Phase 3), but a dedicated phase ensures full coverage, CI stability, and quality signals before any external sharing or publication.

**Delivers:** GitHub Actions matrix (ubuntu-latest, macOS-latest, R release + devel). Full `testthat` suite with network-gated integration tests. Code coverage via `covr` + codecov badge. `pkgdown` site deployment in CI.

**Addresses:** CI/CD; code coverage reporting; contributor guide.

**Research flag:** Standard patterns. No additional research needed.

### Phase Ordering Rationale

- Phase 1 must be first: Rcpp compilation is a hard dependency for `getEigenValues()`, which `run_ipm()` calls. `devtools::load_all()` failing blocks everything.
- Phase 2 before Phase 3: S3 class constructors must exist before the data layer can return typed `ipm_params` objects. Defining contracts first prevents refactoring later.
- Phase 3 before Phase 4: `run_ipm()` calls `get_params()` — cannot orchestrate what doesn't exist.
- Phase 5 throughout and last: Test infrastructure is built incrementally per phase; the dedicated phase hardens coverage and adds CI signals.

### Research Flags

Phases needing deeper research during planning:
- **Phase 3 (Remote Data Layer):** Arrow + non-S3 HTTP host behavior needs proof-of-concept verification before committing to any specific hosting platform or Parquet URL scheme. `httptest2` + arrow CI mocking needs a feasibility test.

Phases with standard, well-documented patterns (no additional research needed):
- **Phase 1:** R package structure is fully documented in r-pkgs.org and Writing R Extensions.
- **Phase 2:** S3 class patterns for scientific R packages are well-established (ipmr, popdemo, tidybayes).
- **Phase 4:** roxygen2, vignettes, pkgdown, and GitHub Actions are fully documented.
- **Phase 5:** GitHub Actions r-lib/actions and testthat edition 3 are well-documented.

---

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | devtools/usethis/roxygen2/testthat toolchain is the unambiguous R standard. Rcpp + RcppEigen already integrated correctly. S3 class system is idiomatic. Arrow is the correct library choice. |
| Features | HIGH | R package infrastructure requirements (CRAN policies, r-pkgs.org) are stable and well-documented. Feature scope aligns with established scientific package patterns (ipmr, popdemo, tidybayes). |
| Architecture | HIGH | Build order, component boundaries, and data flow are derived from direct codebase analysis and standard R package patterns. Parquet partitioning strategy is well-understood. |
| Pitfalls | HIGH (codebase-derived) / MEDIUM (cloud layer) | Six of fifteen pitfalls are confirmed via direct source file analysis. Cloud-related pitfalls (range request support, Parquet chunking) are based on Arrow documentation and standard HTTP behavior — needs runtime verification per host. |

**Overall confidence:** HIGH for the package infrastructure and domain object layers; MEDIUM for the cloud data layer specifics pending host verification.

### Gaps to Address

- **Arrow HTTP range request behavior per host:** Must verify with `curl -I --range 0-100 <url>` against candidate hosts (HuggingFace Hub, Zenodo) before writing fetch code. AWS S3 is confirmed; others are MEDIUM confidence.
- **`httptest2` interceptability of arrow's libcurl calls:** Arrow's internal HTTP is at C++ level. R-level HTTP mocking may not intercept it. Fallback is `skip_on_ci()` for network tests + RDS fixture testing for parsing logic. Needs a proof-of-concept in Phase 3.
- **`purrr::as_vector()` removal timeline:** Verify current purrr changelog to determine whether this is already removed or still a deprecation warning. Replace with `unlist(use.names = TRUE)` regardless.
- **Cloud hosting decision (S3 vs HuggingFace Hub vs Zenodo):** Cost model, institutional constraints, and access patterns (public vs. authenticated) need a decision before Phase 3. This affects the `open_dataset()` URL scheme and auth setup.
- **Companion paper DOI:** `CITATION` file depends on a published DOI. If the paper is not yet published, `CITATION` must be deferred or use a preprint DOI.

---

## Sources

### Primary (HIGH confidence — stable, well-documented standards)

- **R Packages (2e)** (Wickham & Bryan, r-pkgs.org) — package structure, DESCRIPTION, NAMESPACE, roxygen2, testthat, vignettes, pkgdown
- **Writing R Extensions** (CRAN official) — Rcpp/useDynLib, DESCRIPTION fields, `R CMD check` policy
- **CRAN Repository Policy** — mandatory requirements, examples runtime limits
- **Apache Arrow R documentation** (arrow.apache.org/docs/r) — `open_dataset()`, HTTP filesystems, predicate pushdown
- **Rcpp vignettes** (`vignette("Rcpp-package")`) — `compileAttributes`, `src/` structure, `LinkingTo`
- **Direct codebase analysis** (`R/params.R`, `src/eigen.cpp`, `R/kernel.R`, `.planning/codebase/CONCERNS.md`) — pitfalls confirmed from source

### Secondary (MEDIUM confidence — community consensus, multiple sources agree)

- **ipmr package** (Levin et al. 2021, MEE) — prior art for R-packaged IPMs; S3 kernel objects
- **popdemo package** — prior art for population matrix model S3 classes
- **tidybayes package** — prior art for Bayesian posterior access patterns
- **rOpenSci Packages Guide** (devguide.ropensci.org) — peer-review requirements (vignettes, tests, README, CITATION)
- **Arrow Parquet specification** — row group statistics and predicate pushdown mechanics

### Tertiary (LOW-MEDIUM confidence — needs runtime verification)

- **HuggingFace Hub range request support** — likely via CDN; needs `curl` verification
- **Zenodo direct file URL range requests** — supported per documentation; URL format needs testing
- **`httptest2` + arrow interaction** — needs proof-of-concept; R-level mocking may not reach libcurl

---

*Research completed: 2026-02-25*
*Ready for roadmap: yes*
