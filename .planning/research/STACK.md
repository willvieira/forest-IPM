# Technology Stack

**Project:** forest-IPM
**Researched:** 2026-02-25
**Research mode:** Ecosystem — Stack dimension

> **Tool constraint note:** Web search and WebFetch were unavailable during this research session.
> All findings are from training data (knowledge cutoff August 2025). Confidence levels reflect
> this limitation. Versions marked with (*) should be verified against CRAN before committing to
> a DESCRIPTION file.

---

## Recommended Stack

### R Package Infrastructure

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| R | >= 4.1.0 | Minimum R version for DESCRIPTION | Pipe operator `|>` used throughout codebase; 4.1 is safe floor for HPC clusters |
| devtools | >= 2.4.5 | `load_all()`, `install()`, `check()`, `build()` | The standard dev-cycle tool; wraps roxygen2, rcmdcheck, usethis |
| usethis | >= 2.2.0 | Scaffold DESCRIPTION, NAMESPACE, tests, CI | Automates all R package boilerplate; use `usethis::create_package()` to init |
| roxygen2 | >= 7.3.0 | In-source documentation → man/ + NAMESPACE | `@export`, `@param`, `@examples` in R/ files; NAMESPACE generated automatically |
| testthat | >= 3.2.0 | Unit testing framework | Edition 3 (`edition: 3` in DESCRIPTION) is the current standard; parallel test execution |
| rcmdcheck | >= 1.4.0 | `R CMD check` automation for CI | Used internally by devtools::check(); needed for GitHub Actions |
| pkgdown | >= 2.0.7 | HTML documentation site | Optional but valuable for a research package; generates site from roxygen2 docs |

**Confidence:** HIGH — This toolchain (devtools + usethis + roxygen2 + testthat) has been the unambiguous R package standard since ~2016 and shows no signs of displacement. The Tidyverse/r-lib team maintains all four packages.

**What NOT to use:**
- Do not use `package.skeleton()` (base R) — generates outdated structure, not usethis-compatible
- Do not use `documentation` package — roxygen2 is the standard, not an alternative
- Do not manually edit NAMESPACE — always let roxygen2 manage it via `devtools::document()`

---

### C++ / Rcpp Layer

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Rcpp | >= 1.0.12 | R ↔ C++ interface | Already used in `src/eigen.cpp`; keep as-is |
| RcppEigen | >= 0.3.4.0.0 | Eigen linear algebra in C++ | Already used; `SelfAdjointEigenSolver` for lambda computation; keep as-is |

**Key DESCRIPTION fields required:**

```
LinkingTo: Rcpp, RcppEigen
Imports: Rcpp
SystemRequirements: C++11
```

**src/Makevars** (needed for C++11 on some platforms):

```makefile
CXX_STD = CXX11
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)
```

**What NOT to use:**
- Do not switch to `cpp11` (alternative R/C++ bridge by Davis Vaughan) — no benefit here; Rcpp is already working and RcppEigen has no cpp11 equivalent
- Do not use `armadillo` (RcppArmadillo) — Eigen is already integrated; switching would require rewriting eigen.cpp for no gain

**Confidence:** HIGH — Rcpp + RcppEigen is the correct path. The existing `src/eigen.cpp` uses proper Rcpp attributes (`[[Rcpp::export]]`, `[[Rcpp::depends(RcppEigen)]]`) that slot directly into R package structure without modification.

---

### Cloud Parameter Fetching (Arrow + Parquet)

This is the most critical and nuanced stack decision. The goal: ~15 GB of posterior draws stored as Parquet files in cloud, fetched partially (by species, by draw) without downloading the full dataset.

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| arrow | >= 14.0.0 | Parquet read with predicate pushdown; HTTP filesystem | The only mature R library that supports lazy Parquet reads with column/row filtering without full download |
| curl | >= 5.2.0 | HTTP transport layer used by arrow filesystem | arrow's `S3FileSystem` and `GcsFileSystem` depend on libcurl; verify it ships with the binary |

**Confidence:** HIGH that `arrow` is the correct library. MEDIUM on version — verify on CRAN (*).

#### How Arrow partial reads work in R

Arrow's `open_dataset()` + `filter()` + `collect()` pattern enables predicate pushdown: only row groups matching the filter are fetched over HTTP, not the full file. This works because Parquet stores metadata (row group statistics: min/max per column) at the end of the file, and Arrow performs an HTTP range request for the footer first, then only fetches the row groups that pass the filter.

**The pattern to implement in `R/fetch_params.R`:**

```r
# Option A: Arrow Dataset API (recommended for partitioned datasets)
library(arrow)

ds <- open_dataset(
  "s3://your-bucket/forest-ipm-params/",    # or https:// for public HTTP
  format = "parquet",
  partitioning = c("species", "vital_rate")  # Hive partitioning
)

pars <- ds |>
  filter(species == "28728-ACE-RUB", vital_rate == "growth") |>
  filter(iter <= 100) |>   # partial posterior draws
  collect()                # triggers actual HTTP fetch
```

```r
# Option B: read_parquet() for single files (simpler, less efficient for large files)
pars <- arrow::read_parquet(
  "https://your-host.com/params/growth/28728-ACE-RUB.parquet",
  col_select = c("iter", "r", "Beta", "theta", "tau_temp", "tau_prec", "Lmax", "sigma_obs"),
  as_data_frame = TRUE
)
```

**Option A is strongly recommended** for this project because:
1. Partitioning by `species/vital_rate` means each query only touches relevant files
2. Row group filtering by `iter` enables requesting exactly N posterior draws
3. The ~15 GB distributes across ~31 × 4 = 124 files, each reasonably small

**Critical: Parquet file structure design for efficient partial reads**

The Parquet files must be written with correct row group sizing and sort order to make predicate pushdown effective. Recommended write-time settings (when generating the Parquet files from RDS):

```r
# Writing side (done once, when preparing cloud files)
arrow::write_parquet(
  pars_df,
  "growth/28728-ACE-RUB.parquet",
  chunk_size = 100L   # 100 rows per row group = 100 posterior draws per chunk
)
```

With `chunk_size = 100` and `filter(iter <= 100)`, Arrow fetches exactly 1 row group (one HTTP range request) instead of all 4000 draws.

**Confidence:** MEDIUM — The predicate pushdown mechanism is well-documented in Arrow internals. The `open_dataset()` + HTTP pattern works for public S3 and HTTPS URLs. However, the exact behavior for non-S3 HTTP servers (e.g., Zenodo, HuggingFace Hub) depends on whether the server supports HTTP range requests (`Accept-Ranges: bytes`). AWS S3 public buckets: CONFIRMED to support this. Zenodo: MEDIUM confidence (DOI-based URLs use range requests for large files, but chunking behavior needs verification). HuggingFace Hub: MEDIUM confidence (supports range requests for LFS files).

#### Cloud hosting decision

| Platform | HTTP Range Requests | Cost | Auth Required | Parquet-friendly |
|----------|--------------------|----|---------------|-----------------|
| AWS S3 (public bucket) | Yes (confirmed) | ~$0.023/GB/month + egress | No (public) | Best: native S3FileSystem in arrow |
| Zenodo | Yes (likely) | Free | No | OK: HTTPS URL |
| HuggingFace Hub | Yes (for LFS files) | Free tier | No | OK: HTTPS URL |
| GitHub Releases | No (full download only) | Free | No | Bad: no range requests |
| Google Cloud Storage | Yes (confirmed) | ~$0.020/GB/month | No (public) | Good: GcsFileSystem in arrow |

**Recommendation: AWS S3 public bucket** for production. Rationale: Arrow has native `S3FileSystem` support, HTTP range requests are guaranteed, latency is predictable, and free tiers or research credits are available. Zenodo is acceptable if S3 is unavailable (use HTTPS URL in `open_dataset()`).

**What NOT to use for cloud hosting:**
- GitHub Releases — no range request support; full 15 GB would download per request
- GitHub LFS — has bandwidth limits (1 GB/month free); not viable for a public research package

---

### S3 Class System

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Base R S3 | — | Class definitions for `ipm_kernel`, `ipm_population`, `ipm_params` | S3 is the correct choice for a scientific package: lightweight, no dependencies, composable with existing tidyverse patterns |

**What NOT to use:**
- R5/R6 classes — reference semantics add complexity without benefit here; IPM objects are value-typed (you don't want in-place mutation)
- S4 classes — heavy formalism (slot validators, setClass, setGeneric) that adds nothing for three domain objects; S4 is appropriate for Bioconductor-style multi-class hierarchies
- R7 (new OOP system) — as of 2025, R7 (`R7` package, now called `S7`) is stable but adoption is not yet broad enough to justify for a domain-specific scientific package; use S3 now, migrate later if needed

**S3 design pattern for this project:**

```r
# Constructor
new_ipm_params <- function(growth, mort, rec, sizeIngrowth, species, method) {
  structure(
    list(
      growth = growth,
      mort = mort,
      rec = rec,
      sizeIngrowth = sizeIngrowth,
      species = species,
      method = method
    ),
    class = "ipm_params"
  )
}

# Validator
validate_ipm_params <- function(x) {
  required <- c("growth", "mort", "rec", "sizeIngrowth")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0)
    stop("ipm_params missing fields: ", paste(missing_fields, collapse = ", "))
  invisible(x)
}

# User-facing constructor (exported)
ipm_params <- function(...) validate_ipm_params(new_ipm_params(...))
```

**Confidence:** HIGH — S3 is the correct, idiomatic choice for R scientific packages of this scope.

---

### Documentation

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| roxygen2 | >= 7.3.0 | Inline documentation in R/ files | Standard; generates `man/*.Rd` and `NAMESPACE` |
| pkgdown | >= 2.0.7 | HTML website from roxygen2 docs | Optional but recommended; GitHub Pages deployment via `usethis::use_pkgdown_github_pages()` |

**Mandatory roxygen2 tags for this package:**

```r
#' @param sp Species ID as character (e.g., "28728-ACE-RUB")
#' @param method Either "mean" (posterior mean) or "random" (single posterior draw)
#' @param model Demographic model variant. One of "intcpt", "intcpt_plot",
#'   "intcpt_plot_comp", "intcpt_plot_comp_clim" (default). Can be length-1
#'   (same for all vital rates) or length-3 (growth, mortality, recruitment).
#' @return An \code{ipm_params} object.
#' @export
#' @examples
#' \dontrun{
#'   pars <- getPars_sp("28728-ACE-RUB", method = "mean")
#' }
```

Use `\dontrun{}` for examples that require cloud access. Use `\donttest{}` for examples that are slow but don't need auth.

**Confidence:** HIGH.

---

### Testing

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| testthat | >= 3.2.0 | Unit tests in `tests/testthat/` | Edition 3 standard; `expect_snapshot()` for regression testing complex outputs |
| httptest2 | >= 1.0.0 | Mock HTTP requests in tests | Enables testing the Arrow/Parquet fetching code without network access; records real responses as fixtures |
| withr | >= 2.5.0 | Temporary environment state in tests | Clean setup/teardown of temp files, options, env vars |

**httptest2 is critical:** The parameter fetching functions make HTTP requests. Without mocking, tests would require live cloud access (slow, fragile, not CI-friendly). httptest2 records actual Arrow HTTP requests as JSON fixtures, then replays them in CI.

**Alternative: vcr** — vcr is the other R HTTP mocking library. httptest2 is preferred here because it integrates more cleanly with non-httr HTTP clients like arrow's internal curl calls.

**Confidence:** HIGH for testthat. MEDIUM for httptest2 with arrow — arrow's internal HTTP is handled by libcurl at C++ level, which may not be interceptable by R-level HTTP mocking. Fallback: use `skip_on_ci()` for network-dependent tests and test the data transformation logic separately with pre-saved RDS fixtures.

---

### CI/CD

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| GitHub Actions | — | `R CMD check` on push/PR | Standard for R packages; `usethis::use_github_action("check-standard")` scaffolds it |
| r-lib/actions | v2 | Reusable GitHub Actions for R | Official r-lib workflows handle R version matrix, Rcpp compilation, cache |

**Recommended workflow:** `r-lib/actions/check-r-package@v2` with:
- R versions: release + devel
- OS: ubuntu-latest (fast), macOS-latest (for Rcpp on macOS)
- Skip Windows initially (per PROJECT.md constraints)

**Confidence:** HIGH.

---

### Package Metadata (DESCRIPTION skeleton)

```
Package: forestIPM
Type: Package
Title: Bayesian Hierarchical Integral Projection Models for Tree Population Dynamics
Version: 0.1.0
Authors@R: person("Will", "Vieira", email = "...", role = c("aut", "cre"))
Description: Computes population growth rates (lambda), sensitivities, and
    multi-decade projections for 31 tree species in eastern North America using
    Bayesian hierarchical Integral Projection Models. Parameters are fetched
    on demand from cloud-hosted Parquet files.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.x
Depends:
    R (>= 4.1.0)
Imports:
    Rcpp,
    truncnorm,
    dplyr,
    tidyr,
    purrr,
    stringr,
    arrow (>= 14.0.0)
LinkingTo:
    Rcpp,
    RcppEigen
Suggests:
    testthat (>= 3.0.0),
    httptest2,
    withr,
    devtools,
    pkgdown
Config/testthat/edition: 3
```

---

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| HTTP mocking in tests | httptest2 | vcr | vcr works best with httr/httr2; httptest2 is more HTTP-agnostic |
| Cloud format | Parquet via arrow | HDF5 via rhdf5 | HDF5 has no predicate pushdown; entire file must be read |
| Cloud format | Parquet via arrow | RDS via HTTPS | RDS does not support partial reads; full file downloads required |
| Cloud format | Parquet via arrow | SQLite via RSQLite | SQLite requires full DB download; no chunk-level HTTP range reads |
| OOP system | S3 | S4 | S4 is overkill for 3 classes; S3 is idiomatic for scientific R |
| OOP system | S3 | R6 | R6 reference semantics are wrong for immutable parameter objects |
| C++ bridge | Rcpp + RcppEigen | cpp11 + armadillo | Existing code uses Rcpp/RcppEigen; no reason to rewrite |
| Documentation | roxygen2 | sinew | sinew generates stubs but roxygen2 is the standard |
| Testing | testthat | RUnit | testthat Edition 3 is universally adopted; RUnit is legacy |

---

## Installation

```r
# Core runtime dependencies
install.packages(c(
  "Rcpp",
  "RcppEigen",
  "truncnorm",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "arrow"
))

# Development toolchain
install.packages(c(
  "devtools",
  "usethis",
  "roxygen2",
  "testthat",
  "httptest2",
  "withr",
  "pkgdown",
  "rcmdcheck"
))
```

```bash
# System: macOS (Homebrew) — for arrow binary
brew install apache-arrow

# System: Ubuntu — for arrow binary
sudo apt-get install -y libarrow-dev
```

> Note: `arrow` binary installation is often the hardest step in the stack. The CRAN binary
> on macOS/Linux typically includes the C++ Arrow library. If building from source, see
> https://arrow.apache.org/docs/r/articles/install.html for the correct install flow.
> Add `SystemRequirements: Arrow C++ library >= 14.0` to DESCRIPTION.

---

## Sources

All findings are from training data (knowledge cutoff August 2025). External verification
tools were unavailable during this research session.

**Authoritative sources to verify before implementation:**

- Arrow R package: https://arrow.apache.org/docs/r/ (current version, installation)
- Rcpp DESCRIPTION integration: https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf
- R Packages (2e) by Hadley Wickham: https://r-pkgs.org/ (definitive reference for usethis/devtools/roxygen2/testthat)
- testthat Edition 3: https://testthat.r-lib.org/articles/third-edition.html
- httptest2: https://enpiar.com/httptest2/
- r-lib/actions: https://github.com/r-lib/actions

**Confidence summary:**

| Component | Confidence | Notes |
|-----------|------------|-------|
| devtools/usethis/roxygen2/testthat toolchain | HIGH | Verified stable, unambiguous standard since ~2018 |
| Rcpp + RcppEigen integration | HIGH | Existing code already uses correct attributes |
| S3 class design | HIGH | Idiomatic R, no ambiguity |
| arrow for Parquet partial reads | HIGH (library choice) / MEDIUM (HTTP range behavior) | Arrow is correct choice; HTTP range behavior needs verification per hosting platform |
| httptest2 for mocking arrow HTTP | MEDIUM | R-level mocking may not intercept arrow's libcurl calls; fallback strategy documented |
| Cloud hosting (AWS S3) | MEDIUM | Range request support confirmed; exact arrow S3FileSystem URL format needs verification |
| Arrow + non-S3 HTTPS (Zenodo/HuggingFace) | LOW-MEDIUM | Range requests likely work; behavior with Parquet metadata footer needs testing |
