# Architecture Patterns

**Domain:** Scientific R package with remote Bayesian posteriors + Rcpp acceleration
**Researched:** 2026-02-25
**Confidence:** HIGH (R package structure, Arrow/Parquet, Rcpp patterns are well-documented and stable)

---

## Recommended Architecture

The package wraps an existing working simulation engine. The architectural challenge is not the IPM math — that already works — but adding three production-grade layers on top:

1. A **remote data layer** that fetches Parquet-stored posteriors on demand
2. A **domain object layer** (S3 classes) that makes IPM objects first-class R citizens
3. A **package infrastructure layer** (DESCRIPTION, NAMESPACE, C++ compilation) that makes it installable

```
┌─────────────────────────────────────────────────────────────────┐
│                     User-Facing API                             │
│  run_ipm(species, coords, climate)                              │
│  project_population(kernel, n_steps)                            │
└────────────────────┬────────────────────────────────────────────┘
                     │
┌────────────────────▼────────────────────────────────────────────┐
│                  Domain Object Layer (S3)                       │
│  ipm_params   ipm_kernel   ipm_population   ipm_projection      │
└───────┬───────────────┬───────────────┬────────────────────────┘
        │               │               │
┌───────▼──────┐ ┌──────▼──────┐ ┌─────▼──────────────────────┐
│ Remote Data  │ │ IPM Engine  │ │  C++ Acceleration          │
│ (Arrow +     │ │ (R kernel,  │ │  (RcppEigen eigenvalues)   │
│  Parquet)    │ │  vital      │ │                            │
│              │ │  rates,     │ │                            │
│              │ │  competition│ │                            │
└───────────── ┘ └─────────────┘ └────────────────────────────┘
```

---

## Component Boundaries

| Component | Responsibility | Communicates With | Files |
|-----------|---------------|-------------------|-------|
| **Remote Data** (`R/params.R` extended) | Fetch species posterior draws from cloud Parquet; cache locally; validate | User API → returns `ipm_params` | `R/params.R`, `R/cache.R` |
| **Vital Rates** (`R/vital_rates.R`) | Von Bertalanffy growth, logistic survival, ingrowth; pure math functions | Kernel assembly (consumes parameters) | `R/vital_rates.R` |
| **Competition** (`R/BasalArea_competition.R`) | Basal area aggregation from size distributions | Kernel assembly (provides BA_comp vectors) | `R/BasalArea_competition.R` |
| **Kernel Assembly** (`R/kernel.R`) | Build P and F matrices via numerical integration; assemble full K | Vital rates + Competition → `ipm_kernel` | `R/kernel.R` |
| **C++ Eigenvalue Solver** (`src/eigen.cpp`) | Fast eigenvalue computation for lambda | Kernel assembly (consumes K matrix) | `src/eigen.cpp` |
| **S3 Domain Objects** (`R/classes.R`) | Constructor, validator, print/summary methods for all domain types | Used by all layers | `R/classes.R` |
| **User API** (`R/ipm.R`) | `run_ipm()`, `project_population()`, `lambda()` — high-level orchestration | All components | `R/ipm.R` |
| **Visualization** (`R/plot.R`) | S3 `plot()` methods for kernels, projections, size distributions | S3 objects | `R/plot.R` |

---

## Data Flow

### Primary Flow: `run_ipm(species, coords, climate)`

```
1. Parameter fetch
   get_params(species, method = "random", model = "intcpt_plot_comp_clim")
     → arrow::open_dataset(remote_url) |> filter(species_id == sp) |> collect()
     → Returns: ipm_params object (nested named list, now with class attr)

2. Population initialization
   init_population(params, n_size_classes = 100)
     → Returns: ipm_population (meshpts + Nvec + h, with class attr)

3. Per-timestep kernel assembly [loop]
   a. size_to_BAcomp(population_intra, population_inter, plot_size)
        → BA_comp_intra vector, BA_comp_inter vector

   b. build_kernel(params, BA_comp, delta_time, climate)
        → P matrix (growth × survival via P_xEC outer product)
        → F matrix (ingrowth via ingrowth_lk outer product)
        → K = P + F
        → Returns: ipm_kernel (list with K, P, F matrices + class attr)

   c. lambda(kernel)
        → getEigenValues(K) [C++ via Rcpp]
        → Returns scalar, stored in projection record

4. Population projection
   N(t+1) = K %*% N(t)$Nvec * h
     → Returns: ipm_projection (data frame: year, lambda, total_N, BA)

5. Return assembled result
   → list(kernel = ipm_kernel, population = ipm_population,
           projection = ipm_projection, params = ipm_params)
```

### Remote Data Flow: Parquet Fetching

```
get_params(species = "28731-ACE-RUB", method = "random")
  │
  ├─ Check local cache (~/.cache/forestIPM/ or tools::R_user_dir())
  │     HIT  → read_parquet(cache_path) → deserialize → ipm_params
  │     MISS ↓
  │
  ├─ arrow::open_dataset(remote_parquet_url)   ← HTTP range request
  │     |> filter(species_id == species)        ← predicate pushdown
  │     |> filter(vital_rate %in% c("growth", "mort", "rec"))
  │     |> collect()                            ← materialize only matched rows
  │
  ├─ write_parquet(result, cache_path)         ← cache for reuse
  │
  └─ new_ipm_params(result)                    ← wrap in S3 class
```

Arrow's `open_dataset()` with an HTTP URL uses HTTP range requests natively when the server supports them (S3, HuggingFace Hub, most CDNs all do). The predicate pushdown (`filter()` before `collect()`) is the critical optimization: only rows matching the species are transferred, not the full 15 GB file.

---

## S3 Class Design

### Design Principle: Thin Wrappers, Not Heavy OOP

R's S3 system is informal by design. Use it to add `class` attributes to existing list/data structures. This enables `print()`, `summary()`, `plot()`, `[` dispatch without rewriting the math core. Do NOT use R5/Reference Classes or R6 — they are harder to serialize and unnecessary here.

### `ipm_params` — Parameter Container

```r
# Constructor
new_ipm_params <- function(
  species_id,
  model,
  method,     # "mean" or "random"
  draw_id = NULL,  # which posterior draw (if method = "random")
  growth,     # named numeric vector
  mort,       # named numeric vector
  rec,        # named numeric vector
  sizeIngrowth  # named numeric vector
) {
  structure(
    list(
      species_id = species_id,
      model = model,
      method = method,
      draw_id = draw_id,
      growth = growth,
      mort = mort,
      rec = rec,
      sizeIngrowth = sizeIngrowth
    ),
    class = "ipm_params"
  )
}

# Validator: called inside constructor
validate_ipm_params <- function(x) {
  required_growth <- c("Lmax", "r", "Beta", "tau_temp", "tau_prec",
                        "optimal_temp", "optimal_prec", "sigma_obs")
  stopifnot(all(required_growth %in% names(x$growth)))
  invisible(x)
}

# Methods
print.ipm_params <- function(x, ...) { ... }
summary.ipm_params <- function(object, ...) { ... }
```

### `ipm_kernel` — Kernel Matrix Container

```r
new_ipm_kernel <- function(K, P, F_mat, meshpts, h, params_ref) {
  structure(
    list(
      K = K,        # full kernel (n × n matrix)
      P = P,        # growth × survival subkernel
      F = F_mat,    # recruitment subkernel
      meshpts = meshpts,
      h = h,
      n = nrow(K),
      params_ref = params_ref  # ipm_params used to build this kernel
    ),
    class = "ipm_kernel"
  )
}

# Methods
print.ipm_kernel <- function(x, ...) { ... }
plot.ipm_kernel <- function(x, ...) { ... }  # delegates to matrix_image
lambda.ipm_kernel <- function(x, ...) { max(getEigenValues(x$K)) }
```

### `ipm_population` — Size Distribution State

```r
new_ipm_population <- function(meshpts, Nvec, h, t = 0) {
  structure(
    list(
      meshpts = meshpts,  # size class midpoints (mm)
      Nvec = Nvec,        # density per class
      h = h,              # bin width
      t = t               # time step (0 = initial)
    ),
    class = "ipm_population"
  )
}

# Methods
print.ipm_population <- function(x, ...) { ... }
plot.ipm_population <- function(x, ...) { ... }   # size distribution histogram
total_N.ipm_population <- function(x, ...) { sum(x$Nvec) }
```

### `ipm_projection` — Time Series Output

```r
# This is a data.frame with class prepended
new_ipm_projection <- function(df) {
  # df has columns: year, lambda, total_N, BA_intra, BA_inter, ...
  class(df) <- c("ipm_projection", "data.frame")
  df
}

# Methods
plot.ipm_projection <- function(x, ...) { ... }  # ggplot2 time series
summary.ipm_projection <- function(object, ...) { ... }
```

**Why S3 over S4/R5/R6:**
- S3 is idiomatic R. All base R uses it. Bioconductor uses S4 but only because genomic data has strict slot contracts — IPM objects don't need that formalism.
- S3 objects serialize trivially via `saveRDS()` — critical for caching.
- `arrow::Table` objects are NOT directly serializable; always `collect()` to a plain data.frame before wrapping in S3.

---

## Patterns to Follow

### Pattern 1: Arrow Predicate Pushdown for Remote Parquet

**What:** Filter before materializing. Let Arrow push column/row predicates to the Parquet file server via HTTP range requests, returning only the data needed.

**When:** Any time you read from a remote Parquet source. Never `collect()` first then filter.

```r
get_params <- function(species, model, method, n_draws = 1) {
  ds <- arrow::open_dataset(
    sources = remote_parquet_url(species, model),
    format = "parquet"
  )

  # Predicates evaluated server-side — only matching rows transferred
  result <- ds |>
    dplyr::filter(vital_rate %in% c("growth", "mort", "rec", "sizeIngrowth")) |>
    dplyr::filter(if (method == "random") iter %in% sample(1:4000, n_draws) else TRUE) |>
    dplyr::collect()

  new_ipm_params(result)
}
```

**Confidence:** HIGH. Arrow's `open_dataset()` supports HTTP URLs and uses HTTP range requests when the server responds with `Accept-Ranges: bytes`. AWS S3, HuggingFace Hub datasets, and Zenodo all support this.

### Pattern 2: XDG-Compliant Local Cache

**What:** Store fetched parameters in a user-local cache directory. Use `tools::R_user_dir()` (R >= 4.0) for the platform-correct cache location.

**When:** Every remote fetch. Write Parquet to cache after first fetch; subsequent calls hit cache.

```r
cache_dir <- function() {
  dir <- tools::R_user_dir("forestIPM", which = "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

cache_path <- function(species, model) {
  file.path(cache_dir(), paste0(species, "_", model, ".parquet"))
}

# Usage pattern
get_params_cached <- function(species, model, method) {
  path <- cache_path(species, model)
  if (file.exists(path)) {
    result <- arrow::read_parquet(path)
  } else {
    result <- fetch_remote_params(species, model)
    arrow::write_parquet(result, path)
  }
  new_ipm_params(result, method = method)
}
```

**Confidence:** HIGH. `tools::R_user_dir()` is the recommended approach per CRAN policy (avoids writing to package install dir or user home).

### Pattern 3: Rcpp in Package — Compilation at Install Time

**What:** Place all C++ in `src/`. Add `LinkingTo: Rcpp, RcppEigen` in DESCRIPTION. Call `Rcpp::compileAttributes()` to auto-generate `RcppExports.R` and `RcppExports.cpp`. The package build system (`R CMD INSTALL`) compiles the C++ automatically.

**When:** Any package with C++ code. This is the standard path.

```
Package structure:
  DESCRIPTION          ← LinkingTo: Rcpp, RcppEigen
  NAMESPACE            ← useDynLib(forestIPM, .registration=TRUE); importFrom(Rcpp, evalCpp)
  src/
    eigen.cpp          ← existing code (already correct [[Rcpp::export]] annotation)
    Makevars           ← optional: add compiler flags (e.g., -O3)
  R/
    RcppExports.R      ← auto-generated by compileAttributes(); DO NOT edit manually
```

The existing `src/eigen.cpp` already has `// [[Rcpp::export]]` and `// [[Rcpp::depends(RcppEigen)]]` — it will compile cleanly with no changes. Just run `Rcpp::compileAttributes()` after creating DESCRIPTION.

**Confidence:** HIGH. This is fully documented in "Writing R Extensions" and the Rcpp vignettes.

### Pattern 4: usethis-Driven Package Scaffold

**What:** Use `usethis::create_package()`, `usethis::use_rcpp_eigen()`, `usethis::use_testthat()`, `usethis::use_roxygen_md()` to generate all infrastructure files correctly. Avoids manual DESCRIPTION/NAMESPACE editing errors.

**When:** Package initialization phase only.

```r
# Run once to scaffold:
usethis::create_package("forestIPM")
usethis::use_rcpp_eigen()     # adds Imports: Rcpp, LinkingTo: Rcpp, RcppEigen
usethis::use_testthat(3)
usethis::use_mit_license()
usethis::use_github_actions()  # CI/CD with C++ compilation check
```

**Confidence:** HIGH.

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: `source()` Inside Package Functions

**What:** Calling `source("R/vital_rates.R")` from inside another package function.

**Why bad:** Package functions are loaded via the namespace, not `source()`. This causes "file not found" errors when the package is installed. The existing simulation scripts use this pattern — it must NOT be carried into the package.

**Instead:** Proper NAMESPACE with `useDynLib` and roxygen2 `@export` tags. All R files in `R/` are automatically loaded when the package loads.

### Anti-Pattern 2: Hardcoded Paths in Package Functions

**What:** `path = file.path('data', 'output_sim_processed')` as a default argument.

**Why bad:** Relative paths break when package is installed. The existing `getPars_sp()` does this — it works in the script workflow but will fail in a package.

**Instead:** Default to the cache dir (`tools::R_user_dir("forestIPM", "cache")`), or require the caller to provide a URL/path explicitly. Remote Parquet URLs become the new default.

### Anti-Pattern 3: `collect()` Before Filter on Remote Data

**What:** `arrow::open_dataset(url) |> collect() |> filter(species == "...")`

**Why bad:** Downloads the full 15 GB before filtering. Defeats the entire point of Parquet + Arrow.

**Instead:** Always filter, select, and limit BEFORE `collect()`. Arrow evaluates these as pushdown predicates at the file level.

### Anti-Pattern 4: Monolithic Parameter Files Per Vital Rate

**What:** Storing one giant Parquet file per vital rate (e.g., `growth_all_species.parquet`) — all 31 species in one file.

**Why bad:** Even with predicate pushdown, reading a large file to extract one species wastes bandwidth and adds latency. Parquet row groups are the unit of skipping; if species are randomly distributed across row groups, you still scan much of the file.

**Instead:** Partition Parquet by species_id (`hive partitioning`): `s3://bucket/params/species_id=28731-ACE-RUB/growth.parquet`. Arrow automatically uses partition columns as free filters — zero-overhead species selection.

```
Remote layout (recommended):
  params/
    species_id=28731-ACE-RUB/
      growth.parquet   (4k draws × ~10 parameters)
      mort.parquet
      rec.parquet
      sizeIngrowth.parquet
    species_id=28733-ACE-SAC/
      ...
```

### Anti-Pattern 5: S4 Classes for IPM Domain Objects

**What:** Using `setClass()`, `setGeneric()`, `setMethod()` for `ipm_kernel` etc.

**Why bad:** S4 adds dispatch overhead, complicates serialization (S4 objects need special handling with `saveRDS`), and is overkill for these domain objects which have stable, flat structure. S4 is justified for Bioconductor's genomic data (strict slot contracts, multiple inheritance) but not for a numerical list-of-matrices.

**Instead:** S3 with explicit constructors and validators (see "S3 Class Design" section).

---

## Scalability Considerations

| Concern | Current (local RDS) | Package (remote Parquet) | Future (ML decoder) |
|---------|--------------------|--------------------------|--------------------|
| Parameter access | Read full RDS per species | HTTP range → species rows only | Predict from lat/lon, no download |
| Memory | Load full posterior | Load n_draws rows on demand | Negligible (inference output only) |
| First-run latency | Fast (local) | 2-10s per species (HTTP + cache) | <1s (local inference) |
| Subsequent runs | Fast | Instant (cache hit) | Instant |
| New species | Requires refit | Not supported (31 fixed species) | Supported |
| 31 species × 4k draws | ~15 GB total | Never all in memory | Not applicable |

---

## Build Order (Dependencies)

This order matters because each layer depends on the previous being functional.

### Phase 1: Package Skeleton (no dependencies)

```
1. DESCRIPTION file (name, version, license, dependencies)
2. NAMESPACE (initially empty, roxygen2 will manage it)
3. Move R/*.R into package R/ directory
4. Run Rcpp::compileAttributes() → generates RcppExports.R + RcppExports.cpp
5. Verify: devtools::load_all() works; getEigenValues() callable
```

Nothing else can be built until `devtools::load_all()` succeeds. The C++ must compile first because other R code imports `getEigenValues`.

### Phase 2: S3 Classes (depends on: Phase 1)

```
6. new_ipm_params() constructor + validator
7. new_ipm_population() constructor + validator
8. new_ipm_kernel() constructor + validator
9. new_ipm_projection() constructor
10. print/summary methods for each class
11. Unit tests for all constructors
```

S3 classes must exist before the data layer (Phase 3) and high-level API (Phase 4) can return properly typed objects.

### Phase 3: Remote Data Layer (depends on: Phase 2)

```
12. Parquet file generation script (convert existing RDS → Parquet, partitioned by species)
13. Cloud hosting setup (S3 bucket or HuggingFace dataset)
14. get_params() with Arrow open_dataset() + predicate pushdown
15. Local caching with tools::R_user_dir()
16. Integration tests for parameter fetch (can run against local Parquet first)
```

The Parquet layout decision (partitioning scheme) must be finalized before cloud upload — repartitioning after upload is expensive.

### Phase 4: High-Level API (depends on: Phases 1-3)

```
17. run_ipm() orchestration function
18. project_population() simulation loop
19. lambda() generic + ipm_kernel method
20. plot() methods for all S3 classes
21. roxygen2 documentation for all exported functions
22. pkgdown site (optional)
```

### Phase 5: Testing and CI (should run throughout)

```
23. testthat tests for each component as built
24. GitHub Actions workflow (R-CMD-check on ubuntu-latest, macos-latest)
25. C++ compilation check in CI (catches linker issues early)
```

---

## Component Dependency Graph

```
DESCRIPTION/NAMESPACE
       │
       ├──► src/eigen.cpp (compiled at install)
       │         │
       │    getEigenValues()
       │         │
R/vital_rates.R  │   R/BasalArea_competition.R
       │         │              │
       └────►R/kernel.R ◄───────┘
                 │
            (returns matrices)
                 │
           R/classes.R  ←── R/params.R ←── Arrow/Parquet (remote)
                 │
            R/ipm.R (run_ipm, project_population)
                 │
            R/plot.R (plot.ipm_*)
```

Data flows top-to-bottom. Infrastructure (DESCRIPTION, src/) is the foundation. Classes define the contracts. The API orchestrates everything.

---

## Technology Decisions

| Decision | Choice | Rationale | Confidence |
|----------|--------|-----------|------------|
| Remote data format | Parquet + Arrow | Columnar, predicate pushdown, HTTP range requests, R support via `arrow` package | HIGH |
| Parquet partitioning | Hive by `species_id` | Free filter at directory level, no row group scanning | HIGH |
| Local cache | `tools::R_user_dir()` | Platform-correct, CRAN-compliant, persistent across sessions | HIGH |
| S3 class system | S3 (not S4/R6) | Idiomatic R, trivial serialization, no formalism needed | HIGH |
| C++ integration | Rcpp + RcppEigen (existing) | Already working, standard package support, fast eigenvalues | HIGH |
| Package scaffold | `usethis` helpers | Generates correct DESCRIPTION/NAMESPACE/CI files | HIGH |
| Cloud host | TBD (S3/HuggingFace) | Both support HTTP range requests; HuggingFace has free tier for research | MEDIUM |

---

## Sources

- R Packages (2e) by Hadley Wickham: https://r-pkgs.org — package structure, NAMESPACE, roxygen2
- Writing R Extensions (CRAN official): https://cran.r-project.org/doc/manuals/r-release/R-exts.html — Rcpp, useDynLib, DESCRIPTION
- Apache Arrow R documentation: https://arrow.apache.org/docs/r/ — open_dataset, HTTP filesystems, predicate pushdown
- Rcpp vignettes (package vignette, not web): `vignette("Rcpp-package")` — compileAttributes, src/ structure
- Existing codebase analysis: `/Users/wvieira/GitHub/forest-IPM/R/` and `/Users/wvieira/GitHub/forest-IPM/src/`
