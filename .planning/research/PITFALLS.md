# Domain Pitfalls

**Domain:** Scientific R package — Bayesian hierarchical IPM with Rcpp/C++, remote Parquet parameter storage
**Researched:** 2026-02-25
**Confidence note:** Web search and WebFetch were unavailable during this research session. Findings are based on (1) direct analysis of this codebase's CONCERNS.md, STACK.md, TESTING.md, ARCHITECTURE.md, and source files, and (2) well-established R package development patterns from the author's training knowledge (R-pkgs.org, CRAN policies, Arrow docs, Rcpp best practices). Confidence levels reflect this.

---

## Critical Pitfalls

Mistakes that cause rewrites, broken installations, or broken reproducibility.

---

### Pitfall 1: NAMESPACE pollution from implicit tidyverse imports

**What goes wrong:** `params.R` uses `map`, `filter`, `pivot_wider`, `pivot_longer`, `reframe`, `select`, `mutate`, `str_replace`, `case_match`, `as_vector`, `set_names`, and the base pipe `|>` — none declared in a DESCRIPTION `Imports:` field. When the code is packaged, `R CMD check` will produce `no visible binding for global variable` and `no visible global function definition` NOTEs for every undeclared symbol. On CRAN these are errors. On GitHub install they are silent until the user's machine doesn't have tidyverse loaded.

**Why it happens:** Scripts that use `source()` inherit the calling session's loaded packages. The package build process strips that inheritance — every symbol must be resolved through the NAMESPACE.

**Consequences:**
- `devtools::install_github()` succeeds but `library(forestIPM)` throws `Error in map(...): could not find function "map"` on any machine without tidyverse pre-loaded
- `R CMD check --as-cran` fails with multiple NOTEs/ERRORs, blocking any CRAN submission
- The tidyverse pipe vs base pipe mismatch (line 79 of params.R uses `%>%` while the rest uses `|>`) will cause a NOTE about mixing pipe styles

**Prevention:**
1. Add all tidyverse packages to `Imports:` in DESCRIPTION (dplyr, tidyr, purrr, stringr)
2. Use `package::function()` qualified calls for all non-base functions, OR
3. Add `@importFrom dplyr filter select mutate` etc. in roxygen2 blocks to populate NAMESPACE correctly
4. Run `devtools::check()` early (Phase 1) and fix every NOTE — NOTEs compound

**Warning signs:** `R CMD check` output contains `no visible binding for global variable` or `no visible global function definition`

**Phase:** Package Structure (Phase 1) — must be resolved before any other work proceeds

---

### Pitfall 2: S3 method export without NAMESPACE registration

**What goes wrong:** If S3 methods for new classes (`ipm_kernel`, `ipm_population`, `ipm_params`) are defined as `print.ipm_kernel` but not registered in NAMESPACE via `S3method(print, ipm_kernel)`, R will find them in the package namespace by accident in interactive use but fail when the package is loaded fresh. The symptom is that `print(my_kernel)` dispatches to the default method instead of the custom one.

**Why it happens:** roxygen2 requires `@export` on S3 methods, but the tag that generates the correct `S3method()` line in NAMESPACE (not `export()`) requires using `@method generic class` OR relying on roxygen2's automatic detection. If the function is named with a dot (e.g., `print.ipm_kernel`) roxygen2 auto-detects it, but if you name it differently and use `@export`, it generates `export(print.ipm_kernel)` which is wrong — it exports the function but does not register the S3 dispatch.

**Consequences:**
- Users must call `print.ipm_kernel(x)` explicitly instead of `print(x)` — destroys the OOP abstraction
- Intermittent failures depending on whether the package namespace is locked or not
- `R CMD check` will warn: `S3method 'print.ipm_kernel' declared in NAMESPACE but not exported`

**Prevention:**
1. Always use `@export` on S3 methods AND ensure the function is named `generic.class` — roxygen2 will then generate `S3method(generic, class)` automatically
2. Run `devtools::document()` and inspect the generated `NAMESPACE` file — every `S3method()` line must be present
3. Create a minimal S3 test: `x <- new_ipm_kernel(...); stopifnot(inherits(x, "ipm_kernel")); print(x)` — run it in a fresh R session immediately after `devtools::load_all()`

**Warning signs:** `methods::existsMethod()` returns FALSE for methods you defined; `print(obj)` shows the default list/matrix method instead of your custom one

**Phase:** S3 Classes (Phase 3 or wherever classes are introduced)

---

### Pitfall 3: Rcpp compilation breaks on `SelfAdjointEigenSolver` for non-symmetric matrices

**What goes wrong:** `src/eigen.cpp` uses `SelfAdjointEigenSolver` which assumes the input matrix is symmetric (self-adjoint). The IPM kernel K = P + F is a non-negative matrix but is generally NOT symmetric (upper-triangular P dominates smaller size classes, F adds recruitment entries in the small-size rows). Calling `getEigenValues()` on a non-symmetric matrix produces incorrect eigenvalues silently — the solver does not throw an error; it returns mathematically meaningless values.

**Why it happens:** `SelfAdjointEigenSolver` is faster than `EigenSolver` and was likely chosen for performance. The existing simulation scripts use `max(Re(eigen(K)$K$values))` from base R (which handles asymmetric matrices correctly), not `getEigenValues()` from the C++ solver — so the bug may never have manifested because the C++ function is not called on the full K matrix in production paths.

**Consequences:**
- Silent numerical errors in lambda calculations if the C++ solver is wired into the public API
- Results are wrong but no error is raised — the worst category of bug for scientific software

**Prevention:**
1. In the package API, either (a) replace `SelfAdjointEigenSolver` with `EigenSolver<MatrixXd>` for general real matrices and take `es.eigenvalues().real()`, or (b) keep the self-adjoint solver only for the symmetric P sub-kernel (growth × survival) and use base R `eigen()` for the full K
2. Add a compile-time or runtime assertion: `if(!M.isApprox(M.transpose())) Rcpp::stop("Matrix must be symmetric")`
3. Document clearly in the function that `getEigenValues()` requires symmetric input

**Warning signs:** Lambda computed via C++ differs from `max(Re(eigen(K)$values))` by more than floating point tolerance; test both on the same kernel matrix

**Phase:** Package Structure / Rcpp integration (Phase 1)

---

### Pitfall 4: Arrow Parquet partial reads require correct chunking at write time

**What goes wrong:** The plan is to use Arrow HTTP range requests to fetch only the parameters needed for one species/draw without downloading the full 15 GB. This only works if the Parquet files were written with row groups that align to query predicates. If all 31 species are written into a single Parquet file with row groups of 100k rows that span multiple species, a `filter(species == "ABBA")` predicate will still read every row group that contains any ABBA row — which could be all of them.

**Why it happens:** Parquet row group sizes are a write-time decision. Arrow's predicate pushdown skips row groups only when the row group's min/max statistics prove the predicate is false for the entire group. Interleaved species data defeats this.

**Consequences:**
- Users experience ~15 GB downloads despite the promise of partial reads
- Latency for every `run_ipm()` call is dominated by network I/O, making the package unusable without fast internet
- The "no local data setup" promise is broken

**Prevention:**
1. Write one Parquet file per species (31 files), or partition the dataset by species using `write_dataset(df, path, partitioning = "species")` — Arrow will create `species=ABBA/part-0.parquet` directories, and reads against `open_dataset()` will only fetch files for the queried species
2. Within each species file, sort rows by posterior draw index and set row group size to one full draw (all parameters for draw i) so a `filter(draw == 42)` fetches exactly one row group
3. Test actual byte transfer: use `curl -v` with a `Range:` header against your cloud host before committing to a schema

**Warning signs:** `arrow::read_parquet()` with a filter takes as long as without a filter; `httr::GET()` with a range header returns the full file size

**Phase:** Cloud Parameter Storage (Phase 2)

---

### Pitfall 5: Cloud host does not support HTTP range requests

**What goes wrong:** Arrow's S3/GCS/HTTP filesystem backends rely on HTTP 206 Partial Content responses (byte-range requests) to fetch individual row groups. Several academic hosting services (Zenodo, institutional repositories, some S3 presigned URL configurations) disable or do not support range requests. If the server returns `HTTP 200` for all range requests and sends the full file, Arrow will silently download the entire file.

**Why it happens:** Range request support is server-side configuration, not client-side. A server that doesn't implement RFC 7233 simply ignores the `Range:` header.

**Consequences:**
- The entire 15 GB is downloaded on first parameter fetch — the package hangs for minutes on a standard connection
- Users on slow connections or HPC login nodes (where large downloads are prohibited) cannot use the package

**Prevention:**
1. Before committing to a cloud host, verify: `curl -I --range 0-100 <url>` — a correctly configured server returns `HTTP/1.1 206 Partial Content` with `Accept-Ranges: bytes`
2. AWS S3 public buckets reliably support range requests — HIGH confidence
3. HuggingFace Hub dataset files also support range requests via their CDN
4. Zenodo: range requests are supported on individual files but the URL must point directly to the file (not the landing page)
5. Add a startup check in the package: `options(arrow.use_threads = TRUE)` and validate the first request returns a partial response before proceeding

**Warning signs:** First `getPars_sp()` call takes > 30 seconds; `curl -v` shows `HTTP/1.1 200 OK` instead of `206 Partial Content`

**Phase:** Cloud Parameter Storage (Phase 2) — verify before writing any conversion code

---

### Pitfall 6: testthat tests that call `getPars_sp()` will fail on CI and for all users

**What goes wrong:** The current `getPars_sp()` reads from `data/output_sim_processed/*.RDS`. Any test that calls this function directly will fail on GitHub Actions, on any collaborator's machine, and for any user who installs the package — because the 15 GB data is not in the repository. This is the most common cause of failed `R CMD check` on CI for data-heavy scientific packages.

**Why it happens:** Scripts that `source()` R files run in a known working directory with data adjacent. Tests run from the package's `tests/testthat/` directory with no assumption about external data.

**Consequences:**
- `devtools::check()` fails with "cannot open file" errors
- CI pipeline (GitHub Actions) always fails — package appears broken
- Users who run `devtools::test()` after installing from GitHub see immediate test failures

**Prevention:**
1. Create `tests/testthat/fixtures/` with minimal synthetic parameter sets (small named lists matching the structure returned by `getPars_sp()`) — use these for all unit tests of kernel and vital rate functions
2. Gate integration tests behind `testthat::skip_if_not(has_remote_data())` where `has_remote_data()` checks whether Arrow can reach the cloud endpoint
3. Use `withr::with_tempdir()` and `mockery::mock()` to stub file reads in unit tests
4. The `getPars_sp()` function itself should be tested only for its parsing logic using a small in-memory RDS fixture, not for its ability to locate real files

**Warning signs:** `devtools::check()` shows test failures with "No such file or directory"; test suite passes locally but fails on CI

**Phase:** Testing (any phase) — establish fixtures in Phase 1 before any test is written

---

## Moderate Pitfalls

---

### Pitfall 7: `@export` on internal kernel functions inflates the public API

**What goes wrong:** `R/kernel.R` already has `#' @export` on `mkKernel` and `init_pop`. During conversion, the instinct is to export every function that was previously `source()`-accessible. But `vonBertalanffy_lk`, `P_xEC`, `ingrowth_lk`, and the `survival_f`/`vonBertalanffy_f`/`ingrowth_f` trio are implementation details of the kernel — not part of the user-facing API. Exporting them creates a public API surface that must be maintained forever (CRAN compatibility rules discourage unexported → exported changes in patch releases, but exported → unexported is a breaking change).

**Prevention:**
1. Define the public API explicitly: only `run_ipm()`, `mkKernel()`, `init_pop()`, `getPars_sp()`, plus S3 methods for the new classes
2. All probability density helpers (`vonBertalanffy_lk`, `P_xEC`, `ingrowth_lk`) should be internal (no `@export`), documented with `@noRd`
3. Use `@keywords internal` if you want roxygen2 to skip them in the documentation index

**Warning signs:** Package has > 15 exported symbols; users call `package:::vonBertalanffy_lk` in their scripts (fragile by definition)

**Phase:** Package Structure (Phase 1)

---

### Pitfall 8: `pars_to_list()` uses `%>%` pipe but `getPars_sp()` returns with `|>` — both are undeclared

**What goes wrong:** Line 79 of `params.R` uses the magrittr pipe `%>%` with the comment "fuck tidyverse for not wanting to implement a named list argument." The rest of the codebase uses the base pipe `|>`. In a package context, `%>%` requires either `@importFrom magrittr %>%` or `Import: magrittr` in DESCRIPTION. The base pipe `|>` requires R >= 4.1.0. If the package only declares `Depends: R (>= 3.5.0)` out of habit, `R CMD check` will catch the base pipe on the declared minimum version.

**Prevention:**
1. Eliminate `%>%` — rewrite `pars_to_list()` using base R or `|>` with a named-list workaround (`setNames()`)
2. Declare `Depends: R (>= 4.1.0)` if using `|>` throughout
3. Run `R CMD check` with `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE` to catch version-related issues locally

**Warning signs:** `R CMD check` NOTE: `%>%` no visible global function definition; runtime error on R 4.0.x

**Phase:** Package Structure (Phase 1)

---

### Pitfall 9: Global state in `init_pop()` breaks test isolation

**What goes wrong:** Line 195 of `kernel.R` checks `if(exists('fct'))`. This `exists()` call looks in the parent environment — in an interactive session, `fct` might persist from a previous call. In a package context with proper namespacing, the function has its own environment and `fct` will always be `FALSE` on first call unless assigned within the current invocation. But the logic is fragile: if the accuracy loop never executes (because `diff_N <= accuracy` on the first check), `fct` is never assigned, and the function falls through to `N_out <- dbh_den` (correct) — but this is undocumented and relies on `exists()` returning FALSE. A refactor that changes environment scoping (e.g., wrapping in `local()`) could break it silently.

**Prevention:**
1. Replace the `exists('fct')` pattern with explicit initialization: `fct <- NULL` before the loop, and check `if(!is.null(fct))` afterward
2. This is a one-line fix that makes the logic explicit and testable

**Warning signs:** `init_pop()` returns wrong results after being called twice in the same R session (stale `fct` in global env)

**Phase:** Package Structure (Phase 1) — fix before writing any tests

---

### Pitfall 10: `as_vector()` from purrr/vctrs is deprecated and removed

**What goes wrong:** `params.R` lines 44 and 54 call `as_vector()` which was deprecated in purrr 1.0.0 (released 2022) and scheduled for removal. The replacement is `purrr::list_simplify()` or `unlist()` depending on context. On purrr >= 1.0.2, `as_vector()` throws a deprecation warning; in future purrr versions it may be removed entirely.

**Consequences:** Package emits deprecation warnings on every parameter load; after purrr deprecates the function completely, the package breaks silently on new purrr versions.

**Prevention:**
1. Replace `as_vector()` with `unlist(use.names = TRUE)` for named numeric vectors (which is all that's needed here)
2. Test that the replacement produces identical structure: `identical(names(unlist(x)), names(as_vector(x)))`

**Warning signs:** `Warning: as_vector() was deprecated in purrr 1.0.0` during `getPars_sp()` calls

**Phase:** Package Structure (Phase 1)

---

### Pitfall 11: `truncnorm` as a hard dependency for a single function call

**What goes wrong:** `truncnorm::dtruncnorm()` is called exactly once in `kernel.R` line 58-64 for the ingrowth size distribution. This creates a hard package dependency for a function that can be replaced with two lines of base R. If `truncnorm` is not installed, `devtools::install_github()` fails with a missing dependency error. If `truncnorm` ever drops from CRAN, the package breaks.

**Prevention:**
1. Replace `truncnorm::dtruncnorm(x, a, b, mean, sd)` with the base R equivalent:
   ```r
   dnorm(x, mean, sd) / (pnorm(b, mean, sd) - pnorm(a, mean, sd))
   ```
2. This eliminates the dependency entirely — CONCERNS.md already recommends this migration

**Warning signs:** `R CMD check` NOTE: `truncnorm` is in `Imports:` but only used once; user install fails because truncnorm is unavailable on their CRAN mirror

**Phase:** Package Structure (Phase 1)

---

### Pitfall 12: Hardcoded path separators and relative paths in `getPars_sp()`

**What goes wrong:** `getPars_sp()` constructs paths with `paste0(path, '/', ...)` — a forward-slash hardcode that fails on Windows. The default argument `path = file.path('data', 'output_sim_processed')` is relative to the working directory, which is undefined in a package context (it depends on where the user calls the function from). In the post-migration world where parameters come from Arrow/Parquet over HTTP, this entire path construction must be replaced — but if it is left as-is during Phase 1, it will silently fail for Windows users and any user who doesn't set their working directory correctly.

**Prevention:**
1. Replace `paste0(path, '/', ...)` with `file.path(path, ...)` throughout
2. Add a validation step: `if (!file.exists(pars_dir[1])) stop("Parameter file not found: ", pars_dir[1], ". Did you set path= correctly?")`
3. When migrating to Arrow, deprecate the `path` argument and add an `url` argument that defaults to the cloud endpoint

**Warning signs:** Package works on macOS/Linux but fails on Windows with "cannot open file"; tests fail with confusing path errors

**Phase:** Package Structure (Phase 1) for path separator fix; Cloud Storage (Phase 2) for the full migration

---

## Minor Pitfalls

---

### Pitfall 13: `R CMD check` NOTE from non-ASCII characters

**What goes wrong:** The comment "fuck tidyverse..." on line 80 of `params.R` won't cause a NOTE (it's ASCII), but any non-ASCII characters in comments, documentation strings, or strings (e.g., accented characters in author names in DESCRIPTION, or `µ`, `°` symbols in documentation examples) trigger `R CMD check` NOTEs about "non-ASCII characters in package."

**Prevention:** Use `\u00B0` unicode escapes for special characters in strings; keep comments pure ASCII; check DESCRIPTION author names use ASCII-safe encoding

**Phase:** Package Structure (Phase 1)

---

### Pitfall 14: `matrix.image()` name conflicts with base R

**What goes wrong:** The visualization function is named `matrix.image` with a dot. R's dispatch rules mean this looks like an S3 method `image.matrix` spelled backwards — but it is not registered as one. The name `matrix.image` could conflict if `matrix` is ever used as an S3 generic (unlikely but possible in future R versions). More practically, `image.matrix` (the correct S3 name) might be what users expect after reading `?image`.

**Prevention:** Rename to `plot_kernel()` or `image_kernel()` — snake_case or explicit names that don't look like S3 dispatch

**Phase:** Package Structure (Phase 1) — rename before writing any documentation

---

### Pitfall 15: Eigenvalue discrepancy between `getEigenValues` (C++) and base R `eigen()`

**What goes wrong:** The existing simulation scripts compute lambda as `max(Re(eigen(K)$values))` using base R. The C++ `getEigenValues()` returns all eigenvalues. The dominant eigenvalue (lambda) is the maximum real part. If the package API exposes `getEigenValues()` directly and users take `max()` of the full returned vector (which contains complex eigenvalues with small imaginary parts), they may get a slightly different answer than `max(Re(...))` depending on floating point.

**Prevention:**
1. Expose a high-level `compute_lambda(K)` function that wraps the C++ call and applies `max(Re(...))` internally
2. Never expose `getEigenValues()` as a public API function — it's an implementation detail

**Phase:** S3 Classes / high-level API phase

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Package structure setup | NAMESPACE pollution from tidyverse imports | Run `devtools::check()` in first commit; zero NOTEs policy |
| Package structure setup | `exists('fct')` global state bug | Fix before any test is written — tests will give false positives otherwise |
| Package structure setup | `as_vector()` deprecation | Replace with `unlist()` immediately |
| Package structure setup | `%>%` vs `\|>` pipe mixing | Standardize to `\|>` and declare `R >= 4.1.0` |
| Package structure setup | `matrix.image` naming | Rename before documenting |
| Rcpp integration | `SelfAdjointEigenSolver` on asymmetric matrix | Validate solver choice before exposing C++ function in API |
| Cloud parameter storage | Parquet row group alignment | Decide on partitioning strategy before writing first Parquet file |
| Cloud parameter storage | Host range request support | Verify with `curl` before writing any Arrow fetch code |
| S3 class definitions | Method not dispatching | Inspect NAMESPACE after `devtools::document()` for every new method |
| Testing setup | Tests calling real parameter files | Create synthetic fixtures in Phase 1 before writing any test |
| Testing setup | Non-deterministic `init_pop()` in tests | Use `set.seed()` in test setup or refactor to accept seed parameter |

---

## Sources

**Codebase-derived findings (HIGH confidence):**
- Direct analysis of `/Users/wvieira/GitHub/forest-IPM/R/params.R` — tidyverse imports, `%>%` pipe, `as_vector()` usage confirmed
- Direct analysis of `/Users/wvieira/GitHub/forest-IPM/src/eigen.cpp` — `SelfAdjointEigenSolver` confirmed, asymmetry issue identified
- Direct analysis of `/Users/wvieira/GitHub/forest-IPM/R/kernel.R` — `exists('fct')` pattern, `truncnorm::dtruncnorm()` call confirmed
- `/Users/wvieira/GitHub/forest-IPM/.planning/codebase/CONCERNS.md` — tidyverse undeclared dependency, `exists()` global state, truncnorm migration, path construction issues all previously flagged

**Well-established R package patterns (MEDIUM confidence — standard practices, no web verification performed):**
- NAMESPACE S3method registration behavior: documented in Writing R Extensions §1.5.1
- roxygen2 `@export` vs `S3method()` generation: documented in roxygen2 package vignette
- Arrow predicate pushdown requires row group alignment: Arrow columnar format specification
- HTTP 206 Partial Content required for range requests: RFC 7233
- testthat external data pattern: documented in testthat vignette "Testing with external data"

**Flagged for verification when web access is available:**
- `purrr::as_vector()` removal timeline — verify against purrr changelog at https://purrr.tidyverse.org/news/index.html
- Arrow `open_dataset()` with HTTP partitioning — verify against https://arrow.apache.org/docs/r/articles/dataset.html
- HuggingFace Hub range request support — verify with `curl -I --range 0-100 <test-url>`
