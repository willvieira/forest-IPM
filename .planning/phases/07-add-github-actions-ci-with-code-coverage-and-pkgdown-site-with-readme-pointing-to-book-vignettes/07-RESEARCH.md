# Phase 7: Add GitHub Actions CI with code coverage and pkgdown site with README pointing to book vignettes - Research

**Researched:** 2026-03-13
**Domain:** GitHub Actions for R packages (r-lib/actions), pkgdown site generation, codecov.io integration
**Confidence:** HIGH

---

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

**CI matrix:**
- Linux only — Ubuntu runner, R release only to start
- No macOS runner in initial setup (C++ Xcode license complexity deferred)
- Triggers: push + pull request to main
- Use `r-lib/actions` standard templates (`check-standard` or `check-release` workflow)
- Tests run without network access (`plot_random = c(0,0,0)` baseline) — no `skip_on_ci()` needed for main suite

**Code coverage:**
- Service: codecov.io (free for public repos)
- Coverage computed via `covr` package in a separate CI step
- Coverage badge lives in README.md alongside the build status badge
- Target: >= 70% line coverage (from TST-02)

**pkgdown site:**
- Deployment: GitHub Pages from `gh-pages` branch, auto-deployed on every push to main
- Theme: Bootstrap 5, Bootswatch: `materia`
- Home page: README.md rendered as the pkgdown home (single source of truth, no duplicate index.md)
- Content: Function reference + redirect to book documentation; no standalone vignettes in pkgdown (vignettes live in the book)
- `_pkgdown.yml` created with `url`, `template` (bs_version: 5, bootswatch: materia), and `navbar` pointing to book

**README content:**
- Minimal prose intro — brief description of what the package is and what it can do (no code blocks, no API walkthrough)
- Badges: build status (GitHub Actions) + coverage (codecov.io) at the top
- Documentation link: Prominent link to book website `https://willvieira.github.io/book_forest-demography-IPM/` as the primary documentation source
- Install instructions: `devtools::install_github("willvieira/forest-IPM")` (standard, no extras)
- No quick-start code example in README — book vignettes cover this

### Claude's Discretion

- Exact GitHub Actions YAML structure (which `r-lib/actions` steps to chain)
- Whether to use one combined workflow file or separate files for CI check vs pkgdown deploy
- `_pkgdown.yml` navbar sections and reference grouping

### Deferred Ideas (OUT OF SCOPE)

- macOS runner in CI matrix — deferred until Xcode license issue is resolved
- Windows runner — out of scope per PROJECT.md
- Vignettes in pkgdown — book is the documentation; no duplication needed
</user_constraints>

---

## Summary

Phase 7 is a pure infrastructure phase: wire up three orthogonal pieces of CI/CD automation that have no dependency on each other at runtime, only on the package being installable. The r-lib/actions ecosystem provides canonical, battle-tested templates for all three concerns (R CMD check, code coverage via covr + codecov, and pkgdown). The correct approach is to copy those templates verbatim and then narrow them to the locked decisions (Linux only, R release only, no macOS).

Two workflow files are the right split: `R-CMD-check.yaml` (R CMD check + coverage in two separate jobs) and `pkgdown.yaml` (build and deploy pkgdown site). Keeping them separate means pkgdown can deploy even when R CMD check passes but coverage fails, and the two jobs can run concurrently. The pkgdown workflow deploys only on push to main; coverage uploads on both push and PR.

The `_pkgdown.yml` configuration is minimal: url, Bootstrap 5 + Bootswatch materia, a single navbar component pointing to the book, and an automatically-generated reference section. The README is a complete rewrite of the current 2-line placeholder to carry badges, a brief package description, install instructions, and a book link.

**Primary recommendation:** Use two workflow files. Copy r-lib/actions examples verbatim, then delete the macOS and Windows matrix rows. Add `extra-packages: covr, xml2` to the coverage job's setup-r-dependencies step. Add `covr` and `pkgdown` to DESCRIPTION `Suggests`. Enable GitHub Pages on the `gh-pages` branch in the repository settings before the first pkgdown deploy runs.

---

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| r-lib/actions/setup-r | v2 | Installs R on runner | Official r-lib action, handles RSPM mirror |
| r-lib/actions/setup-r-dependencies | v2 | Installs DESCRIPTION deps + extras | Handles pak, system deps on Ubuntu, caching |
| r-lib/actions/check-r-package | v2 | Runs R CMD check | Official action, uploads snapshots on failure |
| covr | CRAN current | Measures R test coverage | De facto standard for R package coverage |
| codecov/codecov-action | v5 | Uploads coverage to codecov.io | Official Codecov GitHub Action |
| pkgdown | CRAN current (2.x) | Generates static documentation site | Official r-lib documentation tool |
| JamesIves/github-pages-deploy-action | v4.5.0 | Deploys docs/ to gh-pages branch | Used by r-lib/actions pkgdown template |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| xml2 | CRAN current | Required by covr for Cobertura XML output | Always paired with covr in coverage workflow |
| actions/checkout | v4 | Checks out repo code | First step in every workflow job |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| JamesIves/github-pages-deploy-action | actions/deploy-pages (GitHub native) | r-lib templates still use JamesIves; switching requires `pages` permission block and `actions/upload-pages-artifact` step — extra complexity with no benefit |
| codecov.io | coveralls.io | Both free for public repos; codecov is more commonly used in R ecosystem and has better badge support |

**Installation (add to DESCRIPTION Suggests):**
```r
# Add to DESCRIPTION Suggests field:
# covr,
# pkgdown
```

---

## Architecture Patterns

### Recommended File Structure

```
.github/
└── workflows/
    ├── R-CMD-check.yaml      # R CMD check + code coverage (two jobs)
    └── pkgdown.yaml          # pkgdown build + gh-pages deploy
_pkgdown.yml                  # pkgdown site configuration (repo root)
README.md                     # Full rewrite: badges + intro + book link
```

### Pattern 1: Separate Workflow Files

**What:** Two `.yaml` files — one for code quality (check + coverage), one for documentation deployment.

**When to use:** Always for this phase. Decoupling means a coverage upload failure does not block the pkgdown deploy, and the two workflows run concurrently on each push.

**R-CMD-check.yaml:**
```yaml
# Source: https://github.com/r-lib/actions/blob/v2/examples/check-release.yaml
# Modified: removed macOS and Windows matrix rows per locked decisions

on:
  push:
    branches: [main, master]
  pull_request:

name: R-CMD-check.yaml

permissions: read-all

jobs:
  R-CMD-check:
    runs-on: ubuntu-latest
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
      R_KEEP_PKG_SOURCE: yes
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::rcmdcheck
          needs: check
      - uses: r-lib/actions/check-r-package@v2
        with:
          upload-snapshots: true
          build_args: 'c("--no-manual", "--compact-vignettes=gs+qpdf")'
```

**test-coverage.yaml (separate job or second workflow):**
```yaml
# Source: https://github.com/r-lib/actions/blob/v2/examples/test-coverage.yaml

on:
  push:
    branches: [main, master]
  pull_request:

name: test-coverage.yaml

permissions: read-all

jobs:
  test-coverage:
    runs-on: ubuntu-latest
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::covr, any::xml2
          needs: coverage
      - name: Test coverage
        run: |
          cov <- covr::package_coverage(
            quiet = FALSE,
            clean = FALSE,
            install_path = file.path(normalizePath(Sys.getenv("RUNNER_TEMP"), winslash = "/"), "package")
          )
          covr::to_cobertura(cov)
        shell: Rscript {0}
      - uses: codecov/codecov-action@v5
        with:
          fail_ci_if_error: ${{ github.event_name != 'pull_request' && true || false }}
          file: ./cobertura.xml
          plugin: noop
          disable_search: true
          token: ${{ secrets.CODECOV_TOKEN }}
```

### Pattern 2: pkgdown Workflow

**What:** Builds pkgdown site and deploys to `gh-pages` branch. Build runs on all triggers; deploy only on push to main.

```yaml
# Source: https://github.com/r-lib/actions/blob/v2/examples/pkgdown.yaml

on:
  push:
    branches: [main, master]
  pull_request:
  release:
    types: [published]
  workflow_dispatch:

name: pkgdown.yaml

permissions: read-all

concurrency:
  group: pkgdown-${{ github.event_name != 'pull_request' || github.run_id }}
  cancel-in-progress: true

jobs:
  pkgdown:
    runs-on: ubuntu-latest
    concurrency:
      group: pkgdown-${{ github.ref }}
      cancel-in-progress: true
    env:
      GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
    permissions:
      contents: write
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-pandoc@v2
      - uses: r-lib/actions/setup-r@v2
        with:
          use-public-rspm: true
      - uses: r-lib/actions/setup-r-dependencies@v2
        with:
          extra-packages: any::pkgdown, local::.
          needs: website
      - name: Build site
        run: pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
        shell: Rscript {0}
      - name: Deploy to GitHub pages
        if: github.event_name != 'pull_request'
        uses: JamesIves/github-pages-deploy-action@v4.5.0
        with:
          clean: false
          branch: gh-pages
          folder: docs
```

### Pattern 3: _pkgdown.yml

**What:** Minimal configuration — url, Bootstrap 5 materia theme, navbar with Book link, no vignettes section.

```yaml
# Source: https://pkgdown.r-lib.org/articles/customise.html

url: https://willvieira.github.io/forest-IPM/

template:
  bootstrap: 5
  bootswatch: materia

navbar:
  structure:
    left:  [intro, reference]
    right: [book, github]
  components:
    book:
      text: "Book"
      href: https://willvieira.github.io/book_forest-demography-IPM/

reference:
  - title: "IPM Constructors"
    desc: "Functions to build IPM components"
    contents:
      - starts_with("stand")
      - starts_with("species_model")
      - starts_with("env_condition")
      - starts_with("competition")
      - starts_with("parameters")
  - title: "IPM Engines"
    desc: "Functions to run IPM calculations"
    contents:
      - starts_with("ipm_lambda")
      - starts_with("ipm_projection")
  - title: "Utilities"
    desc: "Helper functions"
    contents:
      - supported_species
```

### Pattern 4: README Badges

**What:** Two badges at top of README — one for CI status, one for coverage.

```markdown
<!-- Source: GitHub native badge + codecov.io docs -->
[![R-CMD-check](https://github.com/willvieira/forest-IPM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/willvieira/forest-IPM/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/willvieira/forest-IPM/graph/badge.svg)](https://app.codecov.io/github/willvieira/forest-IPM)
```

### Anti-Patterns to Avoid

- **Single workflow file combining check + coverage + pkgdown:** Coverage failures block documentation deploys; concurrently running jobs is cleaner.
- **index.md alongside README.md:** pkgdown uses README.md as home when no index.md exists. Adding index.md creates duplication and breaks the "single source of truth" constraint.
- **Hardcoding CODECOV_TOKEN as required:** For public repos, `secrets.CODECOV_TOKEN` can be set or left as empty string — the action handles it. Marking it required blocks CI for repos where the secret hasn't been set yet.
- **`--no-build-vignettes` in R CMD check args:** The package currently has no vignettes; omitting this flag means check correctly reports "no vignettes found" rather than silently skipping.

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| System dep detection on Ubuntu | apt-get install build-essential manually | r-lib/actions/setup-r-dependencies | pak auto-detects system deps for Rcpp/RcppEigen on Ubuntu; manual apt steps are fragile |
| Coverage measurement | Custom test runner with line counting | covr::package_coverage() | Handles S3, Rcpp, and source-level tracking correctly |
| XML coverage report for codecov | Manual XML generation | covr::to_cobertura() | Codecov requires Cobertura format; to_cobertura() produces the exact schema expected |
| pkgdown home page duplication | Separate index.md with same content | Use README.md directly | pkgdown natively reads README.md as home page when no index.md is present |
| Badge URLs | Custom badge service | GitHub native SVG + codecov.io native | GitHub's own `badge.svg` reflects actual workflow state; no third-party lag |

**Key insight:** All three CI concerns (check, coverage, docs) have official r-lib templates that handle edge cases (Rcpp compilation, snapshot uploads, gh-pages branch creation). Copy templates and narrow to locked constraints rather than building from scratch.

---

## Common Pitfalls

### Pitfall 1: GitHub Pages Not Enabled Before First Deploy

**What goes wrong:** The pkgdown workflow runs successfully, `JamesIves/github-pages-deploy-action` pushes to `gh-pages` branch, but the site returns 404. The branch exists but GitHub Pages is not configured to serve from it.

**Why it happens:** GitHub Pages must be manually enabled in repository Settings > Pages > Source = "Deploy from a branch" > Branch = `gh-pages` / `/(root)`. The workflow cannot do this automatically.

**How to avoid:** Enable GitHub Pages in repository settings before or immediately after the first workflow run. The `usethis::use_pkgdown_github_pages()` helper does this automatically if run interactively, but since the planner will be writing files directly, manual enablement is required.

**Warning signs:** Workflow shows green but `https://willvieira.github.io/forest-IPM/` returns 404 or "There isn't a GitHub Pages site here."

### Pitfall 2: DESCRIPTION Missing `Suggests: covr, pkgdown`

**What goes wrong:** `setup-r-dependencies` with `extra-packages: any::covr, any::xml2` installs these packages outside of DESCRIPTION, but `R CMD check --as-cran` will warn about packages used in CI but not declared.

**Why it happens:** DESCRIPTION is the package's dependency manifest. CI-only packages belong in `Suggests`, not `Imports`.

**How to avoid:** Add to DESCRIPTION:
```
Suggests:
    covr,
    pkgdown
```

### Pitfall 3: `_pkgdown.yml` Missing `url` Field

**What goes wrong:** pkgdown builds the site but all internal links use relative paths; external `pkgdown::in_pkgdown()` checks fail; auto-linking between packages doesn't work.

**Why it happens:** Without `url:`, pkgdown cannot generate absolute URLs for cross-references or the canonical `<link>` tag.

**How to avoid:** Always include `url: https://willvieira.github.io/forest-IPM/` as the first key in `_pkgdown.yml`.

### Pitfall 4: Rcpp/RcppEigen Compilation on Ubuntu

**What goes wrong:** R CMD check fails with "compilation of C++ code failed" because `g++` or `libstdc++-dev` is not installed.

**Why it happens:** Ubuntu runners on GitHub Actions have `build-essential` available but pak's automatic system dependency detection may not catch all necessary headers for RcppEigen.

**How to avoid:** `r-lib/actions/setup-r-dependencies` with `use-public-rspm: true` generally handles this for CRAN packages. If compilation still fails, add an explicit step:
```yaml
- name: Install system dependencies
  run: sudo apt-get install -y build-essential
```
This is the fallback; the standard template usually handles it without manual intervention.

**Warning signs:** CI error containing `g++: not found` or `fatal error: Eigen/Core: No such file or directory`.

### Pitfall 5: Coverage Job Measuring < 70% Due to Rcpp Code

**What goes wrong:** `covr::package_coverage()` only measures R code coverage, not C++ coverage. The `src/eigen.cpp` eigenvalue solver is not counted.

**Why it happens:** covr instruments R source code; C++ coverage requires a separate toolchain (gcov/lcov) that is not part of the standard setup.

**How to avoid:** This is expected behavior. The 70% target (TST-02) is measured against R code only. Document this in comments near the coverage job definition. Do not attempt C++ coverage — it adds significant complexity for marginal gain on a single-function C++ file.

### Pitfall 6: pkgdown Navbar `github` Component Missing

**What goes wrong:** The `github` icon in the navbar does not appear because pkgdown cannot find the GitHub repo URL.

**Why it happens:** pkgdown reads the GitHub URL from DESCRIPTION's `URL` field or the `url` field in `_pkgdown.yml`. If neither points to the repo, the github component is silently dropped.

**How to avoid:** Add the GitHub repo URL to DESCRIPTION:
```
URL: https://github.com/willvieira/forest-IPM
```
Or keep only `[book, github]` in the right navbar and let pkgdown infer from the `url` field — but the DESCRIPTION `URL` field is the most reliable source.

---

## Code Examples

### Verified: Codecov badge URL pattern
```markdown
<!-- Source: https://docs.codecov.com/docs/status-badges -->
[![codecov](https://codecov.io/github/willvieira/forest-IPM/graph/badge.svg)](https://app.codecov.io/github/willvieira/forest-IPM)
```

### Verified: GitHub Actions native badge URL pattern
```markdown
<!-- Source: GitHub documentation — badge is at /actions/workflows/{filename}/badge.svg -->
[![R-CMD-check](https://github.com/willvieira/forest-IPM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/willvieira/forest-IPM/actions/workflows/R-CMD-check.yaml)
```

### Verified: pkgdown Bootstrap 5 + Bootswatch config
```yaml
# Source: https://pkgdown.r-lib.org/articles/customise.html
template:
  bootstrap: 5
  bootswatch: materia
```

### Verified: External navbar link in pkgdown
```yaml
# Source: https://pkgdown.r-lib.org/articles/customise.html
navbar:
  structure:
    right: [book, github]
  components:
    book:
      text: "Book"
      href: https://willvieira.github.io/book_forest-demography-IPM/
```

### Verified: covr + cobertura XML pipeline
```r
# Source: https://github.com/r-lib/actions/blob/v2/examples/test-coverage.yaml
cov <- covr::package_coverage(
  quiet = FALSE,
  clean = FALSE,
  install_path = file.path(normalizePath(Sys.getenv("RUNNER_TEMP"), winslash = "/"), "package")
)
covr::to_cobertura(cov)
```

### Verified: pkgdown build for GitHub Pages
```r
# Source: https://pkgdown.r-lib.org/reference/build_site_github_pages.html
pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Travis CI for R packages | GitHub Actions (r-lib/actions) | 2020-2021 | Travis CI deprecated free tier; r-lib provides drop-in templates |
| pkgdown deploy via `pkgdown::deploy_site_github()` | `JamesIves/github-pages-deploy-action@v4.5.0` | ~2022 | More reliable, works with modern GitHub Pages |
| covr + coveralls.io | covr + codecov.io via Cobertura XML | 2022-2023 | Codecov action v4+ uses Cobertura; codecov more popular in R ecosystem |
| Bootstrap 4 pkgdown themes | Bootstrap 5 themes | pkgdown 2.0 (2022) | BS5 is current; BS4 still works but new setups should use BS5 |
| `check-standard` workflow (matrix) | `check-release` workflow (single runner) | r-lib/actions v2 | check-release = Linux only, R release — exactly what is locked |

**Deprecated/outdated:**
- `r-lib/actions@v1`: superseded by v2; use `@v2` for all steps
- `covr::codecov()`: replaced by `covr::to_cobertura()` + `codecov/codecov-action`; the old direct upload function is fragile
- `pkgdown::deploy_site_github()`: still works but the GitHub Actions approach is preferred; requires SSH key setup which is more complex

---

## Open Questions

1. **CODECOV_TOKEN secret setup**
   - What we know: Public repos don't need a token for badge display. The `codecov-action@v5` accepts `secrets.CODECOV_TOKEN` but won't fail if the secret is absent for public repos.
   - What's unclear: Whether codecov.io requires the repository to be manually "activated" in the codecov.io dashboard before the first upload, or if the first upload auto-activates.
   - Recommendation: The planner should note that the `secrets.CODECOV_TOKEN` secret must be added to GitHub repository settings (Settings > Secrets > Actions) after creating a codecov.io account and activating the repo. This is a one-time manual step.

2. **`_pkgdown.yml` reference section grouping**
   - What we know: The package exports ~10-15 functions based on NAMESPACE; exact function list needs to be read from NAMESPACE to write the reference grouping.
   - What's unclear: Whether all exported functions are documented well enough for pkgdown to render without errors.
   - Recommendation: Read `/Users/wvieira/GitHub/forest-IPM/NAMESPACE` during planning to enumerate all exported functions for the reference section grouping in `_pkgdown.yml`.

3. **`URL` field in DESCRIPTION**
   - What we know: DESCRIPTION currently has no `URL` field. The `github` navbar component and pkgdown cross-referencing both benefit from it.
   - What's unclear: Whether adding `URL` to DESCRIPTION needs to go through `usethis::use_github_links()` or can be added manually.
   - Recommendation: Add manually: `URL: https://github.com/willvieira/forest-IPM, https://willvieira.github.io/forest-IPM/` — both the GitHub repo and the pkgdown site URL, comma-separated.

---

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework | testthat (version from DESCRIPTION Suggests) |
| Config file | `tests/testthat.R` |
| Quick run command | `Rscript -e "devtools::test()"` |
| Full suite command | `Rscript -e "devtools::test()"` |

### Phase Requirements → Test Map

This phase is pure CI/CD infrastructure — no R source code changes. The deliverables are YAML files, a YAML configuration file, and a Markdown file. Standard R unit tests do not apply to these artifacts. Validation is via CI execution.

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| TST-01 | GitHub Actions CI runs on push | smoke (CI execution) | Push to main branch, observe workflow green | N/A — CI execution |
| TST-02 | Coverage >= 70% line coverage | measured by covr | `Rscript -e "covr::package_coverage()"` | ✅ existing tests |
| DOC-04 | pkgdown site renders without errors | smoke (local build) | `Rscript -e "pkgdown::build_site()"` | ❌ Wave 0: `_pkgdown.yml` |

### Sampling Rate

- **Per task commit:** `Rscript -e "devtools::test()"` (verify existing tests still pass after file additions)
- **Per wave merge:** `Rscript -e "pkgdown::build_site()"` + `Rscript -e "covr::package_coverage()"`
- **Phase gate:** Full suite green + pkgdown builds locally before `/gsd:verify-work`

### Wave 0 Gaps

- [ ] `_pkgdown.yml` — required before pkgdown smoke test can run; created in Wave 0 of this phase
- [ ] `.github/workflows/R-CMD-check.yaml` — CI execution requires this file
- [ ] `.github/workflows/pkgdown.yaml` — pkgdown CI deploy requires this file

---

## Sources

### Primary (HIGH confidence)
- `https://github.com/r-lib/actions/blob/v2/examples/check-release.yaml` — R CMD check workflow template; steps and triggers verified
- `https://github.com/r-lib/actions/blob/v2/examples/test-coverage.yaml` — covr + codecov workflow template; steps and env vars verified
- `https://github.com/r-lib/actions/blob/v2/examples/pkgdown.yaml` — pkgdown workflow template; deploy action version verified
- `https://pkgdown.r-lib.org/articles/customise.html` — Bootstrap 5 config, navbar external links, bootswatch keys verified
- `https://pkgdown.r-lib.org/reference/build_site_github_pages.html` — build_site_github_pages() parameters verified
- `https://docs.codecov.com/docs/status-badges` — codecov badge URL pattern verified; public repo token requirement confirmed

### Secondary (MEDIUM confidence)
- `https://github.com/r-lib/actions` (README) — action names and v2 tag confirmed; specific YAML details from example files
- GitHub Actions documentation — native badge URL pattern `/{owner}/{repo}/actions/workflows/{file}/badge.svg` confirmed via multiple sources

### Tertiary (LOW confidence)
- WebSearch results on RcppEigen system deps on Ubuntu — `setup-r-dependencies` auto-handles this for CRAN packages; manual `apt-get` fallback noted but not tested against this specific package

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — r-lib/actions v2 templates verified against official GitHub source; codecov action v5 confirmed current
- Architecture: HIGH — workflow YAML patterns verified from official r-lib/actions examples; pkgdown config verified from official docs
- Pitfalls: MEDIUM — GitHub Pages enablement and Rcpp compilation pitfalls from official docs and known ecosystem behavior; Rcpp-specific Ubuntu interaction not directly tested

**Research date:** 2026-03-13
**Valid until:** 2026-04-13 (r-lib/actions v2 is stable; pkgdown 2.x is stable; codecov API changes infrequently)
