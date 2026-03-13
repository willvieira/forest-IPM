# Phase 7: Add GitHub Actions CI with code coverage and pkgdown site with README pointing to book vignettes - Context

**Gathered:** 2026-03-13
**Status:** Ready for planning

<domain>
## Phase Boundary

Set up CI/CD infrastructure and package documentation site: (1) GitHub Actions CI on Linux with R release + code coverage via codecov.io, (2) pkgdown site deployed to GitHub Pages on every push to main, (3) update README.md with badges (build status, coverage) and a brief package intro linking to the book documentation site.

</domain>

<decisions>
## Implementation Decisions

### CI matrix

- **Linux only** — Ubuntu runner, R release only to start
- No macOS runner in initial setup (C++ Xcode license complexity deferred)
- Triggers: push + pull request to main
- Use `r-lib/actions` standard templates (`check-standard` or `check-release` workflow)
- Tests run without network access (`plot_random = c(0,0,0)` baseline) — no `skip_on_ci()` needed for main suite

### Code coverage

- **Service:** codecov.io (free for public repos)
- Coverage computed via `covr` package in a separate CI step
- Coverage badge lives in README.md alongside the build status badge
- Target: ≥ 70% line coverage (from TST-02)

### pkgdown site

- **Deployment:** GitHub Pages from `gh-pages` branch, auto-deployed on every push to main
- **Theme:** Bootstrap 5, Bootswatch: `materia`
- **Home page:** README.md rendered as the pkgdown home (single source of truth, no duplicate index.md)
- **Content:** Function reference + redirect to book documentation; no standalone vignettes in pkgdown (vignettes live in the book)
- `_pkgdown.yml` created with `url`, `template` (bs_version: 5, bootswatch: materia), and `navbar` pointing to book

### README content

- **Minimal prose intro** — brief description of what the package is and what it can do (no code blocks, no API walkthrough)
- **Badges:** build status (GitHub Actions) + coverage (codecov.io) at the top
- **Documentation link:** Prominent link to book website `https://willvieira.github.io/book_forest-demography-IPM/` as the primary documentation source
- **Install instructions:** `devtools::install_github("willvieira/forest-IPM")` (standard, no extras)
- No quick-start code example in README — book vignettes cover this

### Claude's Discretion

- Exact GitHub Actions YAML structure (which `r-lib/actions` steps to chain)
- Whether to use one combined workflow file or separate files for CI check vs pkgdown deploy
- `_pkgdown.yml` navbar sections and reference grouping

</decisions>

<specifics>
## Specific Ideas

- Book URL: `https://willvieira.github.io/book_forest-demography-IPM/`
- Package installed via: `devtools::install_github("willvieira/forest-IPM")`
- pkgdown theme: `bootstrap: 5`, `bootswatch: materia`
- README should be intentionally minimal — the book is the documentation; README just points there

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets

- `tests/testthat/` — existing test suite (test-constructors.R, test-regression-baselines.R, test-workflow.R) runs without network; safe for CI
- `DESCRIPTION` — all Imports already declared; `Suggests` may need `covr`, `pkgdown` added for CI tooling

### Established Patterns

- Tests use `plot_random = c(0, 0, 0)` — no cloud calls, no `skip_on_ci()` needed
- Rcpp C++ compilation happens at `devtools::load_all()` — Linux runner will need `build-essential` (handled by `r-lib/actions/setup-r-dependencies`)
- `roxygen2::roxygenise(load_code='source')` avoids Xcode requirement for docs — not relevant here since CI uses Linux

### Integration Points

- `.github/workflows/` — new directory; two workflow files likely: one for R CMD check + coverage, one for pkgdown deploy
- `_pkgdown.yml` — new file at repo root
- `README.md` — currently a 2-line placeholder; full rewrite needed

</code_context>

<deferred>
## Deferred Ideas

- macOS runner in CI matrix — deferred until Xcode license issue is resolved
- Windows runner — out of scope per PROJECT.md
- Vignettes in pkgdown — book is the documentation; no duplication needed

</deferred>

---

*Phase: 07-add-github-actions-ci-with-code-coverage-and-pkgdown-site-with-readme-pointing-to-book-vignettes*
*Context gathered: 2026-03-13*
