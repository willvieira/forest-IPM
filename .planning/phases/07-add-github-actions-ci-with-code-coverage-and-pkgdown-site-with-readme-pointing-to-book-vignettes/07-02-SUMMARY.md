---
phase: 07-add-github-actions-ci-with-code-coverage-and-pkgdown-site-with-readme-pointing-to-book-vignettes
plan: "02"
subsystem: documentation
tags: [pkgdown, github-actions, github-pages, documentation-site]
dependency_graph:
  requires: []
  provides: [_pkgdown.yml, pkgdown-ci-workflow]
  affects: [github-pages, pkgdown-site]
tech_stack:
  added: [pkgdown, JamesIves/github-pages-deploy-action]
  patterns: [bootstrap5-materia, reference-grouping, gh-pages-deploy]
key_files:
  created:
    - _pkgdown.yml
    - .github/workflows/pkgdown.yaml
  modified: []
decisions:
  - "pkgdown site URL is https://willvieira.github.io/forest-IPM/ — must be first key in _pkgdown.yml for canonical link tag"
  - "Bootstrap 5 + Bootswatch materia theme (locked from CONTEXT.md)"
  - "Book navbar link points to https://willvieira.github.io/book_forest-demography-IPM/"
  - "JamesIves/github-pages-deploy-action@v4.5.0 used instead of GitHub native deploy-pages action"
  - "Deploy to gh-pages branch conditional on non-PR events only"
  - "GitHub Pages must be manually enabled in Settings > Pages > Source = gh-pages branch"
metrics:
  duration: 15
  completed_date: "2026-03-14"
  tasks_completed: 2
  files_modified: 2
---

# Phase 7 Plan 02: pkgdown Site Configuration Summary

**One-liner:** pkgdown site configured with Bootstrap 5 materia theme, book navbar link, and all 11 exported functions grouped in reference section; GitHub Actions workflow deploys to gh-pages on push to main.

## Tasks Completed

| Task | Name | Commit | Files Modified |
|------|------|--------|----------------|
| 1 | Create _pkgdown.yml | 67ca965 | _pkgdown.yml |
| 2 | Create pkgdown.yaml GitHub Actions workflow | 0c5f2f5 | .github/workflows/pkgdown.yaml |

## What Was Built

### _pkgdown.yml

Created the pkgdown site configuration at the repository root:

- `url: https://willvieira.github.io/forest-IPM/` as first key (required for absolute URLs and canonical link tag)
- Bootstrap 5 + Bootswatch materia theme
- Navbar with `intro`, `reference` on left and `book`, `github` on right
- Book component: `text: "Book"` pointing to `https://willvieira.github.io/book_forest-demography-IPM/`
- Reference section groups all 11 exported functions:
  - **IPM Constructors**: stand, species_model, env_condition, parameters, control
  - **IPM Engines**: lambda, project
  - **Utilities**: supported_species, scale_env, unscale_env, set_random_effects

### .github/workflows/pkgdown.yaml

Created the GitHub Actions workflow for automated pkgdown builds and deployment:

- Triggers: push to main/master, pull requests, releases, and manual `workflow_dispatch`
- Builds site on all triggers using `pkgdown::build_site_github_pages()`
- Deploys to `gh-pages` branch only on non-PR events (push to main, releases, manual dispatch)
- Uses `JamesIves/github-pages-deploy-action@v4.5.0` with `clean: false`
- Sets `GITHUB_PAT` from `secrets.GITHUB_TOKEN` for private dependency access
- Concurrency groups prevent duplicate deploys

## Post-Deployment Manual Step Required

GitHub Pages must be manually enabled before the site is accessible:
- Go to repository Settings > Pages
- Set Source = "Deploy from a branch"
- Set Branch = `gh-pages` / `/(root)`

The site will return 404 until this is done, even if the workflow runs green.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Installed pkgdown package which was missing**
- **Found during:** Task 1 verification
- **Issue:** `pkgdown` not installed in the local R environment, so `pkgdown::build_site()` verification failed with "there is no package called 'pkgdown'"
- **Fix:** Installed pkgdown and all dependencies (sass, knitr, bslib, downlit, httr2, ragg, rmarkdown, yaml) from CRAN
- **Files modified:** None (system-level R package installation)
- **Commit:** N/A (package installation, not a code change)

**Note on build verification:** The `pkgdown::build_site(preview = FALSE)` verification step was attempted but `Rscript` command execution was blocked by the tool permission system during this session. The `_pkgdown.yml` file matches the exact YAML specification from the plan. The GitHub Actions workflow will perform the authoritative build verification in CI.

## Self-Check

- [x] `_pkgdown.yml` exists at repo root — VERIFIED
- [x] `.github/workflows/pkgdown.yaml` exists — VERIFIED
- [x] `_pkgdown.yml` contains `url: https://willvieira.github.io/forest-IPM/` — VERIFIED
- [x] `_pkgdown.yml` contains `bootswatch: materia` — VERIFIED
- [x] `_pkgdown.yml` contains Book navbar component — VERIFIED
- [x] `_pkgdown.yml` lists all 11 exported functions — VERIFIED
- [x] `pkgdown.yaml` contains `JamesIves/github-pages-deploy-action@v4.5.0` — VERIFIED
- [x] `pkgdown.yaml` deploys conditionally on non-PR events — VERIFIED
- [x] `pkgdown.yaml` targets `gh-pages` branch — VERIFIED
- [x] Commits exist: 67ca965 (_pkgdown.yml), 0c5f2f5 (pkgdown.yaml) — VERIFIED
