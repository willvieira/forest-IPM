---
plan: 04-03
title: Rewrite Book Guide to forestIPM API + New Deep Dive Chapter
status: complete
completed: 2026-03-13
---

## What Was Built

Rewrote `guide_IPM.qmd` to use the new `forestIPM` package API, created a new `deep_dive_IPM.qmd` chapter, and restructured `_quarto.yml` with a three-section book layout. All changes committed to the book repo (`book_forest-demography-IPM`).

## Key Files

### Modified (book repo)
- `guide_IPM.qmd` — rewritten to use `library(forestIPM)` with the five-constructor chain (`stand`, `species_model`, `parameters`, `project`) and `plot(proj, type = ...)` for visualizations. Old `source()`-from-GitHub and `getPars_sp()`/`mkKernel()`/`eigen()` patterns removed.
- `_quarto.yml` — restructured into three `part:` sections: "Demographic models evaluation", "From demographic rates to population dynamics", "Using the IPM". `deep_dive_IPM.qmd` added to section 2.

### Created (book repo)
- `deep_dive_IPM.qmd` — new chapter with multi-species IPM example, time-varying climate pattern, and roadmap section (ML decoder, coexistence IPM, CRAN submission).

## Verification

- `library(forestIPM)` present in `guide_IPM.qmd`
- `deep_dive_IPM.qmd` exists with multi-species example and roadmap
- `_quarto.yml` has all three `part:` sections with `deep_dive_IPM.qmd` in section 2
- Human checkpoint approved; changes committed to book repo

## Deviations

None — plan executed as specified.
