---
plan: 04-02
title: Add @examples to All Exported Functions
status: complete
completed: 2026-03-13
---

## What Was Built

Added `@examples` sections to all 11 exported functions in the forestIPM package, enabling `?function_name` to display working examples and `devtools::run_examples()` to pass.

## Key Files

### Modified
- `R/stand.R` — example: construct stand from data.frame
- `R/species_model.R` — example: load ABIBAL parameters
- `R/parameters.R` + `R/set_random_effects` — example: parameters() and RE override
- `R/env_condition.R` — example: env_condition() constructor
- `R/env_scaling.R` — examples: scale_env() / unscale_env() round-trip
- `R/lambda.R` — example: lambda() with stand + env
- `R/project.R` — example: project() with \donttest{} wrapper (slow)
- `R/supported_species.R` — example: supported_species()
- `R/control.R` — example: control() defaults

### Regenerated
- `man/*.Rd` — all 11 man pages now contain `\examples{}` sections

## Verification

- `devtools::run_examples()` passes with 0 errors
- Only pre-existing warnings: roxygen2 version mismatch (7.1.2 vs 7.3.2) and block-quote in stand.R @param

## Deviations

None — plan executed as specified.
