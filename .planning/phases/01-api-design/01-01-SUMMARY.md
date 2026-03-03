---
phase: 01-api-design
plan: "01"
subsystem: api-design
tags: [api, interface-contract, s3-classes, validation, design]

dependency_graph:
  requires: []
  provides:
    - Complete API interface contract for Phase 2 implementation
    - Finalized signatures for all 8 exported functions
    - S3 class specifications for all 6 ipm_* classes
    - Validation error message templates (cli_abort format)
    - Naming conventions (snake_case, ipm_ prefix, dot prefix)
  affects:
    - Phase 2: Package Skeleton and Bug Fixes (all implementation tasks)
    - Phase 3: Remote Data Layer (parameters() fetch interface)
    - Phase 4: Documentation and Tests (roxygen2 docs target this API)

tech_stack:
  added: []
  patterns:
    - Three-layer S3 constructor/validator/helper pattern (Advanced R, Wickham)
    - cli::cli_abort() named-vector error message format with .arg/.field/.val/.run tokens
    - Fail-fast validation at construction time (not deferred to engines)
    - stringdist::amatch() for closest-match species ID suggestions

key_files:
  created:
    - .planning/phases/01-api-design/API_DESIGN.md
  modified: []

decisions:
  - "Five constructors + two engines replaces run_ipm() — fulfills API-03 with composable design"
  - "Single ipm_spModel class for both single and multi-species — avoids API branching"
  - "env_condition() singular is canonical — env_conditions() plural is never used"
  - "plot_size is a column in stand() data argument, not a separate argument"
  - "store_every in control() enables memory management for long (>200 year) runs"
  - "Stochasticity fully resolved in parameters() — lambda() and project() are deterministic"
  - "stand() size constraint is >= 127 mm (minimum DBH threshold for inventory data)"

metrics:
  duration_minutes: 60
  tasks_completed: 2
  tasks_total: 2
  files_created: 1
  files_modified: 1
  completed_date: "2026-03-03"

requirements_addressed:
  - API-01
  - API-02
  - API-03
  - API-04
  - API-05
---

# Phase 1 Plan 01: API Design — Complete Interface Contract Summary

**One-liner:** Self-contained API design document specifying all 8 exported functions, 6 S3 classes, cli_abort() validation templates, and naming conventions for the forestIPM package using the five-constructor + two-engine pattern.

---

## What Was Built

`.planning/phases/01-api-design/API_DESIGN.md` — a 788-line Phase 2 implementation reference covering:

- **8 exported functions** with full signatures, argument tables, return types, and validation rules:
  `stand()`, `species_model()`, `parameters()`, `env_condition()`, `control()`, `lambda()`, `project()`, `supported_species()`
- **6 S3 classes** fully specified: `ipm_stand`, `ipm_spModel`, `ipm_parameters`, `ipm_env`, `ipm_control`, `ipm_projection` — including all field names, types, and print method contracts
- **ipm_projection** complete field specification: `$species`, `$years`, `$lambda`, `$stand_series`, `$summary` (no TBD placeholders)
- **Canonical workflow** replacing `run_ipm()` with the composable five-constructor + two-engine pattern
- **API-03 supersession note** making explicit that `run_ipm()` MUST NOT be implemented
- **Three-layer S3 pattern** (new_/validate_/helper) with ipm_stand code example
- **Validation table** with exact cli_abort() message templates for every constructor
- **Naming conventions**: snake_case functions, ipm_ class prefix, . prefix for internals
- **Design decisions and rationale** for all locked choices
- **Extension points** documented for v2+ (spatial random effects, climate lookup) without implementing them

---

## Tasks Completed

| Task | Name | Commit | Files |
|------|------|--------|-------|
| 1 | Write API_DESIGN.md — function reference, S3 classes, canonical workflow | 7e294f8 | .planning/phases/01-api-design/API_DESIGN.md |
| 2 | Verify API_DESIGN.md is complete and ready for Phase 2 (human-verify checkpoint) | 17e1981 | .planning/phases/01-api-design/API_DESIGN.md |

---

## Deviations from Plan

### Auto-fixed Issues

None — plan executed as written.

### Human-Verify Correction (Task 2)

**Found during:** Task 2 (human verification checkpoint)
**Issue:** The `stand()` size constraint was documented as `> 0 mm` (any positive value) instead of the correct minimum DBH threshold for inventory data.
**Fix:** Updated size constraint to `>= 127 mm` in three locations:
  1. Section 3.1 Required Columns table (constraint column)
  2. Section 3.1 Validation Rules table (condition and error message)
  3. Section 7 Complete Validation Table for `stand()` (condition and error message)
  4. Section 3.1 Usage Example (updated dbh values to be >= 127 mm)
**Files modified:** `.planning/phases/01-api-design/API_DESIGN.md`
**Commit:** 17e1981

---

## Key Decisions Made

1. **Five constructors + two engines replaces `run_ipm()`** — provides composability; researchers can call `lambda()` without full projection; stochasticity is explicit and controllable in `parameters()`. API-03 fulfilled by canonical workflow in Section 2.

2. **Single `ipm_spModel` class** — `species_model()` always returns `"ipm_spModel"` regardless of species count. No `ipm_single_model` or `ipm_community_model` anti-patterns.

3. **`env_condition()` singular is canonical** — locked form; `env_conditions()` plural never used. The function creates one environment condition object.

4. **`plot_size` as a column in data** — CONTEXT.md locked this over the earlier architecture v2 doc that showed it as a separate argument. Column approach keeps plot metadata co-located with tree records.

5. **`store_every` in `control()`** — enables 10x memory reduction for 500+ year runs. Exposed to researchers explicitly.

6. **`stand()` size constraint: >= 127 mm** — corrected from `> 0` during human verification. Reflects the minimum DBH threshold used in forest inventory protocols.

---

## Phase 1 Requirements Coverage

| Requirement | Status | Evidence |
|-------------|--------|---------|
| API-01: All 8 exported functions with finalized signatures + return types | DONE | Section 3 of API_DESIGN.md |
| API-02: supported_species() fully specified (tibble with 6 named columns) | DONE | Section 3.8 of API_DESIGN.md |
| API-03: Canonical workflow; run_ipm() supersession explicit | DONE | Section 2 + Section 9, Decision 1 |
| API-04: Naming conventions section (snake_case, ipm_ prefix, dot prefix) | DONE | Section 5 of API_DESIGN.md |
| API-05: Document IS the Phase 2 deliverable | DONE | API_DESIGN.md self-contained; no other file needed |

---

## Self-Check

**API_DESIGN.md exists:**
- Path: `.planning/phases/01-api-design/API_DESIGN.md`
- Lines: 788

**Commits exist:**
- 7e294f8 — feat(01-01): write API_DESIGN.md
- 17e1981 — docs(01-01): fix stand() size constraint to >= 127 mm

## Self-Check: PASSED
