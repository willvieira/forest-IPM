---
phase: 01-api-design
verified: 2026-03-03T00:00:00Z
status: passed
score: 7/7 must-haves verified
re_verification: false
gaps: []
human_verification:
  - test: "Read API_DESIGN.md Section 3 end-to-end as a Phase 2 implementor"
    expected: "All eight functions feel fully specified — no gaps requiring opening CONTEXT.md"
    why_human: "Self-containment is a subjective judgment of completeness; programmatic checks cover structure but not readability"
---

# Phase 1: API Design Verification Report

**Phase Goal:** Produce a complete, locked API design document that Phase 2 implementors can follow without consulting any other file
**Verified:** 2026-03-03
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|---------|
| 1 | Phase 2 implementors can work from API_DESIGN.md alone — no other file needed | VERIFIED | Audience line on line 5: "read this document alone; no need to open CONTEXT.md, RESEARCH.md, or the architecture doc." All 10 sections are self-contained. Two Section 9 mentions of CONTEXT.md are historical rationale only — they explain why decisions were made, not instructions to read that file. |
| 2 | Every exported function has a finalized signature, argument names, and return type | VERIFIED | All 8 functions have `**Signature:**` lines, argument tables, and return value sections confirmed by grep. Signatures: `stand(data)`, `species_model(x, on_missing = "error")`, `parameters(mod, draw = "random", seed = NULL)`, `env_condition(MAT, MAP)`, `control(years = 100, delta_time = 1, store_every = 1, bin_width = 1)`, `lambda(mod, pars, stand, env, ctrl)`, `project(mod, pars, stand, env, ctrl)`, `supported_species()`. |
| 3 | The ipm_projection object structure is fully specified (field names, types, accessor patterns) | VERIFIED | Section 4.6 provides a complete field table: `$species` (character vector), `$years` (numeric vector, length = ceiling(ctrl$years / ctrl$store_every)), `$lambda` (named list of numeric vectors), `$stand_series` (list of ipm_stand objects), `$summary` (tibble). Access patterns shown with code block. Section header explicitly states "all fields required — no TBD". |
| 4 | Naming conventions are codified — snake_case functions, ipm_ prefix for S3 classes, dot prefix for internal helpers | VERIFIED | Section 5 (line 567) provides a complete table covering: Exported functions (snake_case), S3 class names ("ipm_" prefix + noun), S3 methods (<generic>.<class>), Low-level constructors (new_<class>), Validators (validate_<class>), Internal helpers (. prefix), Type checks (inherits(x, "ipm_stand")). Anti-pattern section also listed. |
| 5 | All validation rules enumerate the exact cli_abort() message template per constructor | VERIFIED | Section 7 (line 629) contains complete validation table for stand(), species_model(), parameters(), env_condition(), control(), and cross-constructor engine rules. All use named vector cli_abort() format with {.arg}, {.field}, {.val}, {.run}, {.cls}, {.code} tokens. 42 token occurrences confirmed by grep. |
| 6 | The canonical workflow supersedes the original run_ipm() concept — API-03 is explicitly addressed | VERIFIED | API-03 note appears at line 84 (Section 2) and line 731 (Section 9, Decision 1). Section 2 note states: "This canonical workflow IS the fulfillment of API-03." Section 9 states: "Phase 2 implementors MUST NOT implement a run_ipm() function." REQUIREMENTS.md API-03 entry is marked [x] with the same language. |
| 7 | env_condition() vs env_conditions() naming conflict resolved to a single canonical form | VERIFIED | Line 280: "**Canonical name:** env_condition() — singular. **NEVER use env_conditions() (plural).**" CONTEXT.md had used the plural form (line 73); API_DESIGN.md explicitly resolves this. Line 753 in Section 9 provides the rationale. |

**Score:** 7/7 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `.planning/phases/01-api-design/API_DESIGN.md` | Complete API interface contract for Phase 2 | VERIFIED | File exists at 787 lines (confirmed by wc -l). Substantially beyond the 150-line threshold in the PLAN's automated verify step. Contains all 10 sections. No TBD placeholders found (the only "TBD" occurrence is in the phrase "no TBD" on line 529). |

**Level 1 (Exists):** The file is present at the declared path.

**Level 2 (Substantive):** 787 lines covering 10 complete sections: Overview and Design Principles, Canonical Workflow, Exported Functions Reference (all 8), S3 Classes (all 6), Naming Conventions, Constructor/Validator/Helper Pattern, Validation Rules (complete table), Internal Functions, Design Decisions and Rationale, Versioning and Extension Points.

**Level 3 (Wired):** This is a documentation artifact. "Wired" means it is self-contained and forward-linked. Confirmed: REQUIREMENTS.md traceability table references it, ROADMAP.md Phase 1 plan entry is marked [x], and SUMMARY.md explicitly records the document path and line count.

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| CONTEXT.md locked decisions | API_DESIGN.md function tables | Transcription of all 8 exported functions with finalized signatures | VERIFIED | All 8 functions from CONTEXT.md "Public API surface" table are present in API_DESIGN.md Section 3 with expanded detail. Two discrepancies from CONTEXT.md were correctly resolved: (1) CONTEXT.md used class name "ipm_model"; API_DESIGN.md uses "ipm_spModel" throughout as explicitly required by the PLAN. (2) CONTEXT.md section heading used "env_conditions()" plural; API_DESIGN.md resolves to "env_condition()" singular with an explicit prohibition note. |
| API_DESIGN.md | Phase 2 implementation tasks | Self-contained reference — no back-references to CONTEXT.md or RESEARCH.md needed | VERIFIED | Audience declaration on line 5 states "read this document alone." The two mentions of CONTEXT.md in Section 9 are retrospective rationale explaining why decisions were made (e.g., "CONTEXT.md locked this"). They do not direct readers to open other files. REQUIREMENTS.md and ROADMAP.md both mark Phase 1 as Complete and reference API_DESIGN.md as the Phase 2 input. |

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|-------------|-------------|--------|---------|
| API-01 | 01-01-PLAN.md | All exported functions listed with finalized signatures, argument names, and return types | SATISFIED | Section 3 contains all 8 functions each with: `**Signature:**` line, arguments table (name, type, default, description), return value table. Confirmed by grep showing one Signature occurrence per function. |
| API-02 | 01-01-PLAN.md | supported_species() — returns data frame of available species IDs, names, and model variants | SATISFIED | Section 3.8 specifies a tibble with exactly 6 columns: species_id, common_name, nom_commun, growth_model, surv_model, recruit_model. Note: PLAN interfaces block described 3 columns (species_id, common_name, model_variant); the executor expanded "model_variant" into 3 separate model-type columns (growth, survival, recruitment). This is a refinement that provides more specificity for Phase 2 — not a contradiction. REQUIREMENTS.md marks API-02 [x]. |
| API-03 | 01-01-PLAN.md | run_ipm() superseded by five-constructor + two-engine pattern; canonical workflow documented | SATISFIED | API-03 supersession note in Section 2 (line 84) and Section 9, Decision 1 (line 721). REQUIREMENTS.md API-03 entry updated to reflect supersession. |
| API-04 | 01-01-PLAN.md | Naming conventions documented and applied consistently | SATISFIED | Section 5 (line 567) with complete table. Conventions applied consistently across all sections — functions use snake_case, classes use ipm_ prefix, internals use dot prefix throughout the document. |
| API-05 | 01-01-PLAN.md | API design document produced that guides Phase 2 | SATISFIED | The artifact itself is the requirement fulfillment. The document is marked "Final — Phase 1 deliverable" and explicitly states the Phase 2 audience on line 5. |

**Orphaned requirements:** None. Only API-01 through API-05 are mapped to Phase 1 in REQUIREMENTS.md. No additional Phase 1 assignments were found.

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| API_DESIGN.md | 755 | "The singular form appears in the public API table in CONTEXT.md and is treated as locked." | Info | Rationale text references another file for historical context only; does not require the reader to open it. No impact on self-containment. |
| API_DESIGN.md | 761 | "CONTEXT.md (the more recent locked-decision document) specifies plot_size as a column..." | Info | Same pattern — Section 9 rationale referencing source documents historically. No action needed. |

No blockers. No TBD placeholders. No stub implementations. No empty validation tables.

---

### Commit Verification

Both commits declared in SUMMARY.md are confirmed real:

| Commit | Hash | Date | Content |
|--------|------|------|---------|
| feat(01-01): write API_DESIGN.md | 7e294f8 | 2026-03-02 | 787-line insertion creating API_DESIGN.md |
| docs(01-01): fix stand() size constraint | 17e1981 | 2026-03-02 | 8-line diff correcting >= 127 mm threshold in 4 locations |

---

### Human Verification Required

#### 1. Self-containment Readability Check

**Test:** Read API_DESIGN.md from Section 1 through Section 10 as a Phase 2 implementor who has never seen CONTEXT.md or RESEARCH.md.
**Expected:** All questions from the PLAN's `<success_criteria>` block are answerable without opening any other file. The document should feel complete enough to begin writing R code.
**Why human:** "Completeness for implementation" is a subjective judgment. Programmatic checks confirm structure and content presence, but cannot assess whether the document is clear enough for a developer who has no prior context.

#### 2. supported_species() Column Count Alignment

**Test:** Confirm with the package author whether `supported_species()` should return 3 columns (species_id, common_name, model_variant as a single column) or 6 columns (species_id, common_name, nom_commun, growth_model, surv_model, recruit_model as separate columns).
**Expected:** Either the 6-column design is confirmed as correct, or the document is corrected to 3 columns.
**Why human:** The PLAN interfaces block specified 3 columns; API_DESIGN.md documents 6. This is likely a deliberate and correct expansion (the actual data has separate model strings per vital rate), but the discrepancy between the PLAN specification and the document should be confirmed by the author before Phase 2 implementors build to it.

---

### Gaps Summary

No gaps found. All 7 observable truths are verified. The single artifact exists, is substantive (787 lines), and is wired into the planning system. Both key links are verified. All 5 requirements (API-01 through API-05) are satisfied with evidence. No orphaned requirements. No blocker anti-patterns.

The only open item is the human verification for `supported_species()` column count (3 vs 6), which is a clarification rather than a gap — the document is internally consistent in documenting 6 columns throughout, and the 6-column design is well-motivated by the actual vital-rate structure of the data.

---

_Verified: 2026-03-03_
_Verifier: Claude (gsd-verifier)_
