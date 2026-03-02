# Phase 1: API Design - Research

**Researched:** 2026-03-02
**Domain:** R package API design — S3 class system, function signatures, documentation contracts, error messaging
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Overall API shape

- **Minimal object types**: five core objects + two engine functions; single- and multi-species unified under one class
- Strictly separate: ecological state (`stand`), structural model (`species_model`), parameter realization (`parameters`), climate drivers (`env_condition`), simulation control (`control`)
- No separate `community_params` class — multi-species is just a special case of the same types
- No unnecessary wrapper functions; researchers compose constructors + engine functions directly
- Existing low-level functions (`mkKernel`, `init_pop`, `getPars_sp`, `size_to_BAcomp`, `dbh_to_sizeDist`) made internal

#### Public API surface (final)

**Constructors:**

| Function | Signature | Purpose |
|---|---|---|
| `stand()` | `stand(data)` | Represent forest plot state; `data` is a dataframe with required `size||dbh` in mm + `species ID` + `plot_size` in m2 columns |
| `species_model()` | `species_model(x)` where `x` is a species vector or a `stand` object | Define species-level structure for the IPM (kernel constructors) |
| `parameters()` | `parameters(mod, draw = "random", seed = NULL)` | Resolve a single parameter realization (population and random effects) for all species in `mod` |
| `env_condition()` | `env_condition(MAT, MAP)` | Represent exogenous climate drivers and contains scaling parameters to internally scale natural values to 0-1 scale |
| `control()` | `control(years = 100, delta_time = 1)` | Simulation configuration |

**Core engines:**

| Function | Signature | Purpose |
|---|---|---|
| `lambda()` | `lambda(mod, pars, stand, env)` | Compute dominant eigenvalue per `mod` species |
| `project()` | `project(mod, pars, stand, env, ctrl)` | Project population dynamics over time |

**Utility:**

| Function | Signature | Purpose |
|---|---|---|
| `supported_species()` | `supported_species()` | Returns tibble of available species IDs, names, model variants |

#### `stand()` object

- Signature: `stand(data)` — `data` must contain `size` or `dbh` (numeric > 0), `species` or `species_id` (must match known IDs), and `plot_size` columns
- Internal structure: `stand$trees` (individual-level data), `stand$species` (unique species vector), `stand$plot_size`
- S3 class: `"ipm_stand"`
- `stand` is immutable: `project()` returns new stand states, never modifies in-place

#### `species_model()` object

- Defines structural IPM(s): kernel constructors and parameters (L, Lmax, h bins)
- Can be initialized from an explicit species vector (`species_model(c("Picea_sp", "Abies_sp"))`) or inferred from a stand (`species_model(stand)`)
- Single- and multi-species share the same class (`"ipm_model"`)
- If species without parameters (sp NOT in `supported_species()`) is listed, three options are available: `error` to return an verbose error to the user; `drop` for dropping the species from the list and removing it to the `stand`; `static` to keep the species for background competition estimation but not `lambda` or `project` calculation. Default of `missing` argument to `error`
- Internal: `mod$species` is a character vector

#### `parameters()` object

- Signature: `parameters(mod, draw = "random", seed = NULL)` — takes `mod` as first arg (not species ID)
- `draw`: `"mean"` uses posterior mean; `"random"` samples one draw; integer (1-1000) selects exact posterior draw
- Parameters always stored as multi-species structure internally, even for single species
- Internal structure: `pars$species_params$<species_id>$fixed`, `$random_effects`, `$draw_id`
- Random effects live inside `parameters()` — they represent a model realization, not ecological state
- Future extension: `parameters(mod, draw = 12, random_effect = spatial_re(lat, lon))`

#### `env_conditions()` naming and behavior

- Named `env_condition()` — avoids collision with base R `environment()`
- Accepts scalars or functions of time: `MAT = function(t) 4 + 0.02 * t`
- MAT/MAP scaling to [0,1] applied internally with internal parameters; never exposed to user

#### Engine function signatures

- `lambda(mod, pars, stand, env)` — `mod` is first argument
- `project(mod, pars, stand, env, ctrl)` — `mod` is first argument
- `plot_effects` argument removed from both engines; random effects are resolved inside `parameters()`
- `lambda()` returns a named numeric vector (one λ per species in `mod`)
- `lambda()` does not require parameters for heterospecific species — competition is computed from stand structure only; no dynamic coupling

#### Multi-species projections

- `species_model(c("ABIBAL", "ACERUB"))` + `parameters(mod, ...)` — unified interface
- `project()` with multiple species: fully coupled dynamics (competition updated each timestep)
- `project()` with single species: independent dynamic
- Validation: `pars` must contain parameters for all species in `mod`; species in `stand` must be included in `mod`

#### Error handling

- Use `cli::cli_abort()`, `cli_warn()`, `cli_inform()` throughout — `cli` is a package dependency
- `species_model()` with unsupported ID → error with closest-match suggestion (similar to usethis "did you mean X?")
- Cloud unreachable → hard error: `"Could not reach cloud host. Check your internet connection or use a previously cached species_model."`
- All constructor inputs validated at construction time (fail fast, not deferred to projection time)

### Claude's Discretion

- Exact output structure of `ipm_projection` object (beyond what the architecture doc specifies)
- S3 class naming and method conventions (print, summary, plot methods)
- Internal validation helper structure
- Naming conventions for internal-only functions

### Deferred Ideas (OUT OF SCOPE)

- lat/lon → climate lookup helper `get_climate(lat, lon)` — requires spatial climate database, deferred to future milestone
- Spatial random effects: `parameters(mod, draw = 12, random_effect = spatial_re(lat, lon))` — architecture supports it, not in v1
- Multi-species community IPM extensions beyond the current unified model — scoped as v3+
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| API-01 | All exported functions listed with finalized signatures, argument names, and return types | Architecture v2 doc + CONTEXT.md locked decisions provide the complete function table; research supplies the return-type documentation pattern |
| API-02 | `supported_species()` — returns data frame of available species IDs, names, and model variants | Locked in CONTEXT.md; research supports the tibble return type convention and documentation pattern |
| API-03 | `run_ipm(species_id, lat, lon, climate, n_draws = 1, ...)` — single entry point (SUPERSEDED by new API shape in CONTEXT.md) | CONTEXT.md explicitly replaced this with the five-constructor + two-engine design; API-03 must be interpreted as the canonical workflow documented, not a single `run_ipm` function |
| API-04 | Naming conventions documented and applied consistently across all exported functions | Research supports R package naming conventions: `snake_case` for functions, `"ipm_"` prefix for S3 classes, underscore-separated internal function names |
| API-05 | API design document produced capturing decisions, rationale, and interface contracts for Phase 2 implementation | This is the primary deliverable of the phase; research provides the template structure and what the document must contain |
</phase_requirements>

---

## Summary

Phase 1 is a **documentation phase**, not an implementation phase. Its deliverable is an API design document — a written specification that fully captures function signatures, argument contracts, return types, S3 class names, error behavior, and naming conventions. The document will guide Phase 2 implementors and serve as the source of truth when CONTEXT.md decisions need operationalization.

The API shape is almost entirely locked in CONTEXT.md and the architecture v2 document. Research confirms that the chosen design patterns (S3 classes with constructor/validator/helper separation, `cli` for error messages, fail-fast validation) align with current R package best practices from the tidyverse ecosystem. The main open area requiring judgment is Claude's Discretion items: the `ipm_projection` return structure, S3 print/summary/plot method conventions, and internal function naming.

Note on API-03: REQUIREMENTS.md defines API-03 as `run_ipm(species_id, lat, lon, climate, ...)`, but CONTEXT.md locks in a fundamentally different design (five constructors + two engines). CONTEXT.md post-dates REQUIREMENTS.md and represents the user's final decision. The API design document must document the new design as authoritative, and should explicitly note that API-03 is satisfied by the canonical workflow `stand() → species_model() → parameters() → env_condition() → control() → lambda()/project()`.

**Primary recommendation:** Produce a single `API_DESIGN.md` document in `.planning/phases/01-api-design/` that formalizes the locked decisions into a written interface contract — function tables, return types, S3 class definitions, validation rules, error messages, naming conventions — structured to be read by Phase 2 implementors without needing to look at any other file.

---

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| R (base) | >= 4.1 | S3 class system, generics | All class/method dispatch is base R S3 |
| cli | >= 3.6 | Formatted error messages, warnings, informational output | Locked decision; tidyverse standard for user-facing conditions |
| rlang | >= 1.1 | Condition signaling, call attribution in errors | cli depends on it; used by all modern R packages for abort/warn |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| stringdist | >= 0.9 | Approximate string matching for "did you mean X?" suggestions | Used in `species_model()` to suggest closest species ID on invalid input |
| tibble | >= 3.2 | Return type of `supported_species()` | Standard tidy return format; already in tidyverse dependency chain |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `cli::cli_abort()` | `stop()` / `rlang::abort()` | `cli` gives formatted, styled output with bullet points; `stop()` is plain text. Locked to `cli`. |
| S3 classes | R6 / S4 | S3 is simpler, idiomatic for modeling packages; S4 is heavier; R6 is reference semantics (mutable) which contradicts the immutability requirement for `stand` |
| `stringdist::amatch()` | Base `agrep()` | `stringdist` provides richer distance metrics; `agrep()` is simpler but less flexible for scientific species names |

**Installation (runtime deps already in project):**
```r
# These will be declared in DESCRIPTION Imports:
install.packages(c("cli", "rlang", "stringdist", "tibble"))
```

---

## Architecture Patterns

### S3 Class Pattern: Constructor / Validator / Helper

Advanced R (Hadley Wickham, https://adv-r.hadley.nz/s3.html) establishes three layers:

1. **Constructor** (`new_ipm_stand()`): Low-level, fast, minimal checks — sets class and attributes
2. **Validator** (`validate_ipm_stand()`): Thorough checks on values — throws errors with `cli_abort()`
3. **Helper** (`stand()`): User-facing function — calls constructor + validator, handles input coercion

For this package, the user-facing helper IS the exported function (e.g., `stand()`), while the low-level constructor is internal.

### Recommended Document Structure for API Design Output

The Phase 1 deliverable (`API_DESIGN.md`) should have these sections:

```
# forestIPM API Design v1

## Overview and Design Principles
## Canonical Workflow
## Exported Functions Reference
  ### Constructors
    #### stand()
    #### species_model()
    #### parameters()
    #### env_condition()
    #### control()
  ### Engine Functions
    #### lambda()
    #### project()
  ### Utility Functions
    #### supported_species()
## S3 Classes
  ### ipm_stand
  ### ipm_model
  ### ipm_parameters
  ### ipm_env
  ### ipm_control
  ### ipm_projection (Claude's Discretion)
## Naming Conventions
## Validation Rules and Error Messages
## Internal Functions (not exported)
## Design Decisions and Rationale
## Versioning and Extension Points
```

### Pattern 1: Unified Single/Multi-Species Interface
**What:** `species_model()` accepts either a character vector or a `stand` object. Single-species and multi-species models share class `"ipm_model"`. No branching at the user API level.
**When to use:** Always — this is the only model class.
**Example:**
```r
# Source: CONTEXT.md locked decisions / IPM_PACKAGE_ARCHITECTURE_v2.md
# Single species
mod <- species_model("ABIBAL")

# Multi-species
mod <- species_model(c("ABIBAL", "ACERUB"))

# Infer from stand
mod <- species_model(stand_obj)
```

### Pattern 2: Fail-Fast Validation at Construction Time
**What:** All validation happens when constructors are called, not at `lambda()` or `project()` time.
**When to use:** Every constructor must validate its inputs before returning.
**Example:**
```r
# Source: CONTEXT.md error handling locked decisions
stand <- stand(data)          # validates column names, types, positive sizes here
pars  <- parameters(mod, ...) # validates draw argument here; confirms all mod species have params
```

### Pattern 3: Named Numeric Vector Return from `lambda()`
**What:** `lambda()` returns a named numeric vector, one element per species in `mod`.
**When to use:** Always for `lambda()`.
**Example:**
```r
# Source: IPM_PACKAGE_ARCHITECTURE_v2.md
lam <- lambda(mod, pars, stand, env)
# Returns: c(ABIBAL = 1.03, ACERUB = 0.97)
```

### Pattern 4: Immutable Stand Objects
**What:** `stand` objects are never modified in-place. `project()` returns new stand states at each timestep wrapped in the projection result object.
**When to use:** Always — no exceptions.

### Anti-Patterns to Avoid
- **Deferred validation:** Do not validate `pars`/`stand` compatibility inside `project()` or `lambda()` — validate in constructors
- **Mutable state:** Do not modify `stand$trees` in-place; always return new objects
- **Exposing climate scaling:** Do not expose the [0,1] scaling parameters to users — internal only
- **`base::environment()` collision:** The function is named `env_condition()` not `environment()` — never deviate
- **Separate single/multi-species classes:** Do not create `ipm_single_model` and `ipm_community_model` — one class only

---

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| User-facing formatted errors | Custom `paste0()` error messages | `cli::cli_abort()` with bullet formatting | cli handles pluralization, inline styling, class display, call attribution |
| "Did you mean X?" species suggestion | Custom Levenshtein implementation | `stringdist::amatch()` or `agrep()` | Edge cases with scientific names (diacritics, underscores, case variants); tested library handles these |
| S3 method dispatch | `if (is(x, "ipm_stand"))` branching | `UseMethod()` + `NextMethod()` | Correct dispatch is not trivial; UseMethod handles inheritance properly |
| Return type for species table | Named list or data.frame | `tibble::tibble()` | `supported_species()` locked to return tibble; consistent with tidyverse expectations |

**Key insight:** This phase produces a document, not code. The "don't hand-roll" list informs what the API design document must specify as implementation constraints for Phase 2 — Phase 2 implementors should not be discovering these patterns themselves.

---

## Common Pitfalls

### Pitfall 1: API-03 Interpretation Conflict
**What goes wrong:** REQUIREMENTS.md API-03 specifies `run_ipm(species_id, lat, lon, climate, n_draws = 1, ...)`. CONTEXT.md locks a completely different design.
**Why it happens:** Requirements were written before the design session; CONTEXT.md reflects post-discussion decisions.
**How to avoid:** The API design document must explicitly state that the five-constructor + two-engine design supersedes the original `run_ipm` concept, and API-03 is satisfied by the canonical workflow. Reference REQUIREMENTS.md in the document's rationale section.
**Warning signs:** If Phase 2 tasks reference `run_ipm` — the API document failed to make this clear.

### Pitfall 2: Incomplete Return Type Specification
**What goes wrong:** The design document specifies what functions do but not exactly what they return — leaving Phase 2 implementors to guess the structure of `ipm_projection`.
**Why it happens:** "Exact output structure" is Claude's Discretion — easy to defer documenting it.
**How to avoid:** The document MUST specify the `ipm_projection` structure completely, including field names, types, and accessor patterns. This is the most important Claude's Discretion item.
**Warning signs:** Phase 2 asks "what does `project()` return?" — the document failed.

### Pitfall 3: `env_condition()` vs `env_conditions()` Naming Inconsistency
**What goes wrong:** CONTEXT.md section "env_conditions() naming and behavior" uses the plural form `env_conditions()`, while the "Public API surface" table uses singular `env_condition()`.
**Why it happens:** The document has a naming inconsistency between two sections.
**How to avoid:** The API design document must resolve this to ONE canonical name (CONTEXT.md public API table uses `env_condition()` — treat that as the locked form; the "naming and behavior" section header appears to be a section title, not the function name).
**Warning signs:** Any code in Phase 2 that uses `env_conditions()` plural instead of `env_condition()` singular.

### Pitfall 4: Missing Validation Contract Details
**What goes wrong:** The document says "validate at construction time" but doesn't specify exactly which conditions trigger which errors.
**Why it happens:** Validation rules are easy to leave abstract in design documents.
**How to avoid:** For each constructor, enumerate the exact validation conditions and the `cli_abort()` message template. This prevents inconsistent error messages across constructors.
**Warning signs:** Phase 2 produces different error message styles per function.

### Pitfall 5: S3 Class Name Collisions
**What goes wrong:** Using generic class names like `"model"` or `"stand"` that might collide with other packages.
**Why it happens:** Short names feel natural.
**How to avoid:** All S3 class names must use the `"ipm_"` prefix: `"ipm_stand"`, `"ipm_model"`, `"ipm_parameters"`, `"ipm_env"`, `"ipm_control"`, `"ipm_projection"`. This is per Advanced R best practice: include package name in class name.
**Warning signs:** `inherits(x, "stand")` appears in any code — should always be `"ipm_stand"`.

### Pitfall 6: `draw` Argument Type Ambiguity
**What goes wrong:** `parameters(mod, draw = 1)` — is `1` a numeric index (select draw #1 from 1-1000) or does it mean "draw 1 replicate"?
**Why it happens:** The `draw` argument accepts three forms: `"mean"`, `"random"`, and integer 1-1000. Without explicit documentation, numeric `1` is ambiguous.
**How to avoid:** The API document must specify: if `draw` is a positive integer, it selects exactly that posterior draw index. If `draw = "random"`, one draw is selected randomly. The documentation must include an example of each form.
**Warning signs:** Any question from Phase 2 about what `draw = 1` means.

---

## Code Examples

Verified patterns from official sources and locked decisions:

### Canonical Workflow (from architecture doc)
```r
# Source: IPM_PACKAGE_ARCHITECTURE_v2.md / CONTEXT.md
stand <- stand(data)
mod   <- species_model(stand)
pars  <- parameters(mod, draw = "random", seed = 42)
env   <- env_condition(MAT = 8, MAP = 1200)
ctrl  <- control(years = 100, delta_time = 1)

lambda(mod, pars, stand, env)
project(mod, pars, stand, env, ctrl)
```

### S3 Constructor Pattern (from Advanced R)
```r
# Source: https://adv-r.hadley.nz/s3.html
# Low-level constructor (internal)
new_ipm_stand <- function(trees, species, plot_size) {
  structure(
    list(
      trees = trees,
      species = species,
      plot_size = plot_size
    ),
    class = "ipm_stand"
  )
}

# Validator (internal)
validate_ipm_stand <- function(x) {
  # check required columns, positive sizes, known species
  x
}

# User-facing helper (exported)
stand <- function(data) {
  # coerce and validate
  validate_ipm_stand(
    new_ipm_stand(
      trees = data,
      species = unique(data$species),
      plot_size = data$plot_size[[1]]
    )
  )
}
```

### cli Error Message Pattern
```r
# Source: https://cli.r-lib.org/reference/cli_abort.html
# Validated in documentation 2026-03-02
cli::cli_abort(c(
  "{.arg species} must be a character vector of known species IDs.",
  "x" = "You supplied {.val {bad_sp}}.",
  "i" = "Did you mean {.val {closest_match}}?",
  "i" = "Run {.run supported_species()} to see all valid species IDs."
))
```

### `supported_species()` Return Contract
```r
# Source: CONTEXT.md locked decisions
# Returns a tibble with columns: species_id (chr), common_name (chr), model_variant (chr)
supported_species()
#> # A tibble: 93 × 3
#>   species_id  common_name          model_variant
#>   <chr>       <chr>                <chr>
#> 1 ABIBAL      Balsam fir           intcpt_plot_comp_clim
#> 2 ABIBAL      Balsam fir           intcpt_plot_comp
#> ...
```

### Time-Varying Climate in `env_condition()`
```r
# Source: CONTEXT.md locked decisions / architecture doc
# Accepts scalar or function of time
env_static  <- env_condition(MAT = 8, MAP = 1200)
env_dynamic <- env_condition(MAT = function(t) 8 + 0.02 * t, MAP = 1200)
```

### `parameters()` draw argument forms
```r
# Source: CONTEXT.md locked decisions
pars_mean   <- parameters(mod, draw = "mean")          # posterior mean
pars_random <- parameters(mod, draw = "random", seed = 42)  # random draw
pars_exact  <- parameters(mod, draw = 314)             # exact draw #314
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `stop()` for errors | `cli::cli_abort()` | ~2020 (cli 2.0) | Rich formatted messages with inline styles, class display, suggestions |
| Global `exists('fct')` check (current `init_pop()` bug BUG-02) | Proper function argument | Pre-existing bug | Must be made internal AND fixed in Phase 2 |
| `purrr::as_vector()` (deprecated, BUG-03) | `unlist()` | purrr >= 1.0.0 (2023) | Current code will warn on purrr >= 1.0.0 |
| `getPars_sp()` reading local RDS paths | `parameters(mod, ...)` fetching from cloud Parquet | Phase 3 | API design must future-proof parameter fetching; `parameters()` hides the data source |
| Single-species only | Unified single/multi-species class | New design | One class handles both; no API branching |

**Deprecated/outdated patterns that appear in current codebase:**
- `getPars_sp()`: Will become internal or removed; `parameters()` replaces the user-facing role
- `mkKernel()`: Made internal; not exported in new API
- `init_pop()`: Made internal; `stand()` and `species_model()` handle initialization
- `size_to_BAcomp()`, `dbh_to_sizeDist()`: Made internal

---

## Claude's Discretion: Recommendations

These areas were left to Claude's judgment. Research supports specific recommendations.

### `ipm_projection` Object Structure

The projection object must support time-series access. Recommended structure:

```r
# ipm_projection S3 class fields:
proj$species        # character vector of species in projection
proj$years          # numeric vector of timesteps
proj$lambda         # named list: lambda per species per timestep
proj$stand_series   # list of ipm_stand objects, one per timestep (or sparse)
proj$summary        # convenience: final-state summary tibble
```

Rationale: Researchers will want to plot dynamics over time, extract final state, compute trajectories. Storing `stand_series` enables full reconstruction without re-running. This is consistent with how `deSolve` and similar simulation packages structure output.

### S3 Method Conventions

Each S3 class should have:
- `print.<class>()`: One-liner summary (class name, key dims, first few species)
- `summary.<class>()`: Detailed structured output
- `plot.<class>()`: Optional but recommended for `ipm_projection` (trajectory plots)

Class naming (following Advanced R "use package name prefix"):
- `"ipm_stand"` — confirmed locked
- `"ipm_model"` — confirmed locked
- `"ipm_parameters"` — confirmed in architecture doc
- `"ipm_env"` — for `env_condition()` return
- `"ipm_control"` — for `control()` return
- `"ipm_projection"` — for `project()` return

### Internal Function Naming Conventions

Use `.` prefix convention for purely internal helpers (not exported, not S3 methods):
- `.validate_stand_data(data)` — validates input dataframe for `stand()`
- `.scale_climate(MAT, MAP, scaling_params)` — internal climate scaling
- `.compute_competition(stand, species)` — internal BA competition

This distinguishes internal helpers from exported functions and low-level constructors.

---

## Open Questions

1. **`env_condition()` vs `env_conditions()` naming**
   - What we know: CONTEXT.md has the plural form in one section header, singular form in the function table
   - What's unclear: Which is canonical?
   - Recommendation: Use `env_condition()` (singular) as locked — it appears in the public API table and avoids any confusion with multi-env objects. The API design document must resolve this definitively.

2. **`stand()` column handling for `plot_size`**
   - What we know: CONTEXT.md says `data` must contain `plot_size` as a column, but the architecture doc shows `plot_size = 400` as a separate argument
   - What's unclear: Is `plot_size` in the dataframe (as a column) or a separate scalar argument?
   - Recommendation: CONTEXT.md is more recent and locked — `plot_size` is a column in `data`. The API design document must specify this clearly and note it diverges from the v2 architecture doc.

3. **`species_model()` missing species handling: argument name**
   - What we know: Three modes exist (`error`, `drop`, `static`) for unsupported species; default is `error`
   - What's unclear: The argument name for this behavior was not specified in CONTEXT.md
   - Recommendation: Name the argument `on_missing = c("error", "drop", "static")` — follows R convention for option-bearing arguments. The API design document should specify this.

4. **`ipm_projection` object: how many timesteps stored**
   - What we know: `project()` runs for `ctrl$years` timesteps
   - What's unclear: Are all intermediate stand states stored (memory-intensive for long runs) or only first/last?
   - Recommendation: Store all timesteps by default for `ctrl$years <= 200`; add a `ctrl$store_every = N` parameter for long runs. Specify in API design document.

---

## Sources

### Primary (HIGH confidence)
- CONTEXT.md — all locked API decisions; direct user specification
- IPM_PACKAGE_ARCHITECTURE_v2.md — primary design reference for object responsibilities and internal structures
- https://adv-r.hadley.nz/s3.html — S3 class design: constructor/validator/helper pattern, naming conventions
- https://r-pkgs.org/man.html — roxygen2 documentation conventions for @param, @return, @examples
- https://cli.r-lib.org/reference/cli_abort.html — cli_abort/cli_warn/cli_inform signatures and usage patterns

### Secondary (MEDIUM confidence)
- REQUIREMENTS.md — requirement definitions (note: API-03 superseded by CONTEXT.md decisions)
- https://blog.r-hub.io/2023/11/30/cliff-notes-about-cli/ — cli package best practices verified against official docs
- https://cran.r-project.org/web/packages/stringdist/stringdist.pdf — stringdist package for "did you mean" matching (verified CRAN existence)

### Tertiary (LOW confidence)
- WebSearch for "did you mean" R patterns — no authoritative single implementation found; stringdist::amatch() is the standard tool but specific usage pattern for species names needs validation in Phase 2

---

## Metadata

**Confidence breakdown:**
- Locked API decisions: HIGH — directly from user in CONTEXT.md
- Architecture patterns: HIGH — verified against official Advanced R documentation
- Error messaging: HIGH — verified against cli official docs
- `ipm_projection` structure: MEDIUM — Claude's Discretion recommendation; plausible but not locked
- `on_missing` argument name: MEDIUM — follows convention but name not specified in locked decisions
- "did you mean" implementation: LOW — stringdist is the right tool but specific implementation not validated

**Research date:** 2026-03-02
**Valid until:** 2026-04-02 (stable domain — R ecosystem patterns change slowly)
