# forestIPM API Design v1

**Status:** Final — Phase 1 deliverable
**Date:** 2026-03-02
**Audience:** Phase 2 implementors — read this document alone; no need to open CONTEXT.md, RESEARCH.md, or the architecture doc.

---

## Table of Contents

1. [Overview and Design Principles](#1-overview-and-design-principles)
2. [Canonical Workflow](#2-canonical-workflow)
3. [Exported Functions Reference](#3-exported-functions-reference)
   - 3.1 [stand()](#31-stand)
   - 3.2 [species_model()](#32-species_model)
   - 3.3 [parameters()](#33-parameters)
   - 3.4 [env_condition()](#34-env_condition)
   - 3.5 [control()](#35-control)
   - 3.6 [lambda()](#36-lambda)
   - 3.7 [project()](#37-project)
   - 3.8 [supported_species()](#38-supported_species)
4. [S3 Classes](#4-s3-classes)
5. [Naming Conventions](#5-naming-conventions)
6. [Constructor / Validator / Helper Pattern](#6-constructor--validator--helper-pattern)
7. [Validation Rules and Error Messages](#7-validation-rules-and-error-messages)
8. [Internal Functions (not exported)](#8-internal-functions-not-exported)
9. [Design Decisions and Rationale](#9-design-decisions-and-rationale)
10. [Versioning and Extension Points](#10-versioning-and-extension-points)

---

## 1. Overview and Design Principles

### Package Description

forestIPM is *"A stochastic dynamic forest Integral Projection Model built on Bayesian species-specific demographic rates with endogenous competition and flexible climate drivers."*

The package computes population growth rates (lambda) and projects size-structured tree population dynamics over time for one or more species in eastern North America. Parameters are drawn from species-specific Bayesian posterior distributions. Competition is endogenous — basal area updated from the modelled population dynamics each timestep. Climate drivers can be static or time-varying.

### Design Principles

1. **Fail-fast validation at construction time.** All input checking happens in constructors — `stand()`, `species_model()`, `parameters()`, `env_condition()`, `control()` — not inside engines. By the time `lambda()` or `project()` is called, all objects are guaranteed valid. Errors point to the exact source.

2. **Immutable stand objects.** A `stand` object is never modified in-place by any function. `project()` returns new stand states as part of the projection result; the original stand is never touched.

3. **Stochasticity fully resolved in `parameters()`.** All randomness lives in one place. Calling `parameters(mod, draw = "random", seed = 42)` resolves a single parameter realization; everything downstream (`lambda()`, `project()`) is deterministic given those parameters. No hidden randomness inside projection.

4. **Single unified class for single- and multi-species models.** `species_model()` returns `"ipm_spModel"` regardless of whether one or twenty species are included. There is no `"ipm_single_model"` or `"ipm_community_model"` class — these are anti-patterns. Multi-species is a special case handled by the same constructors and engines.

5. **No wrapper functions — researchers compose constructors and engines directly.** The five constructors (`stand`, `species_model`, `parameters`, `env_condition`, `control`) and two engines (`lambda`, `project`) are the entire public API surface. There is no convenience wrapper like `run_ipm()`. Researchers compose them explicitly to keep stochasticity, parameters, and control visible in their scripts.

---

## 2. Canonical Workflow

The complete forestIPM workflow from raw inventory data to population dynamics:

```r
library(forestIPM)

# 1. Represent the forest plot
stand <- stand(data)                              # data is a data.frame with tree records

# 2. Define which species to model (inferred from the stand, or explicit vector)
mod   <- species_model(stand)                     # or: species_model(c("ABIBAL", "ACERUB"))

# 3. Resolve one parameter realization from Bayesian posteriors
pars  <- parameters(mod, draw = "random", seed = 42)

# 4. Set climate drivers
env   <- env_condition(MAT = 8, MAP = 1200)       # static climate; see 3.4 for time-varying

# 5. Configure projection
ctrl  <- control(years = 100, delta_time = 1)

# 6. Compute dominant eigenvalue (asymptotic growth rate per species)
lambda(mod, pars, stand, env)
#> <ipm_lambda>  ABIBAL: 1.03  ACERUB: 0.97

# 7. Project population dynamics over time
proj  <- project(mod, pars, stand, env, ctrl)
```

> **API-03 Note:** The original requirements document (REQUIREMENTS.md) specified a single `run_ipm(species_id, lat, lon, climate, n_draws = 1, ...)` entry point. The design session replaced this with the five-constructor + two-engine pattern above. **This canonical workflow IS the fulfillment of API-03.** The composable design gives researchers explicit control over stochasticity (`parameters()`), stand state, and simulation control that a single wrapper function cannot provide.

---

## 3. Exported Functions Reference

### 3.1 `stand()`

**Signature:** `stand(data)`

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `data` | data.frame or tibble | required | Individual-level tree records. Must contain a size column, a species column, and a `plot_size` column (see validation rules). |

#### Required Columns in `data`

| Column | Accepted Names | Type | Constraint |
|--------|---------------|------|-----------|
| Size | `size` OR `dbh` | numeric | > 0 (mm). Either column name is accepted; both in the same data.frame causes an error. |
| Species | `sp` OR `species` OR `species_id` | character | Must match IDs in `supported_species()` (validated later by `species_model()`). |
| Plot area | `plot_size` | numeric scalar | > 0 (m2). Must be the same value for all rows (one plot per call to `stand()`). |

> **Note:** `plot_size` is a column in the data.frame, NOT a separate argument to `stand()`. This is a deliberate design decision — the column keeps plot metadata co-located with tree records.

#### Return Value

Object of S3 class `"ipm_stand"` with three fields:

| Field | Type | Description |
|-------|------|-------------|
| `$trees` | data.frame | Individual-level records with standardized column names: `size_mm` (numeric), `species_id` (character). Original column names are normalized internally. |
| `$species` | character vector | Unique species IDs present in the stand. |
| `$plot_size` | numeric scalar | Plot area in m2, extracted from `data$plot_size`. |

#### Validation Rules

All checks run at construction time via `cli::cli_abort()`:

| Condition | Error message |
|-----------|--------------|
| No `size` or `dbh` column | `"data must contain a column named {.field size} or {.field dbh} (size in mm)."` |
| `size` or `dbh` values <= 0 | `"All tree sizes must be positive (> 0 mm). Found {n} non-positive values."` |
| No species column (`sp`, `species`, or `species_id`) | `"data must contain a column named {.field sp}, {.field species} or {.field species_id}."` |
| No `plot_size` column | `"data must contain a {.field plot_size} column (plot area in m2)."` |

#### Usage Example

```r
# Minimal data.frame
data <- data.frame(
  dbh       = c(12.5, 34.0, 8.2),        # diameter at breast height in mm
  species_id = c("ABIBAL", "ABIBAL", "ACERUB"),
  plot_size = 400                          # 400 m2 plot
)

s <- stand(data)
s$trees      # data.frame with size_mm and species_id
s$species    # c("ABIBAL", "ACERUB")
s$plot_size  # 400
```

---

### 3.2 `species_model()`

**Signature:** `species_model(x, on_missing = "error")`

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `x` | character vector OR `ipm_stand` | required | Species to model. If `ipm_stand`, species are inferred from `x$species`. If character vector, species are used directly. |
| `on_missing` | character | `"error"` | Behavior when a species ID is not in `supported_species()`. One of `"error"`, `"drop"`, or `"static"`. See below. |

#### `on_missing` Behavior

| Value | Behavior |
|-------|----------|
| `"error"` (default) | `cli_abort()` with the unknown IDs and a closest-match suggestion via `stringdist::amatch()`. |
| `"drop"` | Silently remove the unsupported species from the model. If `x` was an `ipm_stand`, the dropped species are also excluded from the stand used in subsequent computations. |
| `"static"` | Keep the unsupported species for competition background estimation only. They do NOT contribute to `lambda()` or `project()` output, and their trees stay static (do not grow, survive, or recruit) over time. |

> **Implementation note:** Use `stringdist::amatch()` for closest-match suggestions, not a hand-rolled Levenshtein distance.

#### Return Value

Object of S3 class `"ipm_spModel"` with fields:

| Field | Type | Description |
|-------|------|-------------|
| `$species` | character vector | Species IDs included in the model (after applying `on_missing` resolution). |
| `$on_missing` | character | The `on_missing` value used at construction. |

#### Validation Rules

| Condition | Error (when `on_missing = "error"`) |
|-----------|-------------------------------------|
| One or more species IDs not in `supported_species()` | See error template below. |

**Error template for unknown species (`on_missing = "error"`):**

```r
cli::cli_abort(c(
  "{.arg x} contains species IDs not found in {.run supported_species()}.",
  "x" = "Unknown: {.val {bad_ids}}.",
  "i" = "Did you mean {.val {closest_match}}?",
  "i" = "Run {.run supported_species()} to see all valid IDs."
))
```

#### Usage Examples

```r
# From a stand object (species inferred)
mod <- species_model(stand_obj)

# From an explicit species vector
mod <- species_model(c("ABIBAL", "ACERUB"))

# Unsupported species dropped silently
mod <- species_model(c("ABIBAL", "FAKEID"), on_missing = "drop")
# mod$species => c("ABIBAL")  # FAKEID removed

# Unsupported species kept as static background
mod <- species_model(c("ABIBAL", "FAKEID"), on_missing = "static")
# mod$species => c("ABIBAL", "FAKEID")
# FAKEID trees contribute to competition but are excluded from lambda()/project() output
```

---

### 3.3 `parameters()`

**Signature:** `parameters(mod, draw = "random", seed = NULL)`

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `mod` | `ipm_spModel` | required | Species model object defining which species need parameters. |
| `draw` | character or positive integer | `"random"` | How to draw from the posterior. See draw forms below. |
| `seed` | integer or NULL | `NULL` | Random seed for reproducibility. Used only when `draw = "random"`. Ignored for `"mean"` or integer draw. |

#### `draw` Argument Forms

| Form | Type | Meaning |
|------|------|---------|
| `"mean"` | character | Use the posterior mean for all parameters. No randomness. |
| `"random"` | character | Sample one draw uniformly at random from the 1-1000 posterior draws. Use `seed` for reproducibility. |
| Integer 1–1000 | positive integer | Select exactly that posterior draw index. `draw = 314` means draw #314 (the 314th stored draw), NOT "draw 314 replicates". |

> **Important:** `draw` as a positive integer is a **draw index**, not a count. `parameters(mod, draw = 1)` selects the first stored posterior draw, not one random replicate.

#### Usage Examples

```r
# All three draw forms must be understood by Phase 2 implementors:

pars_mean   <- parameters(mod, draw = "mean")               # posterior mean — no randomness
pars_random <- parameters(mod, draw = "random", seed = 42)  # reproducible random draw
pars_exact  <- parameters(mod, draw = 314)                  # exact draw #314 from posterior
```

#### Return Value

Object of S3 class `"ipm_parameters"` with fields:

| Field | Type | Description |
|-------|------|-------------|
| `$species_params` | named list | Keys are `species_id` strings. Each element contains `$fixed` (fixed-effect parameters), `$random_effects` (plot-level random effects), and `$draw_id` (the integer index of the draw used). |
| `$draw_type` | character | The resolved draw specification: `"mean"`, `"random"`, or the integer index as a string. |
| `$seed` | integer or NULL | Seed used for the draw. `NULL` unless `draw = "random"`. |

**Internal structure of `$species_params`:**

```r
pars$species_params$ABIBAL$fixed         # named numeric vector of fixed-effect parameter values
pars$species_params$ABIBAL$random_effects # named numeric vector of plot random effects
pars$species_params$ABIBAL$draw_id        # integer: which posterior draw was selected
```

#### Validation Rules

| Condition | Error message |
|-----------|--------------|
| `mod` is not an `ipm_spModel` | `"{.arg mod} must be an object of class {.cls ipm_spModel}. Got {.cls {class(mod)}}."` |
| `draw` is not `"mean"`, `"random"`, or a positive integer | `"{.arg draw} must be {.val \"mean\"}, {.val \"random\"}, or a positive integer (1-1000). Got {.val {draw}}."` |
| `draw` is an integer outside 1–1000 | `"{.arg draw} integer index must be between 1 and 1000. Got {draw}."` |
| Cloud unreachable during parameter fetch | `"Could not reach cloud host. Check your internet connection or use a previously cached species_model."` |

---

### 3.4 `env_condition()`

> **Canonical name:** `env_condition()` — singular. **NEVER use `env_conditions()` (plural).** The plural form was an internal naming inconsistency in earlier drafts; it is not the API name. The function creates one environment condition object.

**Signature:** `env_condition(MAT, MAP)`

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `MAT` | numeric scalar OR `function(t)` | required | Mean Annual Temperature in Celsius. Scalar for static climate; function of time for dynamic (time-varying) climate. |
| `MAP` | numeric scalar OR `function(t)` | required | Mean Annual Precipitation in mm/year. Scalar for static; function of time for dynamic. |

> **Climate scaling:** MAT and MAP are scaled to [0, 1] internally using stored scaling parameters before being passed to demographic models. The scaling parameters are **never exposed** in the returned object or to the user. This is intentional — users work in natural units (°C and mm/year); the package handles normalization.

#### Return Value

Object of S3 class `"ipm_env"` with fields:

| Field | Type | Description |
|-------|------|-------------|
| `$MAT` | numeric or function | The MAT value or function as provided by the user (unscaled). |
| `$MAP` | numeric or function | The MAP value or function as provided by the user (unscaled). |

#### Validation Rules

| Condition | Error message |
|-----------|--------------|
| `MAT` is not a numeric scalar or a function | `"{.arg MAT} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(MAT)}}."` |
| `MAP` is not a numeric scalar or a function | `"{.arg MAP} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(MAP)}}."` |

#### Usage Examples

```r
# Static climate
env_static  <- env_condition(MAT = 8, MAP = 1200)

# Dynamic (time-varying) climate — warming scenario
env_dynamic <- env_condition(MAT = function(t) 8 + 0.02 * t, MAP = 1200)

# Both time-varying
env_both    <- env_condition(
  MAT = function(t) 8 + 0.02 * t,
  MAP = function(t) 1200 - 2 * t
)
```

---

### 3.5 `control()`

**Signature:** `control(years = 100, delta_time = 1, store_every = 1, bin_width = 1)`

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `years` | positive integer | `100` | Total number of simulation timesteps. |
| `delta_time` | positive numeric | `1` | Duration of each timestep in years. Default 1 = annual timesteps. |
| `store_every` | positive integer | `1` | Store stand state every N timesteps. Default 1 = store every year. Use `store_every = 10` for runs > 200 years to reduce memory footprint. |
| `bin_width` | positive integer | `1` | Bin width used in the IPM kernel discretization over the size axis. |

#### Return Value

Object of S3 class `"ipm_control"` with fields corresponding to the arguments: `$years`, `$delta_time`, `$store_every`, `$bin_width`.

#### Validation Rules

| Condition | Error message |
|-----------|--------------|
| `years` not a positive integer | `"{.arg years} must be a positive integer. Got {.val {years}}."` |
| `delta_time` not a positive numeric | `"{.arg delta_time} must be a positive numeric. Got {.val {delta_time}}."` |
| `store_every` not a positive integer | `"{.arg store_every} must be a positive integer. Got {.val {store_every}}."` |
| `bin_width` not a positive integer | `"{.arg bin_width} must be a positive integer. Got {.val {bin_width}}."` |

#### Usage Examples

```r
ctrl <- control(years = 100, delta_time = 1)          # annual, store all years
ctrl <- control(years = 500, store_every = 10)         # 500-year run, store every 10th year
ctrl <- control(years = 200, bin_width = 5)            # wider size bins for faster kernel assembly
```

---

### 3.6 `lambda()`

**Signature:** `lambda(mod, pars, stand, env, ctrl)`

#### Arguments

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `mod` | `ipm_spModel` | yes | Defines which species to compute lambda for. |
| `pars` | `ipm_parameters` | yes | Parameter realization; must contain parameters for all species in `mod`. |
| `stand` | `ipm_stand` | yes | Initial size distribution and competition context (basal area calculation). |
| `env` | `ipm_env` | yes | Climate drivers. |
| `ctrl` | `ipm_control` | no (optional) | Simulation configuration. `years` and `store_every` are not used by `lambda()` — only `bin_width` and `delta_time` apply. May be omitted to use defaults. |

#### Return Value

Object of S3 class `"ipm_lambda"`, which is a **named numeric vector** — one element per species in `mod`.

```r
lam <- lambda(mod, pars, stand, env)
lam
#> <ipm_lambda>
#>  ABIBAL  ACERUB
#>    1.03    0.97
```

Access individual lambdas: `lam["ABIBAL"]` or `lam[["ABIBAL"]]`.

#### Competition Handling

`lambda()` does **NOT** require parameters for heterospecific species. Background competition is computed from the stand structure only — `stand$trees` provides the basal area of all species regardless of whether they are in `mod`. No dynamic coupling is required to compute the dominant eigenvalue.

#### Focal Species Efficiency Note

By default, `parameters(mod, ...)` will resolve parameters for **all species in `mod`**, and `lambda()` will compute lambda for all of them. If you want lambda for only a focal species, limit the model to that species:

```r
# Efficient: compute lambda only for ABIBAL
mod_focal <- species_model("ABIBAL")
pars_focal <- parameters(mod_focal, draw = "random", seed = 42)
lambda(mod_focal, pars_focal, stand, env)

# Avoid: computing lambda for 10 species when you only need 1 is expensive
```

This matters because each `lambda()` call involves a full eigenvalue computation per species.

#### Validation Rules

| Condition | Error message |
|-----------|--------------|
| `pars` does not cover all species in `mod` | `"pars does not contain parameters for all species in mod. Missing: {.val {missing_ids}}."` |
| `stand` contains species not in `mod` | `"stand contains species not in mod: {.val {extra_ids}}."` |

---

### 3.7 `project()`

**Signature:** `project(mod, pars, stand, env, ctrl)`

#### Arguments

Same as `lambda()`, with `ctrl` required:

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `mod` | `ipm_spModel` | yes | Defines which species to project. |
| `pars` | `ipm_parameters` | yes | Parameter realization for all species in `mod`. |
| `stand` | `ipm_stand` | yes | Initial stand state (t = 0). |
| `env` | `ipm_env` | yes | Climate drivers (may be time-varying functions). |
| `ctrl` | `ipm_control` | yes | Simulation configuration: years, timestep, storage frequency. |

#### Return Value

Object of S3 class `"ipm_projection"`. See [Section 4.6](#46-ipm_projection) for the complete field specification.

#### Multi-species vs. Single-species Behavior

| Case | Behavior |
|------|----------|
| Multiple species in `mod` | **Fully coupled dynamics.** Competition (basal area) is recomputed from the projected size distribution at each timestep. Species interact. |
| Single species in `mod` | **Independent dynamic.** Background competition from the initial stand; no within-projection coupling. |

#### Validation Rules

Same as `lambda()`:

| Condition | Error message |
|-----------|--------------|
| `pars` does not cover all species in `mod` | `"pars does not contain parameters for all species in mod. Missing: {.val {missing_ids}}."` |
| `stand` contains species not in `mod` | `"stand contains species not in mod: {.val {extra_ids}}."` |

---

### 3.8 `supported_species()`

**Signature:** `supported_species()`

#### Arguments

None.

#### Return Value

A `tibble` with exactly six columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `species_id` | character | Species ID used throughout the API | `"ABIBAL"` |
| `common_name` | character | English common name | `"Balsam fir"` |
| `nom_commun` | character | French common name | `"Sapin baumier"` |
| `growth_model` | character | Growth model variant string | `"intcpt_plot_comp_clim"` |
| `surv_model` | character | Survival model variant string | `"intcpt_plot_comp_clim"` |
| `recruit_model` | character | Recruitment model variant string | `"intcpt_plot_comp_clim"` |

One row per supported species.

```r
supported_species()
#> # A tibble: 93 × 6
#>   species_id  common_name    nom_commun      growth_model             surv_model               recruit_model
#>   <chr>       <chr>          <chr>           <chr>                    <chr>                    <chr>
#> 1 ABIBAL      Balsam fir     Sapin baumier   intcpt_plot_comp_clim    intcpt_plot_comp_clim    intcpt_plot_comp_clim
#> 2 ACERUB      Red maple      Erable rouge    intcpt_plot_comp_clim    intcpt_plot_comp_clim    intcpt_plot_comp_clim
#> ...
```

---

## 4. S3 Classes

### 4.1 `ipm_stand`

- **Class string:** `"ipm_stand"`
- **Constructor pattern:** See [Section 6](#6-constructor--validator--helper-pattern) for the three-layer pattern.
- **Fields:** `$trees` (data.frame), `$species` (character vector), `$plot_size` (numeric scalar)
- **print.ipm_stand():** One-liner: `<ipm_stand>  N trees | K species | P m2 plot`
- **summary.ipm_stand():** Table of species with tree counts and size statistics.

### 4.2 `ipm_spModel`

- **Class string:** `"ipm_spModel"` — this is the ONLY model class. NEVER use `"ipm_single_model"` or `"ipm_community_model"`.
- **Fields:** `$species` (character vector), `$on_missing` (character)
- **print.ipm_spModel():** `<ipm_spModel>  K species: ABIBAL, ACERUB, ...`

### 4.3 `ipm_parameters`

- **Class string:** `"ipm_parameters"`
- **Fields:** `$species_params` (named list), `$draw_type` (character), `$seed` (integer or NULL)
- **print.ipm_parameters():** `<ipm_parameters>  draw=random seed=42 | K species`

### 4.4 `ipm_env`

- **Class string:** `"ipm_env"` — returned by `env_condition()` (singular, always)
- **Fields:** `$MAT` (numeric or function), `$MAP` (numeric or function)
- **print.ipm_env():** `<ipm_env>  MAT=8°C  MAP=1200mm/yr` (or `MAT=<function>` for dynamic)

### 4.5 `ipm_control`

- **Class string:** `"ipm_control"`
- **Fields:** `$years` (integer), `$delta_time` (numeric), `$store_every` (integer), `$bin_width` (integer)
- **print.ipm_control():** `<ipm_control>  100 years | dt=1 | store_every=1 | bin_width=1`

### 4.6 `ipm_projection`

- **Class string:** `"ipm_projection"`
- **Complete field specification** (all fields required — no TBD):

| Field | Type | Length / Structure | Description |
|-------|------|-------------------|-------------|
| `$species` | character vector | Length = number of species in `mod` | Species IDs included in the projection. |
| `$years` | numeric vector | Length = `ceiling(ctrl$years / ctrl$store_every)` | Timestep indices at which stand state was stored (e.g., 1, 2, ..., 100 for `store_every = 1`; 10, 20, ..., 100 for `store_every = 10`). |
| `$lambda` | named list of numeric vectors | One list element per species; each vector has length = `length($years)` | Per-species lambda computed at each stored timestep. Access: `proj$lambda$ABIBAL`. |
| `$stand_series` | list of `ipm_stand` objects | Length = `length($years)` | Full stand state at each stored timestep. `proj$stand_series[[1]]` is the stand at the first stored timestep. Supports full reconstruction without re-running. |
| `$summary` | tibble | Rows = `length($years) × length($species)` | Convenience summary. Columns: `timestep` (int), `species_id` (chr), `lambda` (dbl), `n_trees` (int). |

**Access patterns:**

```r
proj <- project(mod, pars, stand, env, ctrl)

proj$species                    # c("ABIBAL", "ACERUB")
proj$years                      # 1:100 (or sparse if store_every > 1)
proj$lambda$ABIBAL              # numeric vector, lambda at each stored year
proj$stand_series[[50]]         # ipm_stand at year 50
proj$summary                    # tibble for plotting/analysis
```

- **print.ipm_projection():** `<ipm_projection>  K species | T years | store_every=N`
- **summary.ipm_projection():** Table showing final-timestep lambda and tree counts per species.
- **plot.ipm_projection():** Trajectory plot of lambda and/or tree count over time, one line per species. Uses `ctrl$store_every` on the x-axis.

### S3 Method Summary

Each class must implement:

| Generic | Classes implementing it |
|---------|------------------------|
| `print.<class>()` | All six classes |
| `summary.<class>()` | All six classes |
| `plot.<class>()` | `ipm_projection` (required); others optional |

---

## 5. Naming Conventions

| Category | Convention | Examples |
|----------|-----------|---------|
| Exported functions | `snake_case` | `stand`, `species_model`, `parameters`, `env_condition`, `control`, `lambda`, `project`, `supported_species` |
| S3 class names | `"ipm_"` prefix + noun | `"ipm_stand"`, `"ipm_spModel"`, `"ipm_parameters"`, `"ipm_env"`, `"ipm_control"`, `"ipm_projection"` |
| S3 methods | `<generic>.<class>` | `print.ipm_stand`, `summary.ipm_spModel`, `plot.ipm_projection` |
| Low-level constructors (internal) | `new_<class>` | `new_ipm_stand`, `new_ipm_spModel`, `new_ipm_parameters` |
| Validators (internal) | `validate_<class>` | `validate_ipm_stand`, `validate_ipm_spModel` |
| Internal helper functions | `.` prefix | `.scale_climate`, `.compute_competition`, `.validate_stand_data` |
| Type checks in code | `inherits(x, "ipm_stand")` | ALWAYS use the full class string — never `inherits(x, "stand")` |

### Anti-Pattern: Class Name Collisions

**Never** use short class names without the `ipm_` prefix:
- Wrong: `class(x) == "stand"`, `inherits(x, "model")`, `"spModel"`
- Correct: `inherits(x, "ipm_stand")`, `inherits(x, "ipm_spModel")`

---

## 6. Constructor / Validator / Helper Pattern

All exported constructors MUST follow the three-layer S3 pattern from Advanced R (Wickham, https://adv-r.hadley.nz/s3.html):

1. **`new_<class>(...)`** — Low-level, fast, minimal checks. Sets the class attribute and required attributes. **INTERNAL — never exported.**
2. **`validate_<class>(x)`** — Thorough value-level checks. Throws `cli::cli_abort()` on any failure. Returns `x` unchanged if valid. **INTERNAL — never exported.**
3. **`<helper>(...)`** — User-facing function. Handles input coercion (e.g., normalizing column names), calls the constructor, then the validator. **EXPORTED.**

### Example: `ipm_stand`

```r
# Low-level constructor (internal)
new_ipm_stand <- function(trees, species, plot_size) {
  structure(
    list(trees = trees, species = species, plot_size = plot_size),
    class = "ipm_stand"
  )
}

# Validator (internal)
validate_ipm_stand <- function(x) {
  # Check that $trees has required columns, sizes are positive, etc.
  # Each failure throws cli_abort() — see Section 7
  x
}

# User-facing helper (exported)
stand <- function(data) {
  validate_ipm_stand(
    new_ipm_stand(
      trees     = data,
      species   = unique(data$species_id),
      plot_size = data$plot_size[[1]]
    )
  )
}
```

This pattern is applied identically to `ipm_spModel`, `ipm_parameters`, `ipm_env`, `ipm_control`, and `ipm_projection`.

---

## 7. Validation Rules and Error Messages

All `cli::cli_abort()` calls use the named vector format:

```r
cli::cli_abort(c(
  "Primary message stating what failed.",
  "x" = "What went wrong (specific value or condition).",
  "i" = "How to fix it."
))
```

**Inline formatting tokens:**
- `{.arg arg_name}` — argument name (e.g., `{.arg data}`)
- `{.field col_name}` — column or field name (e.g., `{.field plot_size}`)
- `{.val value}` — a literal value (e.g., `{.val "FAKEID"}`)
- `{.run fn()}` — a runnable expression (e.g., `{.run supported_species()}`)
- `{.cls class_name}` — a class name (e.g., `{.cls ipm_stand}`)
- `{.code expression}` — inline code (e.g., `{.code function(t) ...}`)

### Complete Validation Table

#### `stand()`

| Condition | Message |
|-----------|---------|
| No `size` or `dbh` column | `"data must contain a column named {.field size} or {.field dbh} (size in mm)."` |
| `size` or `dbh` <= 0 | `"All tree sizes must be positive (> 0 mm). Found {n} non-positive values."` |
| No species column | `"data must contain a column named {.field sp}, {.field species} or {.field species_id}."` |
| No `plot_size` column | `"data must contain a {.field plot_size} column (plot area in m2)."` |

#### `species_model()` with `on_missing = "error"`

| Condition | Message |
|-----------|---------|
| Unknown species IDs | ```cli_abort(c("{.arg x} contains species IDs not found in {.run supported_species()}.", "x" = "Unknown: {.val {bad_ids}}.", "i" = "Did you mean {.val {closest_match}}?", "i" = "Run {.run supported_species()} to see all valid IDs."))``` |

#### `parameters()`

| Condition | Message |
|-----------|---------|
| `mod` not `ipm_spModel` | `"{.arg mod} must be an object of class {.cls ipm_spModel}. Got {.cls {class(mod)}}."` |
| `draw` not `"mean"`, `"random"`, or positive integer | `"{.arg draw} must be {.val \"mean\"}, {.val \"random\"}, or a positive integer (1-1000). Got {.val {draw}}."` |
| `draw` integer outside 1–1000 | `"{.arg draw} integer index must be between 1 and 1000. Got {draw}."` |
| Cloud unreachable | `"Could not reach cloud host. Check your internet connection or use a previously cached species_model."` |

#### `env_condition()`

| Condition | Message |
|-----------|---------|
| `MAT` not numeric scalar or function | `"{.arg MAT} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(MAT)}}."` |
| `MAP` not numeric scalar or function | `"{.arg MAP} must be a numeric scalar or a function of time {.code function(t) ...}. Got {.cls {class(MAP)}}."` |

#### `control()`

| Condition | Message |
|-----------|---------|
| `years` not positive integer | `"{.arg years} must be a positive integer. Got {.val {years}}."` |
| `delta_time` not positive numeric | `"{.arg delta_time} must be a positive numeric. Got {.val {delta_time}}."` |
| `store_every` not positive integer | `"{.arg store_every} must be a positive integer. Got {.val {store_every}}."` |
| `bin_width` not positive integer | `"{.arg bin_width} must be a positive integer. Got {.val {bin_width}}."` |

#### Cross-constructor Rules at Engine Call (`lambda()` / `project()`)

| Condition | Message |
|-----------|---------|
| `pars` missing species in `mod` | `"pars does not contain parameters for all species in mod. Missing: {.val {missing_ids}}."` |
| `stand` has species not in `mod` | `"stand contains species not in mod: {.val {extra_ids}}."` |

---

## 8. Internal Functions (not exported)

The following functions exist in the current codebase and **MUST NOT be exported** in Phase 2 NAMESPACE. They are internal implementation details hidden from package users.

| Function | Role | Replacement |
|----------|------|-------------|
| `mkKernel()` | Assembles IPM kernel (P and F matrices) | Called internally by `lambda()` and `project()` |
| `init_pop()` | Generates initial population size distribution | Called internally during stand initialization |
| `getPars_sp()` | Reads species parameters from local RDS paths | Role replaced at user level by `parameters()`; may be used internally during fetch |
| `size_to_BAcomp()` | Computes competition basal area from size distribution | Called internally for competition calculations |
| `dbh_to_sizeDist()` | Converts DBH to size distribution | Called internally during `stand()` processing |

These functions should:
- Be in unexported `.R` files, OR use the `.` prefix convention
- NOT appear in `NAMESPACE` (no `@export` roxygen tag)
- Be accessible within the package via `:::` for testing only

---

## 9. Design Decisions and Rationale

### Decision 1: Five constructors + two engines instead of `run_ipm()`

**Decision:** Replace the originally specified `run_ipm(species_id, lat, lon, climate, n_draws = 1, ...)` with five composable constructors and two engine functions.

**Rationale:**
- Researchers can call `lambda()` without running a full projection — avoids expensive time-series computation when only the asymptotic growth rate is needed.
- `parameters()` cleanly separates stochasticity from deterministic computation — all randomness is explicit and controllable in one place.
- Each constructor is independently testable and debuggable.
- Composability: researchers can reuse `stand` and `env` across multiple `parameters()` draws without reconstructing objects.

**API-03 supersession:** REQUIREMENTS.md API-03 specified `run_ipm(species_id, lat, lon, climate, n_draws = 1, ...)`. The canonical workflow in [Section 2](#2-canonical-workflow) supersedes and fulfills API-03. Phase 2 implementors MUST NOT implement a `run_ipm()` function.

### Decision 2: One `ipm_spModel` class for single and multi-species

**Decision:** `species_model()` always returns `"ipm_spModel"` regardless of the number of species.

**Rationale:** Avoids API branching where `lambda()` and `project()` would need to dispatch on model type. Multi-species is a special case of the same constructors and engines, not a separate concept. This follows the tidyverse principle of minimal surprising behavior.

### Decision 3: Fail-fast validation at construction time

**Decision:** All input validation runs when constructors are called, not deferred to `lambda()` or `project()`.

**Rationale:** Errors with a line number pointing to `stand(data)` are easier to debug than errors thrown deep inside `project()` after 50 years of simulation. Engines receive valid objects and can trust them — no defensive coding inside engines.

### Decision 4: S3 over R6 / S4

**Decision:** Use base R S3 class system throughout.

**Rationale:** S3 is simpler and idiomatic for ecological modeling packages. S4 is heavier and requires formal class declarations. R6 uses reference semantics (mutable) which directly contradicts the immutability requirement for `stand` objects. S3 is consistent with the tidyverse ecosystem that researchers in this domain already use.

### Decision 5: `env_condition()` singular — canonical and final

**Decision:** The function is named `env_condition()` (singular). `env_conditions()` (plural) is NEVER used.

**Rationale:** The function creates one environment condition object — a scalar or time-varying driver pair. The plural form implies a list of multiple conditions, which is not what this function creates. The singular form appears in the public API table in CONTEXT.md and is treated as locked.

### Decision 6: `plot_size` as a column in `data`, not a separate argument

**Decision:** `plot_size` is a column in the `data` data.frame passed to `stand()`, not a separate scalar argument.

**Rationale:** CONTEXT.md (the more recent locked-decision document) specifies `plot_size` as a column, while the earlier architecture v2 doc showed it as a separate argument. The column approach keeps plot metadata co-located with tree records and is more robust when processing multiple plots.

### Decision 7: `store_every` for memory management in long runs

**Decision:** `control()` includes a `store_every = N` argument. Default 1 stores all timesteps; larger values reduce memory for long runs.

**Rationale:** 500-year runs with 1000 size classes × 10 species × annual storage can exceed memory limits. `store_every = 10` provides a 10x memory reduction. This is a performance control exposed to researchers who run long simulations.

---

## 10. Versioning and Extension Points

The following extension points are **designed into the API** but are **NOT implemented in v1**. Phase 2 implementors should be aware of them to avoid blocking future additions.

| Extension | Target version | API surface | Implementation notes |
|-----------|---------------|-------------|---------------------|
| Spatial random effects | v2 | `parameters(mod, draw = 12, random_effect = spatial_re(lat, lon))` | `spatial_re()` would be a new constructor; `parameters()` signature has a `...` slot for future extension args |
| Climate lookup helper | Future | `get_climate(lat, lon)` — returns an `ipm_env` object | Requires spatial climate database (AWS or cloud tile server); deferred |
| Multi-species community extensions | v3+ | Out of scope for v1 API | Current unified model already supports multi-species; v3 refers to coexistence analysis tools built on top |
| ML decoder for novel plots | v3+ | Internal to `parameters()` | Extrapolation for locations outside training data; parameters() hides the data source so this is backward-compatible |

---

*Document version: 1.0*
*Phase: 01-api-design*
*Produced: 2026-03-02*
*Self-contained: Phase 2 implementors need no other file to implement the package API.*
