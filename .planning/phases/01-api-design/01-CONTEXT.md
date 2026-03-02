# Phase 1: API Design - Context

**Gathered:** 2026-03-02
**Status:** Ready for planning

<domain>
## Phase Boundary

Define the public API surface of the forestIPM package — exported function signatures, argument names, return types, and user-facing conventions — before any implementation begins. Output is an API design document that guides Phase 2 implementation.

</domain>

<decisions>
## Implementation Decisions

### Overall API shape

- **Minimal object types**: five core objects + two engine functions; single- and multi-species unified under one class
- Strictly separate: ecological state (`stand`), structural model (`species_model`), parameter realization (`parameters`), climate drivers (`env_condition`), simulation control (`control`)
- No separate `community_params` class — multi-species is just a special case of the same types
- No unnecessary wrapper functions; researchers compose constructors + engine functions directly
- Existing low-level functions (`mkKernel`, `init_pop`, `getPars_sp`, `size_to_BAcomp`, `dbh_to_sizeDist`) made internal

### Public API surface (final)

**Constructors:**

| Function | Signature | Purpose |
|---|---|---|
| `stand()` | `stand(data)` | Represent forest plot state; `data` is a dataframe with required `size||dbh` in mm + `species ID` + `plot_size` in m2 columns |
| `species_model()` | `species_model(x)` where `x` is a species vector or a `stand` object | Define species-level structure for the  IPM (kernel constructors) |
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

### `stand()` object

- Signature: `stand(data)` — `data` must contain `size` or `dbh` (numeric > 0), `species` or `species_id` (must match known IDs), and `plot_size` columns
- Internal structure: `stand$trees` (individual-level data), `stand$species` (unique species vector), `stand$plot_size`
- S3 class: `"ipm_stand"`
- `stand` is immutable: `project()` returns new stand states, never modifies in-place

### `species_model()` object

- Defines structural IPM(s): kernel constructors and parameters (L, Lmax, h bins)
- Can be initialized from an explicit species vector (`species_model(c("Picea_sp", "Abies_sp"))`) or inferred from a stand (`species_model(stand)`)
- Single- and multi-species share the same class (`"ipm_model"`)
- If species without parameters (sp NOT in `supported_species()`) is listed, three options are available: `error` to return an verbose error to the user; `drop` for dropping the species from the list and removing it to the `stand`; `static` to keep the species for background competition estimation but not `lambda` or `project` calculation. Default of `missing` argument to `error`
- Internal: `mod$species` is a character vector

### `parameters()` object

- Signature: `parameters(mod, draw = "random", seed = NULL)` — takes `mod` as first arg (not species ID)
- `draw`: `"mean"` uses posterior mean; `"random"` samples one draw; integer (1-1000) selects exact posterior draw
- Parameters always stored as multi-species structure internally, even for single species
- Internal structure: `pars$species_params$<species_id>$fixed`, `$random_effects`, `$draw_id`
- Random effects live inside `parameters()` — they represent a model realization, not ecological state
- Future extension: `parameters(mod, draw = 12, random_effect = spatial_re(lat, lon))`

### `env_conditions()` naming and behavior

- Named `env_conditions()` — avoids collision with base R `environment()`
- Accepts scalars or functions of time: `MAT = function(t) 4 + 0.02 * t`
- MAT/MAP scaling to [0,1] applied internally with interal parameters; never exposed to user

### Engine function signatures

- `lambda(mod, pars, stand, env)` — `mod` is first argument
- `project(mod, pars, stand, env, ctrl)` — `mod` is first argument
- `plot_effects` argument removed from both engines; random effects are resolved inside `parameters()`
- `lambda()` returns a named numeric vector (one λ per species in `mod`)
- `lambda()` does not require parameters for heterospecific species — competition is computed from stand structure only; no dynamic coupling

### Multi-species projections

- `species_model(c("ABIBAL", "ACERUB"))` + `parameters(mod, ...)` — unified interface
- `project()` with multiple species: fully coupled dynamics (competition updated each timestep)
- `project()` with single species: independent dynamic
- Validation: `pars` must contain parameters for all species in `mod`; species in `stand` must be included in `mod`

### Error handling

- Use `cli::cli_abort()`, `cli_warn()`, `cli_inform()` throughout — `cli` is a package dependency
- `species_model()` with unsupported ID → error with closest-match suggestion (similar to usethis "did you mean X?")
- Cloud unreachable → hard error: `"Could not reach cloud host. Check your internet connection or use a previously cached species_model."`
- All constructor inputs validated at construction time (fail fast, not deferred to projection time)

### Claude's Discretion

- Exact output structure of `ipm_projection` object (beyond what the architecture doc specifies)
- S3 class naming and method conventions (print, summary, plot methods)
- Internal validation helper structure
- Naming conventions for internal-only functions

</decisions>

<specifics>
## Specific Ideas

- Architecture document at `.planning/IPM_PACKAGE_ARCHITECTURE_v2.md` is the primary design reference — the API design document produced by Phase 1 should build on and formalize it, not contradict it
- The package is described as: *"A deterministic dynamic forest simulator built on Bayesian species-specific IPMs with endogenous competition and flexible climate drivers"*
- Canonical workflow (from architecture v2, section "Full Workflow"):
  ```r
  stand <- stand(data)
  mod   <- species_model(stand)
  pars  <- parameters(mod, draw = "random", seed = 42)
  env   <- env_condition(MAT = 8, MAP = 1200)
  ctrl  <- control(years = 100, delta_time = 1)
  lambda(mod, pars, stand, env)
  project(mod, pars, stand, env, ctrl)
  ```
- Stochasticity is fully resolved in `parameters()` — everything downstream is deterministic; no hidden randomness inside projection
- Reproducibility: all stochasticity controlled through `parameters(draw, seed)`

</specifics>

<deferred>
## Deferred Ideas

- lat/lon → climate lookup helper `get_climate(lat, lon)` — requires spatial climate database, deferred to future milestone
- Spatial random effects: `parameters(mod, draw = 12, random_effect = spatial_re(lat, lon))` — architecture supports it, not in v1
- Multi-species community IPM extensions beyond the current unified model — scoped as v3+

</deferred>

---

*Phase: 01-api-design*
*Context gathered: 2026-03-02*
