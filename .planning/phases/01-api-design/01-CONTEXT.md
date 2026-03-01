# Phase 1: API Design - Context

**Gathered:** 2026-03-01
**Status:** Ready for planning

<domain>
## Phase Boundary

Define the public API surface of the forestIPM package — exported function signatures, argument names, return types, and user-facing conventions — before any implementation begins. Output is an API design document that guides Phase 2 implementation.

</domain>

<decisions>
## Implementation Decisions

### Overall API shape

- **Low-level + new API pattern**: export only the new high-level API; make existing low-level functions (`mkKernel`, `init_pop`, `getPars_sp`, `size_to_BAcomp`, `dbh_to_sizeDist`) internal
- `run_ipm()` from requirements is replaced by `project()` + `lambda()` — the architecture document (`.planning/IPM_PACKAGE_ARCHITECTURE.md`) defines the full design and should be treated as the primary specification
- No unnecessary wrapper functions; researchers compose constructors + engine functions directly

### Public API surface (final)

**Constructors:**

| Function | Signature | Purpose |
|---|---|---|
| `species_model()` | `species_model(species_id)` | Fetch full posterior from cloud (4000 draws) |
| `parameters()` | `parameters(sp, draw = "random", seed = NULL)` | Resolve a single parameter realization |
| `stand()` | `stand(size, species, plot_size)` | Represent forest plot state (immutable) |
| `env_conditions()` | `env_conditions(MAT, MAP)` | Represent exogenous climate drivers |
| `control()` | `control(years = 100, delta_time = 1)` | Simulation configuration |

**Core functions:**

| Function | Signature | Purpose |
|---|---|---|
| `project()` | `project(params, stand, env, plot_effects = rep(0, 3), control = NULL)` | Project population dynamics over time |
| `lambda()` | `lambda(params, stand, env, plot_effects = rep(0, 3))` | Single-step dominant eigenvalue |
| `supported_species()` | `supported_species()` | Returns tibble of available species IDs, names, model variants |

### Climate inputs

- Raw `MAT` (mean annual temperature) and `MAP` (mean annual precipitation) values only — no lat/lon lookup
- `env_conditions()` accepts scalars or functions of time: `MAT = function(t) 4 + 0.02 * t`
- All inputs converted to functions internally; scaling (0–1 normalization) applied inside `project()`, never exposed to user
- lat/lon → climate lookup deferred (requires spatial climate database, out of scope for v1)

### Multi-species projections

- `project(params, stand, env)` where `params` is a **named list** keyed by species ID when running multi-species:
  ```r
  project(
    params = list(ABIBAL = params_ab, ACERU = params_ac),
    stand = st,
    env = env
  )
  ```
- Single-species: pass a single `species_params` object (not a list)

### Plot-level random effects

- Separate argument on `project()` and `lambda()`: `plot_effects = rep(0, 3)` (default = zeros = no plot effect)
- Vector of three values: `[1]` growth offset, `[2]` mortality offset, `[3]` recruitment offset

### `env_conditions()` naming

- Named `env_conditions()` — avoids collision with base R `environment()` function
- Not `environment()`, not `conditions()`, not `forest_env()`

### Error handling

- Use `cli::cli_abort()`, `cli_warn()`, `cli_inform()` throughout — `cli` is a package dependency
- `species_model()` with unsupported ID → error with closest-match suggestion (similar to usethis "did you mean X?")
- Cloud unreachable → hard error: `"Could not reach cloud host. Check your internet connection or use a previously cached species_model."`
- All constructor inputs validated at construction time (fail fast, not deferred to projection time)

### Constructor validation rules

- `stand()`: `size` must be numeric and > 0; `species` must match known IDs; `plot_size` must be > 0
- `env_conditions()`: `MAT` and `MAP` must be numeric or single-argument functions
- `parameters()`: `draw` must be `"mean"`, `"random"`, or integer in 1–4000 range
- `control()`: `years` > 0; `delta_time` > 0

### Claude's Discretion

- Exact output structure of `ipm_projection` object (beyond what the architecture doc specifies)
- S3 class naming and method conventions (print, summary, plot methods)
- Internal validation helper structure
- Naming conventions for internal-only functions

</decisions>

<specifics>
## Specific Ideas

- Architecture document at `.planning/IPM_PACKAGE_ARCHITECTURE.md` is the primary design reference — the API design document produced by Phase 1 should build on and formalize it, not contradict it
- The package is described as: *"A deterministic dynamic forest simulator built on Bayesian species-specific IPMs with endogenous competition and flexible climate drivers"*
- Default workflow from architecture (section 8) should appear verbatim in the API design doc as the canonical example
- `parameters(draw = "mean")` uses posterior mean; `draw = "random"` samples one draw; `draw = 42` selects exact posterior draw 42
- `stand` is immutable: `project()` returns new stand states, never modifies in-place
- Stochasticity is fully resolved in `parameters()` — everything downstream is deterministic

</specifics>

<deferred>
## Deferred Ideas

- lat/lon → climate lookup helper `get_climate(lat, lon)` — requires spatial climate database, deferred to future milestone
- Multi-species community IPM (`community_project()` or similar) — scoped as v3+
- Named `run_ipm()` convenience wrapper — replaced by `project()` + `lambda()` in this design; if later needed, scope as a separate request

</deferred>

---

*Phase: 01-api-design*
*Context gathered: 2026-03-01*
