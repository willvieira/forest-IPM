# Forest Dynamic IPM R Package

## Architecture & API Design Specification

---

# 1. Purpose of the Package

This package implements a **dynamic, species-specific forest Integral Projection Model (IPM)** fitted for 31 species using Bayesian inference (4000 posterior draws per parameter).

The model supports:

1. Estimating population growth rate (λ) under specific stand and climate conditions.
2. Projecting single-species dynamics over time.
3. Projecting full forest community dynamics with endogenous competition.
4. Running climate change scenarios via time-varying environmental drivers.
5. Propagating parameter uncertainty via posterior draw selection.

The package is designed to be:

* Reproducible
* Modular
* Low-level but coherent
* Extensible to more complex ecological processes
* Suitable for research-grade simulation

---

# 2. Core Design Principles

## 2.1 Separation of Concerns

Inputs are separated into four conceptual domains:

| Domain             | Object                            | Responsibility                               |
| ------------------ | --------------------------------- | -------------------------------------------- |
| Model              | `species_model`, `species_params` | Encapsulate fitted Bayesian parameters       |
| Stand              | `stand`                           | Represent the biological state of the forest |
| Environment        | `environment`                     | Represent exogenous climate drivers          |
| Simulation Control | `control`                         | Represent simulation configuration           |

The projection engine (`project()`) orchestrates these components but does not contain biological or scaling logic beyond coordination.

---

## 2.2 Deterministic Core

All stochasticity (posterior sampling) occurs in:

```
parameters()
```

After parameter realization, the projection engine is fully deterministic.

---

## 2.3 Endogenous Competition

Competition is:

* A deterministic function of stand structure
* Recomputed at every time step
* Not manually specified as matrices by users

Time-varying competition emerges naturally from stand dynamics.

Optional modifiers may be applied for experimental scenarios.

---

## 2.4 Extensibility

The architecture must allow future inclusion of:

* Disturbance regimes
* Harvesting modules
* Recruitment pulses
* Fire modules
* Stochastic climate
* Community-level extensions

Therefore, input objects must remain structured and extensible.

---

# 3. Object Definitions

---

# 3.1 `species_model`

Represents the fitted Bayesian IPM model for a species.

## Constructor

```r
species_model(species_id)
```

## Contains

* Posterior parameter matrix (4000 draws)
* Parameter names
* Climate scaling metadata (MAT, MAP scaling ranges)
* Species metadata
* Model fitting version

## Does NOT:

* Collapse posterior uncertainty
* Select parameter draws

---

# 3.2 `species_params`

Represents a single parameter realization derived from `species_model`.

## Constructor

```r
parameters(sp, draw = "random", seed = NULL)
```

### Arguments

* `draw`

  * Integer index (1–4000)
  * `"random"`
  * `"mean"`

* `seed`

  * Optional seed used when `draw = "random"`

## Behavior

* If `draw = integer`, select exact posterior draw
* If `"random"`, sample one draw
* If `"mean"`, compute posterior mean

## Output

* Numeric parameter list ready for kernel construction
* Metadata:

  * draw_id
  * seed (if applicable)
  * selection strategy

After this stage, simulations are deterministic.

---

# 3.3 `stand`

Represents the biological state of a forest plot.

## Constructor

```r
stand(
  size,
  species,
  plot_size
)
```

### Arguments

* `size`: numeric vector of individual tree sizes
* `species`: species identifier vector (same length)
* `plot_size`: area of plot

## Responsibilities

* Store individual-level tree data
* Validate inputs
* Provide methods to compute:

  * Basal area
  * Conspecific competition
  * Heterospecific competition

Competition metrics are computed internally during projection.

## Important Design Rule

`stand` is treated as immutable.
Projection returns updated stand states rather than modifying in-place.

---

# 3.4 `environment`

Represents exogenous climate drivers.

## Constructor

```r
environment(
  MAT,
  MAP
)
```

## Accepted Inputs

Each driver may be:

* Numeric constant
* Function of time `function(t)`

Examples:

```r
environment(MAT = 4, MAP = 800)
```

```r
environment(
  MAT = function(t) 4 + 0.02 * t,
  MAP = 800
)
```

## Internal Behavior

* All inputs converted to functions internally
* Raw values preserved
* Scaling (0–1 normalization) applied during projection, not at construction

## Scaling

Scaling uses metadata stored in `species_model`.

Scaling is performed inside projection:

```
MAT_scaled <- scale_01(MAT_raw, scaling_info)
```

---

# 3.5 `control`

Represents simulation configuration.

## Constructor

```r
control(
  years = 100,
  delta_time = 1
)
```

## Behavior

* Defaults applied if not specified
* Passed into `project()`
* If `control = NULL`, `project()` creates default internally

`delta_time` is fixed at 1 by default and discouraged from modification unless biologically justified.

---

# 4. Projection Engine

---

# 4.1 `project()`

## Signature

```r
project(params, stand, env, control = NULL)
```

## Behavior

If `control` is `NULL`, create default via:

```r
control()
```

## Internal Algorithm

For each time step:

1. Evaluate environmental drivers:

   ```
   MAT_raw <- env$MAT(t)
   MAP_raw <- env$MAP(t)
   ```

2. Scale climate variables.

3. Compute competition from current stand:

   ```
   comp <- compute_competition(stand, focal_species)
   ```

4. Build IPM kernel:

   ```
   kernel <- build_kernel(params, MAT_scaled, MAP_scaled, comp)
   ```

5. Update stand state.

6. Store results.

---

# 4.2 Competition Handling

Competition is:

```
comp = f(stand_state)
```

Time variation emerges via stand dynamics.

Optional experimental modifier may be introduced:

```r
competition_modifier = function(t, comp) comp * 0.8
```

But competition matrices are never user-specified manually.

---

# 5. Lambda Calculation

## Function

```r
lambda(params, stand, env)
```

Computes dominant eigenvalue for given conditions.

Uses same internal pipeline but without iterative stand updating.

---

# 6. Community Projection

Community is represented as:

* A single `stand` object with multiple species.
* Species-specific kernels applied to appropriate subsets.

Competition is fully endogenous.

No manual heterospecific matrix input allowed.

---

# 7. Output Object

`project()` returns an `ipm_projection` object containing:

* Initial stand
* Parameter metadata
* Environment specification
* Control configuration
* Stand history (or summary statistics)
* Possibly λ estimates

Full provenance is stored for reproducibility.

---

# 8. Default Workflow Example

```r
sp      <- species_model("Picea")
params  <- parameters(sp, draw = 42)

st      <- stand(size = ..., species = ..., plot_size = 400)

env     <- environment(MAT = 4, MAP = 800)

out <- project(params, st, env)
```

Advanced climate scenario:

```r
env <- environment(
  MAT = function(t) 4 + 0.02 * t,
  MAP = 800
)

ctrl <- control(years = 200)

out <- project(params, st, env, control = ctrl)
```

---

# 9. Invariants & Constraints

* Posterior uncertainty resolved before projection.
* Climate scaling never exposed to user.
* Competition always derived from stand.
* `project()` never mutates input objects.
* Defaults applied via `control()` internally.

---

# 10. Future Extension Points

Potential expansions supported by this design:

* Disturbance modules
* Thinning events
* Fire dynamics
* Recruitment pulses
* Stochastic climate drivers
* Parallel posterior propagation
* Multi-plot simulation
* Integration with spatial modules

---

# 11. Summary of Final Public API

### Constructors

* `species_model()`
* `parameters()`
* `stand()`
* `environment()`
* `control()`

### Core Functions

* `project()`
* `lambda()`

No unnecessary wrappers (e.g., `lambda_posterior()`).

---

# 12. Architectural Identity

This package is:

> A deterministic dynamic forest simulator built on Bayesian species-specific IPMs with endogenous competition and flexible climate drivers.

It is not merely an IPM kernel wrapper — it is a structured ecological simulation engine.
