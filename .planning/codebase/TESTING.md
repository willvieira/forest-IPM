# Testing Patterns

**Analysis Date:** 2026-02-24

## Test Framework

**Runner:**
- No formal test framework detected (no `testthat`, `RUnit`, or `tinytest`)
- Testing performed via exploratory R scripts in simulation directories

**Assertion Library:**
- Manual assertions using base R conditionals
- `stop()` for explicit failures in core functions
- No formal assertion/expectation library

**Run Commands:**
```bash
Rscript simulations/old/test_envStoch.R     # Run exploratory test script
Rscript simulations/model_lambdaPlot/fake_test.R  # Simulation validation
```

## Test File Organization

**Location:**
- Tests co-located with source in `/simulations/` directories (not separate `tests/` folder)
- Source code in `/R/` and `/src/`
- Test/exploratory scripts in `/simulations/old/`, `/simulations/model_lambdaPlot/`, etc.
- Pattern: separate directory structure rather than `src/*.test.R` co-location

**Naming:**
- Exploratory scripts: `test_envStoch.R`, `fake_test.R`
- Simulation runners: `run_ipm.R`, `run_randomForest.R`
- Not consistent with standard naming patterns (e.g., `*.test.R` or `test_*.R`)

**Structure:**
```
simulations/
├── old/
│   └── test_envStoch.R        # Environmental stochasticity validation
├── model_lambdaPlot/
│   └── fake_test.R            # Fake data generation and model validation
├── lambda_plotv0/
│   ├── run_ipm.R
│   ├── run_randomForest.R
│   └── ipm_i.R
└── species_pair/
    ├── run_ipm.R
    └── ipm_i.R
```

## Test Structure

**Suite Organization:**

From `simulations/old/test_envStoch.R`:
```r
# Setup section
library(tidyverse)
library(ggdist)

source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')

# Load external data
clim_rg <- readRDS('../TreesDemography/data/climate_range.RDS')
treeData = readRDS('../TreesDemography/data/treeData.RDS')

# Test scenario 1: Generate average population growth rate
# [test code block]

# Test scenario 2: Temperature stochasticity
count = 1
out_tb = tibble()
for(Temp in seq(0, 1, .05)) {
  for(clim_sd in seq(0, 2.25, 0.25)) {
    for(i in 1:10) {
      # test logic
    }
  }
}

# Test scenario 3: Competition stochasticity
# [similar nested loop structure]

# Test scenario 4: Combined temperature and competition stochasticity
# [similar nested loop structure]
```

**Patterns:**
- **Setup:** Load libraries, source functions, load external data
- **Scenario structure:** Nested loops for parameter sweeps
- **Progress tracking:** `cat()` with percentage updates
- **Output collection:** Accumulate results in tibbles via `rbind()`
- **Visualization:** Generate plots with `ggplot()` for manual inspection

From `simulations/model_lambdaPlot/fake_test.R`:
```r
# Fake data generation pattern
tibble(
  temp = runif(n, -3, 20)
) |>
bind_cols(
  beta = 0.1,
  Int = 0.2
) |>
mutate(
  # transformations
)

# Model fitting
out <- stan(
  data = data_stan,
  file = "model.stan",
  chains = 2,
  cores = 2
)

# Posterior analysis and visualization
out |>
  as_draws_df() |>
  # transformation pipeline
  ggplot() +
  # visualization code
```

**Key characteristics:**
- Data generation → Model fitting → Posterior analysis → Visualization pipeline
- Tests are exploratory/demonstrative, not automated assertions
- Manual visual inspection of plots is primary validation method

## Mocking

**Framework:**
- No mocking framework detected (no `mockery`, `unittest.mock` equivalent)

**Patterns:**
- Fake data generated inline using base R and tidyverse functions:
  ```r
  tibble(
    temp = runif(n, -3, 20)
  ) |>
  bind_cols(
    beta = 0.1,
    Int = 0.2,
    sigma_int = -2
  )
  ```
- External data loaded from RDS files (no stubbing)
- Mock parameters created by sampling from posterior:
  ```r
  pars_sp1 <- getPars_sp(
    sp = sp1,
    method = 'random',  # Sample from posterior instead of mean
    path = file.path('data', 'parameters')
  )
  ```

**What to Mock:**
- External parameter files (use `getPars_sp(method = 'random')` for sampling)
- Climate ranges (load from RDS or generate synthetic ranges)
- Tree inventory data (generate synthetic population distributions via `init_pop()`)

**What NOT to Mock:**
- Core IPM kernel functions (`mkKernel`, `vital_rates`) - these are tested directly
- Eigenvalue computation - tested via actual kernel matrix evaluation
- Mathematical transformations in `BasalArea_competition.R` - tested on synthetic data

## Fixtures and Factories

**Test Data:**

Population initialization factory in `kernel.R`:
```r
init_pop <- function(
  params,
  L,
  h,
  N,
  accuracy = 0.001,
  meanSize = 130,
  sdSize = 1.8
){
  # Generates smooth initial size distribution from lognormal
  # Returns: list(meshpts, Nvec, h)
}
```

Usage in tests:
```r
N0_sp1 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = 1
)
N0_sp2 <- init_pop(
  params = pars_sp1,
  L = 127,
  h = 1,
  N = 0
)
```

Parameter factory in `params.R`:
```r
getPars_sp <- function(
  sp,
  method,  # 'mean' or 'random'
  model = 'intcpt_plot_comp_clim',
  path = file.path('data', 'output_sim_processed')
)
# Returns posterior parameter estimates as named list
```

Fake data generation (from fake_test.R):
```r
tibble(
  temp = runif(n, -3, 20)
) |>
bind_cols(
  beta = 0.1,
  Int = 0.2,
  sigma_int = -2
) |>
mutate(
  # compute derived values
) ->
dat
```

**Location:**
- Factories are exported functions from core modules (`R/kernel.R`, `R/params.R`)
- No separate fixture directory
- Test data files in external data repository (`../TreesDemography/data/`)

## Coverage

**Requirements:**
- No coverage targets enforced
- No `.coverage` or coverage configuration files detected

**View Coverage:**
- Not applicable; formal coverage metrics not tracked

## Test Types

**Unit Tests:**
- Not formally structured, but core functions tested indirectly
- Example: `size_to_BAind()`, `size_to_BAcomp()`, `dbh_to_sizeDist()` called in test scripts
- Scope: Individual vital rate functions (growth, survival, ingrowth)
- Approach: Direct function calls with synthetic parameters, result inspection

**Integration Tests:**
- Primary test type in codebase
- Scope: Full IPM kernel assembly and eigenvalue computation
- Approach: Generate population vectors and parameters, call `mkKernel()`, extract eigenvalues
- Example from test_envStoch.R:
  ```r
  K0_sp1 = mkKernel(
    Nvec_intra = N0_sp1,
    Nvec_inter = N0_sp2,
    delta_time = pars$deltaTime,
    plotSize = pars$plotSize,
    Temp = temp_sample_scl,
    Prec = prec,
    pars = pars_sp1,
    plot_random = rep(0, 3)
  )
  lambda = max(Re(eigen(K0_sp1$K)$values))
  ```

**E2E Tests:**
- Not formally structured
- Bayesian inference validation in fake_test.R (Stan model fitting with fake data)
- Parameter uncertainty propagation tests (sampling from posterior in loops)

## Common Patterns

**Async Testing:**
- Not applicable (R is single-threaded, no async patterns)

**Error Testing:**

From `kernel.R`:
```r
if(N < 0) {
  stop('Argument `N` must be larger or equal than zero.')
}
else if(N == 0) {
  N_out <- rep(0, length(msh))
}
else {
  # normal computation
}
```

Error handling tested implicitly:
- Negative N values should trigger `stop()`
- Zero N should return empty distribution
- Positive N should return proper distribution

Manual error testing pattern in test scripts:
```r
# Set N = 0 and verify empty distribution
N0_sp2 <- init_pop(params = pars_sp1, L = 127, h = 1, N = 0)

# Set stochastic value that could be negative, clamp to 0
if(BA_stoch < 0)
  BA_stoch = 0
```

**Stochastic/Randomness Testing:**

Pattern in test_envStoch.R:
```r
for(i in 1:10) {  # 10 replicates
  # generate climate value from stochastic dist
  temp_sample = rnorm(1, mean = temp_unscl, sd = clim_sd)

  # run model with random parameter draw
  pars_sp1 <- getPars_sp(sp = sp1, method = 'random', path = ...)

  # compute lambda
  lambda = max(Re(eigen(K0_sp1$K)$values))

  # collect result
  out_tb = rbind(out_tb, tibble(temp = Temp, lambda = lambda, ...))
}
```
- No `set.seed()` for reproducibility (stochasticity intentional)
- Multiple replicates (10) to assess variability
- Posterior sampling via `getPars_sp(method = 'random')` for parameter uncertainty
- Results aggregated for visual inspection

**Visualization as Validation:**

All test scripts include plotting for manual inspection:
```r
out_tb |>
  ggplot() +
  aes(temp, lambda) +
  stat_pointinterval() +
  geom_smooth()
```
- Point-interval plots show distribution of outputs
- Smoothing reveals trends
- Facet plots compare scenarios
- Visual patterns indicate correct/incorrect behavior

---

*Testing analysis: 2026-02-24*
