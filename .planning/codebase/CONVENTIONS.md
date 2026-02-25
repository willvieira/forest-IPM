# Coding Conventions

**Analysis Date:** 2026-02-24

## Naming Patterns

**Files:**
- Descriptive camelCase or snake_case with underscores
- Examples: `vital_rates.R`, `BasalArea_competition.R`, `matrix_image.R`, `kernel.R`, `params.R`
- Pattern: descriptive of main functionality/domain

**Functions:**
- camelCase with underscores for multi-word functions
- Examples: `vonBertalanffy_f`, `vonBertalanffy_lk`, `P_xEC`, `ingrowth_lk`, `init_pop`, `size_to_BAind`, `size_to_BAplot`, `size_to_BAcomp`, `dbh_to_sizeDist`, `getPars_sp`, `pars_to_list`, `matrix.image`
- Internal helper functions use suffix `_f` for mathematical functions (e.g., `vonBertalanffy_f`) or `_lk` for likelihood functions (e.g., `vonBertalanffy_lk`)
- Exported functions may use dots (e.g., `matrix.image`) though underscore is more prevalent

**Variables:**
- camelCase preferred: `size_t0`, `size_t1`, `delta_time`, `BA_comp_intra`, `BA_comp_inter`, `pars`, `plot_random`, `Temp`, `Prec`
- Common abbreviations: `sp` (species), `BA` (basal area), `dbh` (diameter at breast height), `N` (population size), `K` (kernel matrix)
- Matrix/vector variables use capitalized abbreviations: `P`, `F`, `K`, `Nvec`, `N0_sp1`
- Mathematical parameters use mathematical notation as keys: `pars['r']`, `pars['Beta']`, `pars['Lmax']`, `pars['sigma_obs']`

**Types/Classes:**
- No explicit type declarations in R (dynamic typing)
- Named lists used for parameter objects: `pars[['growth']]`, `pars[['mort']]`, `pars[['rec']]`, `pars[['sizeIngrowth']]`

## Code Style

**Formatting:**
- Standard R formatting (no enforced linter in repository)
- Indentation: 1-2 spaces typical, inconsistent across files
- Line length: Generally under 80-100 characters
- Operators: Spaces around `=` and `<-` (mostly consistent)

**Linting:**
- No linting configuration detected (no `.Rprofile`, `lintignore`, or ESLint equivalent)
- No automated style enforcement

**Function Definition Pattern:**
```r
# Example from vital_rates.R
vonBertalanffy_f <- function(
  pars, delta_time, size_t0, BA_comp_intra, BA_comp_inter, Temp, Prec, plot_random
){
  # computation
  return( mu_obs )
}
```
- Opening brace on same line as function definition
- Parameters listed one per line or comma-separated
- `return()` statement explicit (not implicit)

## Import Organization

**Order:**
1. Base R functions (implicitly available)
2. External packages loaded via `library()` or `source()`

**Examples from simulation files:**
```r
library(tidyverse)        # Data manipulation
library(rstan)            # Bayesian inference
library(posterior)        # Posterior analysis
library(ggdist)           # Probability visualizations
```

**Path Aliases:**
- None detected; absolute/relative paths used in simulation scripts
- External data path stored in `_data.path` file (per CLAUDE.md)
- File paths hardcoded in scripts: `file.path('R/', 'vital_rates.R')`, `file.path('data', 'parameters')`

**Source Pattern:**
```r
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')
```

## Error Handling

**Patterns:**
- Use `stop()` with descriptive error messages for invalid inputs
- Example from `kernel.R`:
  ```r
  if(N < 0) {
    stop('Argument `N` must be larger or equal than zero.')
  }
  ```
- Example from `params.R`:
  ```r
  if(!MD %in% c('intcpt', 'intcpt_plot', ...)) {
    stop('model parameter should be one of the following character:\n...')
  }
  ```
- No try/catch (tryCatch) patterns found
- Validation occurs at function entry for parameter constraints
- Silent failures via `if(exists('fct'))` check in `kernel.R` (line 195) - checks for variable existence before use

## Logging

**Framework:** `console` (base R print/cat functions)

**Patterns:**
- Progress reporting via `cat()` with `\r` for in-place updates (test_envStoch.R):
  ```r
  cat('Progress ', round(count/2100 * 100, 1), '%\r')
  ```
- No structured logging framework (no `log4r`, `lgr`, etc.)
- Informal console output during interactive development/simulation runs

## Comments

**When to Comment:**
- Header comments with metadata: file purpose, author, date
  ```
  ##############################################################
  # Vital rate functions for growth, mortality and reproduction
  # Will Vieira
  # January 20, 2019
  # Last updated: March, 18, 2023
  ##############################################################
  ```
- Section headers for major function groups:
  ```r
  # 1. Growth
  # 2. Mortality
  # 3. Ingrowth
  ```
- Inline comments explain mathematical concepts and parameters
  Example from vital_rates.R:
  ```r
  # Compute r
  rPlotInd = exp(...)
  # pre calculate component of the model
  rPlotTime = exp(-rPlotInd * delta_time)
  ```

**JSDoc/TSDoc:**
- roxygen2 format used for exported functions (marked with `#' @export`)
- Examples from codebase:
  ```r
  #' @export
  mkKernel = function(...)

  #' @export
  init_pop <- function(

  #' Function to compute BA competition in function of population size vectors
  #' N_intra: size distribution for focal species; output of `init_pop` function
  #' N_inter: size distribution of competition species
  size_to_BAcomp <- function(...)
  ```
- Parameter documentation inline with roxygen comments rather than formal docstring blocks
- No return value documentation in roxygen comments

## Function Design

**Size:**
- Small functions (5-50 lines) preferred
- Examples: `size_to_BAind` (7 lines), `vonBertalanffy_f` (22 lines)
- Complex orchestration functions up to 70+ lines: `init_pop` (69 lines), `mkKernel` (45 lines)

**Parameters:**
- Functions accept explicit parameters, not `...` (variadic args)
- Parameters are densely packed: `mkKernel` accepts 8 parameters
- Complex objects (like `pars` list, `Nvec`) passed as single arguments containing multiple values
- Default parameters used sparingly: `accuracy = 0.001`, `meanSize = 130`, `sdSize = 1.8` in `init_pop`

**Return Values:**
- Explicit `return()` statements (not implicit final expression)
- Single return: scalar or simple object
- Multiple values returned as named list:
  ```r
  return(list(K = K, P = P, F = F))  # from mkKernel
  return(list(meshpts = msh, Nvec = N_out, h = h))  # from init_pop
  ```

## Module Design

**Exports:**
- roxygen2 `#' @export` marks public functions
- Exported functions in `kernel.R`: `mkKernel`, `init_pop`
- Exported functions in `BasalArea_competition.R`: `size_to_BAplot`, `size_to_BAcomp`, `dbh_to_sizeDist`
- Exported functions in `matrix_image.R`: `matrix.image`
- Parameter functions in `params.R`: `getPars_sp`, `pars_to_list` (marked exportable)
- Non-exported helper functions: `vonBertalanffy_f`, `survival_f`, `ingrowth_f`, `vonBertalanffy_lk`, `P_xEC`, `ingrowth_lk`, `size_to_BAind`

**Barrel Files:**
- No barrel/index files detected
- Each source file in `R/` is independent module with specific domain
- Simulation scripts source files individually as needed

## Data Flow Patterns

**Parameter Structure:**
- Demographic parameters grouped by vital rate type:
  ```r
  pars[['growth']]  # Growth parameters
  pars[['mort']]    # Mortality parameters
  pars[['rec']]     # Recruitment parameters
  pars[['sizeIngrowth']]  # Size at ingrowth distribution
  ```

**Population Vector Structure:**
- `Nvec` object contains:
  ```r
  list(
    meshpts = mesh_points_vector,
    Nvec = population_distribution_vector,
    h = bin_width
  )
  ```

**Climate Variables:**
- Standardized parameter names: `Temp` (temperature), `Prec` (precipitation)
- Typically passed as scaled values (0-1 range) derived from climate data

---

*Convention analysis: 2026-02-24*
