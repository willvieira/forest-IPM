# Architecture

**Analysis Date:** 2025-02-24

## Pattern Overview

**Overall:** Modular Bayesian Hierarchical IPM (Integral Projection Model)

**Key Characteristics:**
- Demographic kernel construction via numerical integration over continuous size distributions
- Separation of vital rate functions (growth, survival, recruitment) from kernel assembly logic
- Competition metrics computed from population size vectors and converted to basal area
- Hybrid R/C++ architecture with R orchestrating models and C++ accelerating eigenvalue computation
- Stochastic simulation driven by parameter draws from posterior distributions

## Layers

**Vital Rates Layer:**
- Purpose: Define demographic functions that predict individual transitions (growth, survival, reproduction)
- Location: `R/vital_rates.R`
- Contains: Von Bertalanffy growth model, logistic survival model, ingrowth/recruitment functions
- Depends on: Species parameters (Bayesian posterior), environmental covariates (temperature, precipitation), competition metrics
- Used by: Kernel assembly layer to build transition probability densities

**Kernel Construction Layer:**
- Purpose: Assemble full IPM kernels (P and F matrices) using numerical integration
- Location: `R/kernel.R`
- Contains: Probability density functions for growth (von Bertalanffy likelihood), survival, ingrowth; kernel assembly via outer product integration; population initialization
- Depends on: Vital rates functions, competition metrics, species parameters
- Used by: Simulation layer for population projection

**Competition Metrics Layer:**
- Purpose: Convert tree size distributions to basal area competition values
- Location: `R/BasalArea_competition.R`
- Contains: Size-to-basal-area conversion, per-tree basal area computation, intra/interspecific competition aggregation
- Depends on: Population size vectors (mesh points and counts), plot size
- Used by: Kernel assembly and vital rates evaluation

**Visualization Layer:**
- Purpose: Display kernel matrices as heatmaps
- Location: `R/matrix_image.R`
- Contains: Generic matrix visualization with optional contours and legends
- Depends on: Kernel matrix (numeric matrix)
- Used by: Diagnostic and publication workflows

**Parameter Loading Layer:**
- Purpose: Access Bayesian posterior parameters for species
- Location: `R/params.R`
- Contains: Functions to load posterior distributions (mean or random draw), convert parameter formats for kernel use
- Depends on: External parameter files (RDS format), posterior summary statistics
- Used by: Simulation initialization

**Eigenvalue Computation Layer:**
- Purpose: Fast eigenvalue/eigenvector computation for population growth rate (lambda)
- Location: `src/eigen.cpp`
- Contains: RcppEigen wrapper for self-adjoint eigenvalue solver
- Depends on: RcppEigen C++ library
- Used by: Simulation tracking of population growth rates

## Data Flow

**Simulation Execution (primary flow):**

1. Parameter Initialization
   - Load species-specific parameters from posterior distribution via `getPars_sp()`
   - Convert to nested list format via `pars_to_list()`
   - Species parameters: growth (Lmax, r, Beta, tau_temp, tau_prec, optimal_temp, optimal_prec, sigma_obs), mortality (psi, etc.), recruitment (mPop_log, sigma_BA, optimal_BA, p_log, beta_p), ingrowth size (size_int, phi_time, sigma_size)

2. Population Initialization
   - Create initial size distribution via `init_pop()` using log-normal distribution
   - Output: Size vector with meshpoints (size classes) and Nvec (density per class)
   - Meshpoints define continuous size domain discretized into integration bins

3. Kernel Assembly Loop (per timestep)
   - Compute competition metrics: `size_to_BAcomp()` for intra/interspecific basal area above focal size class
   - Build P kernel (growth × survival): outer product of transition probabilities via `P_xEC()`
   - Build F kernel (recruitment): outer product of recruitment probabilities via `ingrowth_lk()`
   - Full kernel K = P + F via `mkKernel()`

4. Population Projection
   - Multiply kernel by population vector: N(t+1) = K * N(t)
   - Track population metrics: total abundance, basal area
   - Compute growth rate: lambda = max(Re(eigen(K)))
   - Stop when population reaches stable distribution (consecutive timestep changes < 1e-4)

**State Representation:**
- Population state: Named list with `meshpts` (size classes), `Nvec` (density per class), `h` (bin width)
- Kernel state: List with `K` (full kernel), `P` (growth × survival), `F` (recruitment)
- Temporal dynamics: Data frame with year, population counts, basal area, lambda

## Key Abstractions

**Population Vector (N):**
- Purpose: Represents continuous size distribution discretized into size classes
- Examples: `N0_sp1`, `Ntp1_sp2` in `simulations/coexistence/ipm_i.R`
- Pattern: List with numeric meshpoints vector, numeric Nvec density vector, and scalar bin width h

**IPM Kernel (K):**
- Purpose: Markov transition operator mapping current size distribution to next generation
- Examples: `K0_sp1`, `Kt_sp2` in `simulations/coexistence/ipm_i.R`
- Pattern: Numeric matrix (dim: size classes × size classes) assembled via numerical integration

**Basal Area Competition (BA_comp):**
- Purpose: Aggregated competitive effect of larger individuals on focal size class
- Examples: `BA_comp_intra`, `BA_comp_inter` in `R/kernel.R`
- Pattern: Numeric vector (length: size classes) computed from cumulative basal area above each meshpoint

**Species Parameters (pars):**
- Purpose: Container for all demographic and environmental parameters for a species
- Examples: `pars_sp1[['growth']]`, `pars[['rec']]` in `simulations/coexistence/ipm_i.R`
- Pattern: Nested named list with keys: 'growth', 'mort', 'rec', 'sizeIngrowth' mapping to named numeric vectors

## Entry Points

**Simulation Entry Point:**
- Location: `simulations/*/ipm_i.R`
- Triggers: SLURM array job (cluster) with environment variables BATCH and SLURM_ARRAY_TASK_ID
- Responsibilities: Read simulation parameters from CSV, initialize populations and kernels, execute multi-decade population projections with multiple sub-scenarios (single-species, invaded, coexistence), compute demographic outcomes

**Parameter Setup Entry Point:**
- Location: `simulations/*/run_ipm.R`
- Triggers: Manual execution on cluster head node
- Responsibilities: Generate parameter expansion grid (species pairs, replicates), call `getPars_sp()` for all species, write SLURM submission scripts with array jobs

**Interactive/Analysis Entry Point:**
- Location: User scripts (not in version control or under `simulations/*/` directories)
- Triggers: Manual R console execution
- Responsibilities: Source `R/*.R` files directly, call functions like `init_pop()`, `mkKernel()`, `size_to_BAplot()` for custom analyses

## Error Handling

**Strategy:** Input validation with informative error messages; no try-catch recovery patterns observed

**Patterns:**
- `getPars_sp()`: Validates `model` parameter is length 1 or 3, checks model names against whitelist
- `init_pop()`: Validates `N >= 0`, throws error for invalid population size
- Silent failures: If parameter files missing, `readRDS()` fails with file-not-found error

## Cross-Cutting Concerns

**Logging:** Inline progress reporting to console via `cat()` showing percentage complete in simulation loop (e.g., `cat(' Calculating', round(t/pars$n_time * 100, 2), '%\r')`)

**Validation:** Size class boundaries validated via `cut()` in `dbh_to_sizeDist()`; parameter bounds not explicitly enforced (assumed valid from Bayesian posterior)

**Authentication:** Not applicable (single-user research codebase)

**Environment Configuration:** External data path stored in `_data.path` file; species parameter files loaded from `data/output_sim_processed/` relative to working directory or via `path` argument

**Random Seed Management:** Set via `set.seed(pars$seed)` in simulation initialization for reproducibility; seed value stored in simulation parameters CSV

---

*Architecture analysis: 2025-02-24*
