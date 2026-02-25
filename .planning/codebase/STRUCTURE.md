# Codebase Structure

**Analysis Date:** 2025-02-24

## Directory Layout

```
forest-IPM/
├── R/                           # Core package functions
│   ├── vital_rates.R            # Demographic vital rate models
│   ├── kernel.R                 # IPM kernel assembly and population initialization
│   ├── BasalArea_competition.R  # Competition metric computation
│   ├── params.R                 # Parameter loading and transformation
│   └── matrix_image.R           # Visualization utilities
├── src/                         # C++ source code
│   └── eigen.cpp                # RcppEigen wrapper for eigenvalue computation
├── simulations/                 # Simulation experiments directory
│   ├── coexistence/             # Two-species coexistence scenarios
│   ├── lambda_plot/             # Population growth rate analysis
│   ├── model_lambdaPlot/        # Model comparison lambda studies
│   ├── sensAnalysis_v3/         # Sensitivity analysis studies
│   ├── uncertainty_sim/         # Uncertainty quantification simulations
│   └── [other study dirs]/      # Additional analyses
├── _data.path                   # Plain text file with external data directory path
├── CLAUDE.md                    # Project guidelines for Claude Code
├── README.md                    # Project description
└── LICENSE                      # License file
```

## Directory Purposes

**R/**
- Purpose: Core IPM package implementation
- Contains: Vital rate functions, kernel assembly, competition metrics, parameter utilities, visualization
- Key files: `vital_rates.R`, `kernel.R`, `BasalArea_competition.R`

**src/**
- Purpose: C++ implementation for performance-critical computations
- Contains: Eigenvalue solvers via RcppEigen
- Key files: `eigen.cpp`

**simulations/**
- Purpose: Research experiments using the IPM package
- Contains: Subdirectories for each research question (coexistence, lambda dynamics, sensitivity, etc.)
- Pattern: Each subdirectory has entry points `run_ipm.R` (parameter setup) and `ipm_i.R` (simulation execution)
- Generated files: `simulation_pars.csv` (parameter grid), `*.RDS` (simulation outputs under `output/`)

## Key File Locations

**Entry Points:**

**Simulation Cluster Entry:**
- `simulations/coexistence/run_ipm.R`: Generates parameter grid, writes SLURM job submission scripts
- `simulations/coexistence/ipm_i.R`: Executes individual simulation replicate (called by SLURM array job)
- Pattern: `ipm_i.R` sources all `R/*.R` files and reads CSV parameters via environment variables

**Core Implementation:**
- `R/vital_rates.R`: Growth model (vonBertalanffy_f), survival model (survival_f), ingrowth model (ingrowth_f)
- `R/kernel.R`: Probability density wrappers (vonBertalanffy_lk, ingrowth_lk), kernel assembly (mkKernel), population initialization (init_pop)
- `R/BasalArea_competition.R`: Size conversions (size_to_BAind), plot metrics (size_to_BAplot), competition (size_to_BAcomp), size distribution conversion (dbh_to_sizeDist)
- `R/params.R`: Parameter loading (getPars_sp), format transformation (pars_to_list)
- `R/matrix_image.R`: Kernel visualization (matrix.image)

**C++ Performance:**
- `src/eigen.cpp`: Fast eigenvalue computation via getEigenValues() exported to R

## Naming Conventions

**Files:**
- Function files: `snake_case.R` (e.g., `vital_rates.R`, `BasalArea_competition.R`)
- Configuration/metadata: lowercase (e.g., `_data.path`)
- Simulation scripts: `[descriptor]_[version].R` or simple descriptive names (e.g., `run_ipm.R`, `ipm_i.R`, `fig1.R`)

**Functions:**
- camelCase for user-facing exports (e.g., `mkKernel()`, `init_pop()`, `size_to_BAplot()`)
- snake_case mixed with camelCase in vital rate models (e.g., `vonBertalanffy_f`, `survival_f`, `ingrowth_f`)
- snake_case for internal/helper functions (e.g., `size_to_BAind()`, `size_to_BAcomp()`)
- Short abbreviations in simulation code (e.g., `sp1`, `sp2`, `BA_comp_intra`, `Ntp1` for N at time plus 1)

**Variables in Simulation Context:**
- Species identifiers: `sp1`, `sp2` (strings)
- Population vectors: `N0_sp1`, `Nt_sp2`, `Ntp1_sp2` (list objects with meshpts, Nvec, h)
- Kernels: `K0_sp1`, `Kt_sp2`, `Kt_sp1` (matrix objects)
- Temporal: `t` for timestep counter, `delta_time` for time intervals in years
- Data frames: `year_summ` for annual summaries with columns year, N_sp*, BA_sp*, lambda*

**Directories:**
- kebab-case for simulation study names (e.g., `coexistence`, `lambda_plot`, `model_lambdaPlot`)
- Versioning with underscores and numbers (e.g., `sensAnalysis_v3`)

## Where to Add New Code

**New Demographic Model:**
- Primary code: `R/vital_rates.R` - Add new function with signature matching existing (e.g., `logistic_growth_f()` following `vonBertalanffy_f()` pattern)
- Usage: Update `R/kernel.R` to use new model in `P_xEC()` or ingrowth call
- Tests: Create validation script in `simulations/` directory if testing against empirical data

**New Analysis/Simulation Study:**
- Create new directory: `simulations/[study_name]/`
- Implement: `run_ipm.R` (parameter grid generation, SLURM script writing), `ipm_i.R` (simulation execution)
- Outputs: Save to `simulations/[study_name]/output/sim_[id].RDS` following coexistence pattern
- Visualizations: Create `fig1.R`, `fig2.R` etc. in study directory for results processing

**New Utility Function:**
- Shared helpers: `R/params.R` if parameter-related, or create new `R/[utilities].R` file if general-purpose
- Export pattern: Add `#' @export` roxygen comment above function if public, omit if internal
- Internal helpers: Use simple names without export comments

**New C++ Optimization:**
- Add functions to: `src/eigen.cpp`
- Compilation: Automatically compiled via Rcpp on package load/installation
- Export pattern: Add `// [[Rcpp::export]]` comment above C++ function
- Signature: Map R vectors/matrices via Eigen types (e.g., `Map<MatrixXd>`)

## Special Directories

**simulations/**
- Purpose: Contains all research experiments
- Generated: Yes - simulation output RDS files and parameter CSVs created at runtime
- Committed: Experiment scripts and figures committed; output/ and CSV files typically not committed (vary by .gitignore)

**.history/**
- Purpose: Automatic editor backup/version tracking (likely from Positron/RStudio)
- Generated: Yes - auto-created by IDE
- Committed: No - excluded via .gitignore

**.planning/**
- Purpose: GSD planning documents and codebase analysis
- Generated: Yes - created by GSD tools
- Committed: Varies by workflow (typically excluded or documented separately)

**.git/**
- Purpose: Git version control metadata
- Generated: Yes - git repository state
- Committed: No - excluded by .gitignore

## Loading and Sourcing Pattern

**In Simulation Scripts:**
```r
# Standard sourcing pattern in simulations/*/ipm_i.R
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')
```

**Parameter Access Pattern:**
```r
# Load from external path
pars_sp1 <- getPars_sp(sp = sp1, method = 'mean', path = file.path('..', 'data', 'parameters'))
```

**Typical Usage Sequence:**
1. Load parameters via `getPars_sp()` and optionally convert via `pars_to_list()`
2. Initialize populations via `init_pop()`
3. Build initial kernel via `mkKernel()` with population vectors
4. Project population: `Ntp1$Nvec <- Kt$K %*% Nt$Nvec`
5. Update kernel with new population for next iteration
6. Extract outcomes: population abundance via `sum(Nvec)`, basal area via `size_to_BAplot()`, growth rate via `eigen(K)$values`

---

*Structure analysis: 2025-02-24*
