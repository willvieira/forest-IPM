# Technology Stack

**Analysis Date:** 2026-02-24

## Languages

**Primary:**
- R (version unspecified) - Core package with statistical functions, demographic models, and IPM kernels
- C++ - Low-level performance-critical numerical computation for eigenvalue solving

**Secondary:**
- Bash - SLURM cluster job submission scripts

## Runtime

**Environment:**
- R runtime (base R)

**Package Manager:**
- CRAN (R package manager)
- Lockfile: Not detected

## Frameworks

**Core:**
- Rcpp 1.0+ - C++ integration for R via RcppEigen
- RcppEigen - Eigen linear algebra library bindings for C++

**Data Processing:**
- tidyverse 1.3+ (dplyr, tidyr, ggplot2, purrr, stringr, readr) - Data manipulation, transformation, and visualization

**Statistical/Mathematical:**
- Base R statistics (dnorm, cutoff, density, lognormal distributions)
- Eigen (via RcppEigen) - Eigenvalue decomposition and linear algebra

**Testing:**
- Not detected

**Build/Dev:**
- Not detected

## Key Dependencies

**Critical:**
- truncnorm 1.0+ - Truncated normal distribution sampling; used in `R/kernel.R` for ingrowth size distribution modeling
- purrr - Functional programming tools (map, map_dbl) used in `R/params.R` for parameter processing
- tidyr - Data reshaping functions (pivot_longer, pivot_wider) used throughout for parameter parsing
- dplyr - Data frame manipulation (filter, select, mutate) used in simulation scripts
- stringr - String manipulation for parameter name parsing in `R/params.R`

**Infrastructure:**
- RcppEigen - C++ linear algebra acceleration for eigenvalue solver in `src/eigen.cpp`

## Configuration

**Environment:**
- External data path configured via `_data.path` (plain text file pointing to local data directory)
- SLURM cluster variables: BATCH, SLURM_ARRAY_TASK_ID (environment variables used in simulation scripts)
- Timeout configuration for downloads (options(timeout = 10000))

**Build:**
- Rcpp attributes used: `// [[Rcpp::export]]` and `// [[Rcpp::depends(RcppEigen)]]` in `src/eigen.cpp`
- No explicit build configuration files (CMakeLists.txt, Makefile) detected
- Package loads C++ code on runtime compilation via Rcpp

## Platform Requirements

**Development:**
- R runtime environment
- Rcpp toolchain (requires compiler: clang/gcc on macOS/Linux, MSVC on Windows)
- C++11 or later support for Eigen

**Production:**
- Deployment target: SLURM HPC cluster (uses `SBATCH` directives, array job submission)
- Requires R runtime with Rcpp/RcppEigen compilation support
- External data stored remotely (download from doc.ielab.usherbrooke.ca)

## External Data

**Remote:**
- treeData.RDS - Tree inventory dataset downloaded from HTTPS during simulation setup
- Location: https://doc.ielab.usherbrooke.ca/s/qQjXSatrVpGvKzg/download
- Cached locally after first download

**Local:**
- Plot parameters: RDS files stored in `plot_pars/` directory (species-specific: `pars_*.RDS`, `plotls_*.RDS`)
- Parameters for growth, mortality, recruitment models: stored in `data/output_sim_processed/` subdirectories
- Simulation initialization: `N_init.RDS`, `simulation_pars.csv`

---

*Stack analysis: 2026-02-24*
