# forest-IPM

## What This Is

An R package for running Bayesian hierarchical Integral Projection Models (IPMs) of tree population dynamics in eastern North America. It models growth (von Bertalanffy), survival, and recruitment as functions of tree size, intra/interspecific competition (basal area), and climate — covering 31 species across 3 demographic model variants. The package lets researchers compute population growth rates (lambda), sensitivities, and multi-decade projections without managing 15 GB of Bayesian posteriors locally.

## Core Value

Any researcher can run a full IPM for any supported species at any location by installing the package — no local data setup, no manual sourcing of R files, no cluster required.

## Requirements

### Validated

<!-- Existing capabilities confirmed in codebase -->

- ✓ Von Bertalanffy growth model (size × competition × climate) — existing
- ✓ Logistic survival model (size × competition × climate) — existing
- ✓ Ingrowth/recruitment model (plot-level, climate-driven) — existing
- ✓ IPM kernel assembly: P matrix (growth × survival) + F matrix (recruitment) via numerical integration — existing
- ✓ Competition metrics: intra/interspecific basal area from size distributions — existing
- ✓ Fast eigenvalue computation (lambda) via RcppEigen C++ — existing
- ✓ Parameter loading from local RDS files (posterior mean or random draw) — existing
- ✓ Population initialization with log-normal size distributions — existing
- ✓ Multi-decade population projection simulation — existing
- ✓ Kernel matrix visualization — existing

### Active

<!-- Building toward these — hypotheses until shipped -->

- [ ] Proper R package structure (DESCRIPTION, NAMESPACE, LICENSE, roxygen2 docs, installable via `devtools::install_github()`)
- [ ] Cloud parameter storage in Parquet format (31 species × 3 models hosted remotely)
- [ ] On-demand parameter fetching: partial reads from cloud without downloading full 15 GB
- [ ] Tiered parameter access: known plot locations served from cloud Parquet (full posterior fidelity); unobserved locations served via ML decoder (random effects generation from lat/lon)
- [ ] High-level workflow API: `run_ipm(species, coords, climate)` wrapping current internal functions
- [ ] S3 class definitions for core objects: `ipm_kernel`, `ipm_population`, `ipm_params`
- [ ] roxygen2 documentation for all exported functions with examples
- [ ] Unit tests (testthat) covering core functions and parameter fetching
- [ ] ML decoder for spatial random effects at unobserved locations

### Out of Scope

- Fitting new species to inventory data — this package consumes posteriors, not fits them
- Real-time climate data fetching — user provides climate inputs
- GUI or Shiny app — programmatic R API only
- Windows binary support initially — C++ compilation via Rcpp targets Linux/macOS first

## Context

The codebase currently lives in `R/` as a collection of standalone functions — not a proper R package (no DESCRIPTION, no NAMESPACE, no documentation infrastructure). Simulation scripts in `simulations/` source these files directly and run on SLURM clusters. The C++ eigenvalue solver in `src/eigen.cpp` is already Rcpp-compatible and will slot into the package naturally.

The parameter storage problem is the core packaging challenge: 31 species × 3 models × ~4k posterior draws × spatial random effects (thousands of plots per species) = ~15 GB. The spatial random effects dominate storage. The tiered approach — Parquet for known plots, ML decoder for new locations — solves both storage and the scientifically important extrapolation problem (predicting dynamics at unobserved sites).

Cloud hosting platform for Parquet files is TBD. Requirements: supports HTTP range requests for partial reads, accessible without authentication, reliable for research use. Candidates: AWS S3 (public bucket), Zenodo (with careful chunking), HuggingFace Hub.

## Constraints

- **Tech stack**: R + Rcpp/RcppEigen — must remain in this ecosystem, no Python or Julia core
- **Parameter size**: ~15 GB total — cloud hosting must support partial/lazy reads (Arrow + Parquet)
- **Scientific fidelity**: Full posterior distribution (4k draws) must be preservable for known plot locations
- **Reproducibility**: Package must be installable from GitHub with no local data setup

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Parquet for cloud parameter storage | Supports partial reads via Arrow; columnar format efficient for selecting species/draw subsets | — Pending |
| Tiered access: Parquet + ML decoder | Full fidelity for known plots; extrapolation for new locations; solves both storage and scientific gap | — Pending |
| Keep RcppEigen for eigenvalue solver | Already implemented and working; C++ speedup is critical for iterative projection | — Pending |
| Cloud hosting platform TBD | Needs HTTP range request support; decision depends on storage cost and access patterns | — Pending |

---
*Last updated: 2026-02-25 after initialization*
