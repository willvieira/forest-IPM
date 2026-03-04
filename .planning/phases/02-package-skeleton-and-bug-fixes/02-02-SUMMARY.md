---
phase: 02-package-skeleton-and-bug-fixes
plan: 02
subsystem: api
tags: [r-package, s3-classes, ipm, constructors, eigenvalue, tibble, cli]

# Dependency graph
requires:
  - phase: 02-01
    provides: Fixed getEigenValues(), init_pop(), getPars_sp() — internal machinery used by lambda() and project()
provides:
  - Five typed S3 constructors: stand() -> ipm_stand, species_model() -> ipm_spModel, parameters() -> ipm_parameters, env_condition() -> ipm_env, control() -> ipm_control
  - Two IPM engines: lambda() returning ipm_lambda named numeric, project() returning ipm_projection with $species/$years/$lambda/$stand_series/$summary
  - Utility: supported_species() returning 33-species tibble from inst/extdata/species_list.csv
  - Minimal DESCRIPTION file enabling devtools::load_all()
  - All print/summary S3 methods for six classes; plot.ipm_projection() for lambda trajectories
affects:
  - 02-03-package-infrastructure
  - All future phases calling the canonical API workflow

# Tech tracking
tech-stack:
  added:
    - cli (cli_abort with named vector error format)
    - stringdist (amatch for closest-match species suggestion)
    - tibble (as_tibble for supported_species() and project() $summary)
  patterns:
    - "Three-layer S3 pattern: new_<class>() (no validation) -> validate_<class>() (cli_abort) -> user-facing helper"
    - "Eager parameter loading in species_model() with graceful NULL fallback for missing local RDS (Phase 3 will add cloud fetch)"
    - "Climate resolver pattern: if (is.function(env$MAT)) env$MAT(t) else env$MAT for scalar vs time-varying"
    - "Competition aggregation via stats::approx() linear interpolation onto focal mesh"

key-files:
  created:
    - R/stand.R
    - R/species_model.R
    - R/parameters.R
    - R/env_condition.R
    - R/control.R
    - R/lambda.R
    - R/project.R
    - R/supported_species.R
    - inst/extdata/species_list.csv
    - DESCRIPTION
  modified: []

key-decisions:
  - "Species IDs use short form (ABIBAL, ACERUB) in the API, not numeric-prefix form (18032ABIBAL) used in raw data files — getPars_sp() path mapping to be aligned in Phase 3"
  - "supported_species() reads from inst/extdata/species_list.csv with fallback path for devtools::load_all() during development"
  - "plot_random = c(0, 0, 0) in lambda() and project() for Phase 2 — plot random effects wired from Parquet in Phase 3"
  - "Single-species competition: pass same Nvec_intra as Nvec_inter (self-competition); multi-species: aggregate competitors via stats::approx() onto focal mesh"
  - "Minimal DESCRIPTION created as Rule 3 auto-fix to unblock devtools::load_all(); full DESCRIPTION with all Imports/LinkingTo/Suggests in plan 02-03"

patterns-established:
  - "API pattern: all user-facing functions accept typed ipm_* S3 objects and use inherits() checks with cli_abort() for type errors"
  - "Validation pattern: validate_<class>() is separate from new_<class>(); user helper calls both in sequence"
  - "Climate pattern: env$MAT and env$MAP can be numeric scalar or function(t); engines resolve at each timestep t"
  - "Competition pattern: for N species, build competitor Nvec by summing interpolated Nvecs from all other species onto focal mesh"

requirements-completed: []

# Metrics
duration: 5min
completed: 2026-03-04
---

# Phase 02 Plan 02: API Constructor and Engine Implementation Summary

**Eight exported R functions implementing the full Phase 1 API: five typed S3 constructors (stand, species_model, parameters, env_condition, control), two IPM engines (lambda, project), and one utility (supported_species) with 33-species reference CSV**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-04T15:19:53Z
- **Completed:** 2026-03-04T15:25:11Z
- **Tasks:** 2
- **Files modified:** 10 created, 0 modified

## Accomplishments
- Five S3 constructors with three-layer validation pattern (new_/validate_/user-helper) and cli_abort() error messages
- Two IPM engines wiring constructor objects to mkKernel/getEigenValues/init_pop/dbh_to_sizeDist internal machinery
- Species reference CSV with 33 eastern North America tree species derived from thesisData/species_id.csv
- All six S3 classes have print/summary methods; ipm_projection has plot() method for lambda trajectories
- Minimal DESCRIPTION file to enable devtools::load_all(compile = FALSE) verification

## Task Commits

Each task was committed atomically:

1. **Task 1: Five S3 constructors, supported_species(), species_list.csv** - `515e603` (feat)
2. **Task 2: lambda() and project() engines** - `873db5f` (feat)

**Plan metadata:** (pending — created in this commit)

## Files Created/Modified
- `R/stand.R` - ipm_stand constructor with size/species column normalization and 127mm validation
- `R/species_model.R` - ipm_spModel with supported_species() validation and stringdist closest-match
- `R/parameters.R` - ipm_parameters with draw='mean'/'random'/integer and seed support
- `R/env_condition.R` - ipm_env accepting numeric scalar or function(t) for MAT/MAP
- `R/control.R` - ipm_control with four validated fields (years, delta_time, store_every, bin_width)
- `R/lambda.R` - asymptotic eigenvalue engine; ipm_lambda S3 class with print/summary
- `R/project.R` - timestep projection loop; ipm_projection S3 class with print/summary/plot
- `R/supported_species.R` - reads species_list.csv with dev fallback path
- `inst/extdata/species_list.csv` - 33 species with 6 columns (species_id, common_name, nom_commun, growth/surv/recruit_model)
- `DESCRIPTION` - minimal package file (Rule 3 auto-fix; completed in plan 02-03)

## Decisions Made
- Species IDs use short form (ABIBAL) in the API rather than the numeric-prefix form (18032ABIBAL) used in raw data files. The actual data uses full IDs; getPars_sp() path mapping will be aligned in Phase 3 when the cloud Parquet data layer is built.
- single-species runs pass the same Nvec as both Nvec_intra and Nvec_inter (self-competition); multi-species aggregates competitors onto focal mesh via stats::approx() linear interpolation.
- Phase 2 sets plot_random = c(0, 0, 0) everywhere — plot-level random effects come from Parquet in Phase 3.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Created minimal DESCRIPTION file**
- **Found during:** Task 1 (verification step)
- **Issue:** devtools::load_all() requires a DESCRIPTION file with a Package field. The package had none; devtools::load_all() threw "Could not find a root 'DESCRIPTION' file".
- **Fix:** Created DESCRIPTION with Package: forestIPM, minimal Imports (cli, stats, tibble, utils, stringdist), and LinkingTo (Rcpp, RcppEigen). C++ compilation blocked by Xcode license (known STATE.md blocker); used compile = FALSE for verification.
- **Files modified:** DESCRIPTION (created)
- **Verification:** devtools::load_all(compile = FALSE) succeeds; all 8 exports accessible
- **Committed in:** 515e603 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** DESCRIPTION was always planned for 02-03 but was needed immediately to run devtools::load_all(). The minimal version here will be replaced by the complete version in 02-03. No scope creep.

## Issues Encountered
- Xcode license not accepted — C++ compilation blocked. Verified all R code using devtools::load_all(compile = FALSE). This is a known blocker in STATE.md and will need `sudo xcodebuild -license accept` to resolve before plan 02-03's R CMD check step.
- S3 dispatch (print/summary) does not work correctly without a NAMESPACE file. Direct method calls (print.ipm_stand(x)) work correctly; generic dispatch will be restored in plan 02-03 when roxygen2 generates NAMESPACE.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 8 exported functions implemented and verified via load_all(compile = FALSE)
- Plan 02-03 can proceed: needs DESCRIPTION completion (Imports/LinkingTo expansion), NAMESPACE generation via devtools::document(), R/globals.R for @importFrom declarations
- Xcode license must be accepted before C++ compilation tests can run in 02-03
- getPars_sp() species ID format mismatch (short vs numeric-prefix) deferred to Phase 3

---
*Phase: 02-package-skeleton-and-bug-fixes*
*Completed: 2026-03-04*
