# Codebase Concerns

**Analysis Date:** 2026-02-24

## Tech Debt

**Deprecated Parameter in Vital Rates:**
- Issue: Size-dependent mortality (dbh parameter in `survival_f`) is marked deprecated but code remains
- Files: `R/vital_rates.R` (lines 54-55)
- Impact: Dead code path that should either be removed or its deprecation properly handled; creates confusion for maintenance
- Fix approach: Either remove the commented-out line 68 or implement proper deprecation warning if functionality is conditionally used elsewhere

**Magic Numbers in IPM Integration:**
- Issue: Hardcoded numeric constants scattered throughout calculations without explanation
- Files:
  - `R/kernel.R` (line 60: `a = 127`) - minimum size threshold for ingrowth
  - `R/BasalArea_competition.R` (line 12: `1e-3`, line 25/49: `1e4`) - unit conversions
- Impact: Values lack clear documentation; unit conversions (mm to m², ha scaling) are opaque; makes code fragile to modifications
- Fix approach: Extract as named constants with units documentation; create documentation of unit assumptions (DBH in mm, BA in m²/ha)

**Unsafe Global State Dependency:**
- Issue: `init_pop` function uses `exists()` to check variable state (line 195 in `R/kernel.R`)
- Files: `R/kernel.R` (line 195)
- Impact: Function behavior depends on whether `fct` variable exists in parent environment; non-deterministic behavior; difficult to test
- Fix approach: Refactor to always compute/assign `fct` or make logic explicit without relying on variable existence

**Non-deterministic Population Initialization:**
- Issue: `init_pop` function uses random sampling (`runif`, `sample`) in a binary search loop without seed control
- Files: `R/kernel.R` (lines 165, 184)
- Impact: Same input parameters produce different results across calls; difficult to reproduce simulations; loop may not converge to specified accuracy
- Fix approach: Accept `seed` parameter; implement deterministic convergence (use binary search bounds more carefully)

## Known Bugs

**Infinite Loop Risk in init_pop:**
- Symptoms: Population initialization may hang if accuracy threshold cannot be reached with randomized scaling
- Files: `R/kernel.R` (lines 183-192 binary search loop)
- Trigger: Edge cases where the density distribution shape prevents reaching specified accuracy
- Workaround: Use very low accuracy values (e.g., 0.01 instead of 0.001); manually verify convergence

**Parameter Random Sampling Lacks Bounds:**
- Symptoms: Drawing random parameters from posterior can select iterations at extremes without validation
- Files: `R/params.R` (line 51: `sample(1:4000, 1)`)
- Trigger: When posterior has poor mixing or extreme values at tail iterations
- Workaround: Manually inspect sampled parameter sets for biological plausibility

## Security Considerations

**Unrestricted File Path Construction:**
- Risk: Data paths constructed from user input without validation; potential path traversal
- Files: `R/params.R` (lines 28-32)
- Current mitigation: Hardcoded directory structure limits exposure
- Recommendations: Validate `path` argument; use `normalizePath()` with `mustWork=TRUE` to prevent directory traversal

**Model File Dependency:**
- Risk: `getPars_sp` function depends on external RDS files that must exist; no validation of file integrity
- Files: `R/params.R` (lines 36-54)
- Current mitigation: `readRDS()` will error if file missing
- Recommendations: Add explicit `file.exists()` checks; add MD5 checksums for validation; handle corrupted files gracefully

## Performance Bottlenecks

**Inefficient Kernel Construction with outer():**
- Problem: `mkKernel` uses `outer()` with function application, creating full matrix for both P and F kernels; potentially 10,000+ elements
- Files: `R/kernel.R` (lines 107-117)
- Cause: Matrix dimensions scale with discretization granularity; outer product evaluates all combinations even where contributions are small
- Improvement path:
  1. Implement sparse matrix representation (Matrix package)
  2. Profile to identify which kernel regions have negligible values
  3. Consider approximate methods for far-field size classes

**Repeated Basal Area Calculations:**
- Problem: `mkKernel` recomputes competition metrics from size vectors on every call; no caching
- Files: `R/kernel.R` (lines 86-104)
- Cause: Competition calculations called fresh each time kernel is evaluated
- Improvement path: Accept pre-computed competition metrics as parameters; memoize if called repeatedly

**Population Initialization Convergence:**
- Problem: Binary search in `init_pop` uses random `runif()` sampling instead of systematic bisection
- Files: `R/kernel.R` (lines 182-192)
- Cause: Random search requires many iterations; statistical convergence rather than guaranteed convergence
- Improvement path: Implement deterministic binary search; use bracketing method (Brent's method or Newton-Raphson)

## Fragile Areas

**Vital Rates Functions - Assumption Violations:**
- Files: `R/vital_rates.R` (lines 27-48, 60-80, 91-113)
- Why fragile:
  - `vonBertalanffy_f`: Assumes linear combination of competition terms is appropriate; no validation of parameter signs
  - `survival_f`: Uses logistic with arbitrary intercept; symmetric around `psi`; no bounds checking on probability outputs
  - `ingrowth_f`: Combines exponential and truncated normal; discontinuities possible at boundaries
- Safe modification: Add assertions checking that survival_prob stays in [0,1]; add parameter validation for expected signs (e.g., `Beta < 0` for competition effect)
- Test coverage: No unit tests; only functional validation through simulation equilibria

**Kernel Discretization - Resolution Dependency:**
- Files: `R/kernel.R` (lines 146-149, 202-208)
- Why fragile:
  - Mesh point calculations depend on `h` (integration bin size) and hardcoded `Lmax`
  - Small changes to `h` can cascade to sensitivity analyses
  - No validation that mesh captures distribution appropriately
- Safe modification: Add warnings if population distribution extends beyond mesh; document `h` selection rationale
- Test coverage: No tests for discretization error or convergence in `h`

**Parameter Loading - Model Specification Coupling:**
- Files: `R/params.R` (lines 7-58)
- Why fragile:
  - Hardcoded model names and directory structure
  - If file naming convention changes, code breaks silently (path-based rather than metadata-based)
  - `time_truc` appears typo-prone (should be "truncated"?)
- Safe modification: Use configuration file for model-to-file mapping; add explicit validation of loaded parameters (check required keys)
- Test coverage: No validation that loaded parameters have expected structure

**Competition Metrics - Unit Consistency:**
- Files: `R/BasalArea_competition.R` (all functions)
- Why fragile:
  - Multiple unit conversions (mm → m, individual → per hectare) without clear documentation
  - Magic constants `1e-3`, `1e4` appear in calculations
  - No assertions on output ranges (basal area should be non-negative)
- Safe modification: Document all unit assumptions; add sanity checks on outputs (BA_comp >= 0, <= realistic maximum)
- Test coverage: No unit tests; only validated through simulation outcomes

## Scaling Limits

**Kernel Matrix Memory:**
- Current capacity: Discretization over ~100 size classes typical = 10,000 matrix elements per kernel
- Limit: 3-4 species × 4 kernels each = 12-16 matrices in memory simultaneously; dense matrices become limiting >1000 size classes
- Scaling path: Implement sparse matrix representation; consider separating P and F evaluations to reduce simultaneous memory use

**Population Initialization Loop:**
- Current capacity: Handles populations up to ~10,000 individuals with `1e4` samples
- Limit: Exponentially worse convergence with larger initial populations or finer discretization
- Scaling path: Implement adaptive sampling; use numerical integration instead of Monte Carlo for density estimation

**Parameter Space:**
- Current capacity: Posterior samples stored in memory (assuming 4000 iterations per model × 3 vital rates × multiple species)
- Limit: Large species sets or higher-dimensional posteriors will exhaust memory
- Scaling path: Stream from disk; implement efficient posterior summarization methods

## Dependencies at Risk

**truncnorm Package - Narrow Purpose:**
- Risk: Single use in `ingrowth_lk` (line 58-63 in `R/kernel.R`); could be replaced with qnorm + manual truncation
- Impact: Loss of package breaks ingrowth calculations; package maintenance critical
- Migration plan: Implement custom truncated normal using base functions: `pnorm()` for normalization, `qnorm()` for quantile

**RcppEigen Compilation:**
- Risk: C++ eigenvalue solver (`src/eigen.cpp`) compiled on load; compilation failures block entire package
- Impact: Installation fails on systems without C++ toolchain
- Migration plan: Pre-compile shared library; fall back to `eigen()` from base R (slower)

**tidyverse Ecosystem (dplyr, tidyr, purrr):**
- Risk: Soft dependency in `params.R`; functions use pipe operator and tidyverse functions but not declared in DESCRIPTION
- Impact: Silent failures if tidyverse not installed; inconsistent with rest of package which uses base R
- Migration plan: Declare explicit dependencies in DESCRIPTION; or rewrite `pars_to_list` to use base R only

## Missing Critical Features

**No Input Validation:**
- Problem: Functions accept parameter structures with no validation of required keys, numerical ranges, or data types
- Blocks: Difficult error messages when users pass malformed parameters; silent failures possible

**No Logging or Debugging Support:**
- Problem: No mechanism to trace execution flow or convergence diagnostics
- Blocks: Difficult to diagnose why `init_pop` or `mkKernel` fail in production; simulation results cannot be traced to parameter combinations

**No Sensitivity Analysis Framework:**
- Problem: Sensitivity analyses relegated to `simulations/` directory; no programmatic way to perturb parameters and track outcomes
- Blocks: Cannot easily explore parameter uncertainty; lambda elasticities not implemented

**No Documentation of Model Assumptions:**
- Problem: Von Bertalanffy parameterization, competition interaction form, recruitment model all lack formal justification
- Blocks: Users cannot assess whether model is appropriate for their system; difficult to extend or adapt

## Test Coverage Gaps

**Vital Rates Functions - Untested:**
- What's not tested: `vonBertalanffy_f`, `survival_f`, `ingrowth_f` have no unit tests
- Files: `R/vital_rates.R`
- Risk: Parameter changes silently produce biologically implausible outputs (e.g., negative growth, survival > 1); regressions in refactoring undetected
- Priority: High

**Kernel Construction - Untested:**
- What's not tested: `mkKernel` assembly, numerical integration accuracy, P + F = K identity
- Files: `R/kernel.R` (lines 78-122)
- Risk: Kernel construction errors propagate to lambda calculations; discretization artifacts undetected
- Priority: High

**Population Initialization - Untested:**
- What's not tested: `init_pop` convergence, edge cases (N=0, N=very large), accuracy threshold behavior
- Files: `R/kernel.R` (lines 137-209)
- Risk: Silent failures when accuracy cannot be achieved; population distributions not smooth as expected
- Priority: Medium

**Parameter Loading - Untested:**
- What's not tested: `getPars_sp` with missing files, malformed parameter structures, model name validation
- Files: `R/params.R`
- Risk: Silent failures when files missing; inconsistent parameter structures passed to vital rate functions
- Priority: Medium

**Basal Area Calculations - Untested:**
- What's not tested: Unit conversions, negative/zero input handling, cumulative competition calculations
- Files: `R/BasalArea_competition.R`
- Risk: Unit conversion errors accumulate in simulations; incorrect scaling in competition effects
- Priority: Medium

**Visualization - Untested:**
- What's not tested: `matrix.image` with extreme values, empty matrices, missing axes
- Files: `R/matrix_image.R`
- Risk: Plotting errors halt simulation workflows; undocumented parameter combinations break
- Priority: Low

---

*Concerns audit: 2026-02-24*
