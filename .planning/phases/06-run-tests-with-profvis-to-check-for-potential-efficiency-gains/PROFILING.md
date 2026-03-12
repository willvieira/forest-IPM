# Phase 6 — Profiling Results

**Date:** 2026-03-11
**Package version:** forestIPM (development, post-06-01 tidyverse revamp)
**Benchmark:** QUERUB single-species, 5 timesteps via project() (profvis) + component microbenchmark (10 reps each)

---

## Mesh Size

- QUERUB mesh: **2235 points** (bin_width = 1 mm, DBH range 127 mm to 2362 mm, Lmax = 2361 mm)
- Kernel matrix dimensions: **2235 × 2235 = 4,995,225 cells** per kernel build
- Kernel builds per 100-year projection: 100 (one per timestep)
- Estimated cost for 100-year run: ~76 seconds (100 × 761 ms median mkKernel time)

---

## Top Hotspots (profvis, 5-year single-species project())

Total profiled time: ~3075 ms (615 samples × 5 ms interval).

| Rank | Function | Inclusive % | Notes |
|------|----------|-------------|-------|
| 1 | `project()` | 86.0% | Entry point — nearly all time is inside mkKernel() |
| 2 | `mkKernel()` | 85.9% | Kernel assembly: almost the entirety of project() cost |
| 3 | `outer()` | 73.8% | Both P and F matrix outer() calls combined |
| 4 | `FUN` (outer callback) | 63.3% | Inner function dispatched per cell by outer() |
| 5 | `truncnorm::dtruncnorm` | 46.0% | **Single largest leaf hotspot** — called inside ingrowth_lk per cell |
| 6 | `vonBertalanffy_lk` | 16.7% | Growth probability — called inside P_xEC per cell |
| 7 | `dnorm` | 13.2% | R's base dnorm inside vonBertalanffy_lk |
| 8 | `<GC>` | 15.3% | Garbage collection — 4.99M temporary cell values per outer() call |
| 9 | `size_to_BAcomp` | 9.6% | Competition metrics — called 2× per mkKernel() |

---

## Call Stack Breakdown

The dominant hot path is:

```
project()
  └─ mkKernel()                        85.9% inclusive
       ├─ outer(P_xEC, ...)             ~197ms median  [P matrix]
       │   └─ P_xEC() per cell
       │       ├─ vonBertalanffy_lk()   16.7% inclusive
       │       │   └─ dnorm()           13.2%
       │       └─ survival_f()          (fast — not visible in profvis)
       ├─ outer(ingrowth_lk, ...)        ~491ms median  [F matrix — dominant]
       │   └─ ingrowth_lk() per cell
       │       └─ truncnorm::dtruncnorm  46.0% inclusive  <-- biggest leaf hotspot
       └─ size_to_BAcomp / size_to_BAplot  ~77ms (called 4× total)
```

The F matrix (`outer` + `ingrowth_lk` + `truncnorm::dtruncnorm`) accounts for ~64% of mkKernel() time. The truncated normal distribution evaluation inside `ingrowth_lk` is called once per matrix cell — that is 4,995,225 R function calls per timestep.

---

## mkKernel() Component Microbenchmark (10 reps, QUERUB single-species)

| Component | Median (ms) | % of Full mkKernel |
|-----------|-------------|-------------------|
| `mkKernel()` full | 761 ms | 100% |
| F matrix only (`outer` + `ingrowth_lk` incl. `truncnorm`) | 491 ms | 64% |
| P matrix only (`outer` + `P_xEC` incl. `dnorm`) | 197 ms | 26% |
| BA metrics only (4× calls: `size_to_BAcomp` + `size_to_BAplot`) | 77 ms | 10% |
| **Estimated 100-year run cost** | **~76 s** | — |

Note: BA metrics cost is measured separately (with the BA calls pre-computed for P/F benchmarks) — the sum P+F+BA ≈ 765ms matches the full mkKernel() median of 761ms, confirming the decomposition is accurate.

---

## Optimization Targets

### Target 1: Vectorize outer() calls in mkKernel() [HIGHEST PRIORITY — RECOMMENDED]

**Evidence:** `outer()` accounts for 73.8% of total runtime. Both `P_xEC` and `ingrowth_lk` are called ~5M times per timestep as scalar R functions via `outer()`. This is n² R function dispatch overhead — the most expensive pattern in R.

**Strategy:** Replace `outer(meshpts, meshpts, FUN, ...)` with vectorized computation:
1. Use `expand.grid()` or `rep()`/vectorized operations to build the full input matrix at once.
2. Call `P_xEC` / `ingrowth_lk` on full vectors (n² element vectors) instead of n² individual calls.
3. `matrix()` the result back to n × n shape.

Since both `vonBertalanffy_lk` / `survival_f` / `ingrowth_lk` already operate element-wise (no internal loops), this is a pure call overhead fix — no math changes required.

**Expected gain:** 20–40% reduction in mkKernel() time. The GC overhead (15.3%) from 5M temporary objects per outer() call should also drop substantially.

**Risk:** LOW — math is identical; regression tests at 1e-10 tolerance in `test-regression-baselines.R` guard against numeric drift.

**Implementation:** `R/kernel.R` only (mkKernel function). No API changes.

---

### Target 2: Replace truncnorm::dtruncnorm with native R truncated normal [HIGH PRIORITY]

**Evidence:** `truncnorm::dtruncnorm` is the single largest leaf hotspot at 46% inclusive time. It is called once per cell of the F matrix (4.99M calls per timestep). The `truncnorm` package uses C, but the function has per-call R dispatch overhead that accumulates at this scale.

**Strategy:** Replace with a direct R implementation using `dnorm` + `pnorm` (the definition of a truncated normal PDF):

```r
# Replace truncnorm::dtruncnorm(x, a=127, b=Inf, mean=mu, sd=sig) with:
dnorm(x, mean = mu, sd = sig) / (1 - pnorm(127, mean = mu, sd = sig))
```

This uses only base R `dnorm` and `pnorm` which vectorize natively — critical when called with n² element vectors after Target 1's vectorization is applied.

**Expected gain:** This target compounds with Target 1. If Target 1 passes n²-element vectors to `ingrowth_lk`, then `dtruncnorm` will also receive n²-element vector input. Base R `dnorm`/`pnorm` vectorize efficiently; `truncnorm::dtruncnorm` has additional overhead at the R-to-C boundary per call. Combined with Target 1, this could eliminate most of the 46% hotspot.

**Risk:** LOW — math is numerically equivalent. The truncation lower bound (a=127) and upper bound (b=Inf) are fixed, simplifying the formula. Must verify output matches original to 1e-10 tolerance using the existing regression fixtures.

**Implementation:** `R/kernel.R` — `ingrowth_lk()` function only (2-line change).

---

### Target 3: Cache competition metrics across timesteps [CONDITIONAL]

**Evidence:** `size_to_BAcomp` / `size_to_BAplot` account for ~10% of mkKernel() time (~77ms per call, 4 calls per mkKernel). For a 100-year single-species run, total BA cost ≈ 7.7 seconds.

**Strategy:** Cache BA metrics when the size distribution (Nvec) has not changed significantly between timesteps. Since the population distribution does change every timestep (matrix multiplication updates Nvec), this caching is only useful if the change is below a tolerance threshold — or if the cache is keyed to the exact Nvec vector.

**Expected gain:** 10% reduction in mkKernel() time if applied.

**Risk:** MEDIUM — must not cache inter-species metrics when competitor Nvec changes (stale Nvec_inter bug risk in multi-species runs). The implementation complexity and the relatively small gain (10%) make this lower priority.

**Recommendation:** Apply only after Targets 1 and 2. Evaluate residual profile after those optimizations to see if BA metrics remain significant.

**Implementation:** `R/kernel.R` (mkKernel) + potentially `R/project.R` (pass cached values).

---

### Target 4: C++ migration of P_xEC or ingrowth_lk [CONDITIONAL — LOWER PRIORITY]

**Evidence:** After Targets 1 and 2, the main remaining per-element cost would be the arithmetic inside `vonBertalanffy_lk` (von Bertalanffy growth formula) and `survival_f` (logistic survival). These are scalar-valued but if vectorized, R's JIT compiler handles them reasonably well.

**Strategy:** Move `P_xEC` or `ingrowth_lk` to C++ via Rcpp. The package already depends on RcppEigen, so adding Rcpp functions is low dependency overhead.

**Risk:** MEDIUM — requires `Rcpp::compileAttributes()` and Xcode license acceptance (`sudo xcodebuild -license accept` — see STATE.md blocker). Do not attempt until the Xcode blocker is resolved.

**Recommendation:** Defer until after Targets 1–3 are applied and benchmarked. If the 100-year run time after Targets 1 and 2 drops from ~76s to under 20s, C++ migration may not be necessary.

---

## Unexpected Finding: F Matrix Dominates, Not P Matrix

The plan anticipated that the P matrix (`outer` + `P_xEC`) would be the dominant cost. Profiling shows the **F matrix (`outer` + `ingrowth_lk` + `truncnorm::dtruncnorm`) is 2.5× more expensive than the P matrix**. This changes the priority order: Target 2 (replacing `truncnorm::dtruncnorm`) should be addressed alongside Target 1 (vectorizing `outer`), not treated as secondary.

---

## Artifacts

- `profvis_output.html` — single-species (QUERUB, 5 yr) flamegraph (open in browser)
- `profvis_multi_output.html` — 2-species (QUERUB + QUEPRI, 3 yr) flamegraph
- `profile_benchmark.R` — script to reproduce all results

---

## Summary Recommendation

Apply **Target 1 (vectorize outer())** and **Target 2 (replace truncnorm::dtruncnorm)** together in Plan 03. These two changes address the 73.8% outer() overhead and the 46% truncnorm hotspot. Since Target 2 depends on having n²-element vector input (which Target 1 provides), they must be implemented together. Combined expected gain: **50–70% reduction** in mkKernel() time, bringing a 100-year single-species run from ~76 seconds to an estimated 23–38 seconds. Target 3 (BA caching) can be applied as a secondary optimization if BA metrics remain visible after Targets 1+2. Target 4 (C++) should be deferred until the Xcode license blocker is resolved.
