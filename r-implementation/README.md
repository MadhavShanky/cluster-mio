# Sparse Cluster Selection (SCS) — R/C++ implementation

An R/C++ implementation of sparse cluster selection in linear mixed models: an
$\ell_0$ (cardinality-constrained) selection of which cluster random effects are
held nonzero, with the rest set to exactly zero. This is the scalable companion to
the Julia prototype in the repository root (`clust_mio.jl`, `simul.jl`).

The solver profiles out the fixed effects once and exploits the disjoint support of
the cluster dummies, so each candidate swap has a closed-form gain that costs
$O(p^2)$ to evaluate, independent of the sample size. It depends only on a C++
linear-algebra library (RcppEigen) and requires no commercial optimizer.

## Layout

- `scs.R` — R API: `scs_fit()`, `stability()`, `audit_budget()`, `cluster_ranking()`.
- `src_scs.cpp` — the C++ local-search solver (RcppEigen).
- `sim/` — simulation study: data-generating processes (`dgp.R`, `dgp_hard.R`),
  drivers (`run_sims.R`, `run_hard.R`, `run_demo_phase.R`, `run_smoke_realistic.R`),
  the L0Learn head-to-head (`bench_l0learn.R`), and figure/table builders
  (`figures.R`, `tables.R`, `theme_scs.R`).
- Real-data analysis: `pull_sparcs_hf.R`, `tx_hf_analysis.R`, `add_inference.R`,
  `realdata_sensitivity.R`, `figures_realdata.R` (hospital length-of-stay profiling,
  New York SPARCS and Texas PUDF).
- Validation / benchmarks: `validate_cpp.R`, `verify_fwl.R`, `verify_localsearch.R`,
  `bench_speed.R`, `test_deliverables.R`.
- `*.sh` — array-job scripts for running the simulation grid on a cluster.

## Install

```r
install.packages(c("Rcpp", "RcppEigen"))
Rcpp::sourceCpp("src_scs.cpp")
source("scs.R")
```

## Quick start

```r
fit  <- scs_fit(y, X, cluster, lambda = 15)   # select up to 15 deviating clusters
stab <- stability(y, X, cluster, B = 100)     # bootstrap watch-list + expected-false-flag bound
aud  <- audit_budget(y, X, cluster)           # deviation captured vs number of units reviewed
```

## Reproducing the paper

The simulation and analysis scripts reproduce every figure and table. Note that
several scripts set an absolute working directory near the top (`setwd(...)`) and
read/write under `data/`, `sim/results/`, and `figs/`; adjust the path to your
checkout. Those generated directories are gitignored — the pull scripts
(`pull_sparcs_hf.R`, `tx_hf_analysis.R`) fetch the public discharge data, and
`sim/run_sims.R` regenerates `sim/results/`, after which `sim/figures.R` and
`sim/tables.R` build the outputs. Both hospital datasets are public; for the Texas
PUDF, follow its data-use-agreement terms.
