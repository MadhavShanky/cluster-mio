# clusterMIO — Sparse Cluster Selection for Linear Mixed Models

An R/C++ package for **sparse cluster selection** in linear mixed models: an
$\ell_0$ (cardinality-constrained) selection of which cluster random effects are
held nonzero, with the rest set to exactly zero. The scalable companion to the
Julia prototype in the repository root (`clust_mio.jl`, `simul.jl`).

The solver profiles out the fixed effects once and exploits the disjoint support
of the cluster dummies, so each candidate swap has a closed-form gain that costs
$O(p^2)$ to evaluate, independent of the sample size. It depends only on a C++
linear-algebra library (RcppEigen) and requires no commercial optimizer.

## Install

```r
install.packages(c("Rcpp", "RcppEigen"))   # build-time dependencies
# install.packages("remotes")
remotes::install_github("MadhavShanky/cluster-mio", subdir = "r-implementation")
```

The package lives in the `r-implementation/` subdirectory of the repo, hence the
`subdir` argument. A C++ toolchain is required to compile from source (Rtools on
Windows, the standard build tools on macOS/Linux). Optional: `rpart` (new-cluster
prediction), `lme4` (the BLUP comparison in the analysis scripts).

## Quick start

```r
library(clusterMIO)

set.seed(1)
K <- 30; nk <- 20; n <- K * nk
cluster <- factor(rep(seq_len(K), each = nk))
X <- matrix(rnorm(n * 2), n, 2)
g <- numeric(K); g[c(5, 22)] <- c(5, -5)          # two deviating clusters
y <- as.numeric(X %*% c(1, -1)) + g[as.integer(cluster)] + rnorm(n)

fit <- scs_fit(X, y, cluster, lambda = 2)          # select up to 2 deviating clusters
fit$selected                                       # -> "5", "22"

stability(fit, B = 200L, pairing = "complementary")# CPSS watch-list + false-flag bound
audit_budget(fit)                                  # deviation captured vs. units reviewed
cluster_ranking(fit)                               # budget-path importance ordering
predict(fit, newdata = matrix(rnorm(4), 2, 2))     # effect of a brand-new cluster (CART map)
```

A runnable version is installed at
`system.file("examples", "quickstart.R", package = "clusterMIO")`.

## API

| Function | Purpose |
|---|---|
| `scs_fit()` | Fit the L0 cluster-selection model (symmetric `lambda`, or asymmetric `c(worst, best)` tail budgets). |
| `predict()` | New-cluster prediction via a supervised CART map from cluster-level features to the cluster effect. |
| `stability()` | Subsampling stability selection: Meinshausen–Bühlmann or Shah–Samworth complementary-pairs (CPSS), with the expected-false-flag bound. |
| `scs_path()` | Solve the symmetric problem along the budget path `lambda = 0..lambda_max`. |
| `audit_budget()` | Decision-theoretic budget profile: deviation captured vs. number of clusters reviewed. |
| `cluster_ranking()` | LARS-like importance ranking by budget at which each cluster first enters the active set. |
| `scs_eps_sensitivity()` | Sensitivity of the selected set to the relative-ridge scale `eps`. |

## Reproducing the paper

The simulation and real-data analysis scripts live at the top level of
`r-implementation/` and in `sim/`. They are kept in the repository but are
**excluded from the installed package** (`.Rbuildignore`): several set an absolute
working directory and read/write under `data/`, `sim/results/`, and `figs/`, so
they predate the package and load the source directly — adjust the paths to your
checkout before running. Both hospital datasets (New York SPARCS, Texas PUDF) are
public; for the Texas PUDF follow its data-use-agreement terms.

See `PERFORMANCE.md` for the performance roadmap (profiled hot spots and the
optimizations deferred behind a selected-set regression gate).

## Reference

A distribution-free mixed-integer optimization approach to hierarchical modelling
of clustered and longitudinal data, arXiv:2302.03157.
