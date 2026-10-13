# clusterMIO quickstart -- runs on simulated data, no external files or paths.
# install.packages(c("Rcpp", "RcppEigen"))
# remotes::install_github("MadhavShanky/cluster-mio", subdir = "r-implementation")
library(clusterMIO)

## ---- simulate K clusters, two of them deviating -----------------------------
set.seed(1)
K  <- 30          # clusters
nk <- 20          # observations per cluster
n  <- K * nk
cluster <- factor(rep(seq_len(K), each = nk))
X <- matrix(rnorm(n * 2), n, 2)              # two fixed-effect covariates
gamma_true <- numeric(K)
gamma_true[c(5, 22)] <- c(5, -5)             # clusters 5 (high) and 22 (low) deviate
y <- as.numeric(X %*% c(1, -1)) + gamma_true[as.integer(cluster)] + rnorm(n)

## ---- fit: select up to 2 deviating clusters ---------------------------------
fit <- scs_fit(X, y, cluster, lambda = 2)
print(fit)
fit$selected                                  # should recover "5" and "22"

## ---- asymmetric tail budgets: <=1 worst, <=1 best ---------------------------
fit_asym <- scs_fit(X, y, cluster, lambda = c(worst = 1, best = 1))
fit_asym$selected_pos; fit_asym$selected_neg

## ---- stability selection with the CPSS complementary-pairs bound ------------
stab <- stability(fit, B = 200L, pi_thr = 0.6, pairing = "complementary")
print(stab)

## ---- audit-budget profile and budget-path cluster ranking -------------------
print(audit_budget(fit))
print(cluster_ranking(fit))

## ---- predict the effect of a brand-new cluster (CART soft map) --------------
Xnew <- matrix(rnorm(5 * 2), 5, 2)
predict(fit, newdata = Xnew, type = "gamma")
