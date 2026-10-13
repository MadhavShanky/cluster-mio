# greedy_vs_brute.R -- Supplement S2 (Table S1): greedy forward / marginal |z| / SCS against brute force.
# Objective: unpenalized refit loss R(S) = y'M y - gain(S), gain(S) = r_S' H_SS^{-1} r_S with
# H = D'M D, r = D'M y, M the exact projection off [1, X]. Part A: brute force on enumerable instances.
# Part B: K = 100, 300 against SCS (no brute force), incl. designs built to stress greedy.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM")); source("scs.R"); source("sim/dgp.R")
}))
reduce <- function(X, y, cl) {
  Xf <- cbind(1, X); XtXi <- solve(crossprod(Xf)); G <- rowsum(Xf, cl); nk <- as.numeric(table(cl))
  Xty <- crossprod(Xf, y)
  list(H = diag(nk) - G %*% XtXi %*% t(G), r = as.numeric(rowsum(y, cl) - G %*% XtXi %*% Xty),
       K = length(nk), lev = levels(cl))
}
gain <- function(R, S) if (!length(S)) 0 else { S <- as.integer(S); sum(R$r[S] * solve(R$H[S, S, drop = FALSE], R$r[S])) }
greedy <- function(R, lam) { S <- integer(0)
  for (s in seq_len(lam)) { cand <- setdiff(seq_len(R$K), S)
    g <- vapply(cand, function(j) gain(R, c(S, j)), 0); S <- c(S, cand[which.max(g)]) }
  S }
topz <- function(R, lam) order(R$r^2 / diag(R$H), decreasing = TRUE)[seq_len(lam)]
brute <- function(R, lam) { cmb <- combn(R$K, lam); g <- apply(cmb, 2, function(S) gain(R, S)); cmb[, which.max(g)] }
scs_sel <- function(d, lam) match(scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = 1L)$selected, levels(d$cluster))
# stress design: add a cluster-level covariate (constant within cluster) correlated with the active set
add_cluster_cov <- function(d, rho = 0.7) {
  z <- as.numeric(d$gamma != 0); w <- rho * scale(z) + sqrt(1 - rho^2) * rnorm(length(z))
  d$X <- cbind(d$X, w[as.integer(d$cluster)]); d }

designs <- list(
  base      = list(rho_x = 0.5, unequal = FALSE, delta = 2, clcov = FALSE),
  weak      = list(rho_x = 0.5, unequal = FALSE, delta = 0.75, clcov = FALSE),
  unequal   = list(rho_x = 0.5, unequal = TRUE,  delta = 1, clcov = FALSE),
  highcorr  = list(rho_x = 0.95, unequal = TRUE, delta = 1, clcov = FALSE),
  clustercv = list(rho_x = 0.5, unequal = TRUE,  delta = 1, clcov = TRUE))
set.seed(20261013)
rowsA <- list()
for (dn in names(designs)) for (K in c(14L, 18L)) for (lam in c(3L, 4L)) for (rep in 1:25) {
  ds <- designs[[dn]]
  d <- gen_scs(K = K, nk = 20L, p = 5L, regime = "sparse", pi_active = lam / K, delta = ds$delta,
               rho_x = ds$rho_x, unequal_nk = ds$unequal)
  if (ds$clcov) d <- add_cluster_cov(d)
  R <- reduce(d$X, d$y, d$cluster); gb <- gain(R, brute(R, lam))
  gg <- gain(R, greedy(R, lam)); gz <- gain(R, topz(R, lam)); gs <- gain(R, scs_sel(d, lam))
  rowsA[[length(rowsA) + 1]] <- data.frame(design = dn, K, lam, gap_greedy = gb - gg, gap_topz = gb - gz, gap_scs = gb - gs)
}
A <- do.call(rbind, rowsA); tol <- 1e-8
cat("Part A (brute force): share of instances with suboptimal objective, and max gain gap\n")
print(aggregate(cbind(greedy = gap_greedy > tol, topz = gap_topz > tol, scs = gap_scs > tol) ~ design, A, mean), digits = 3)
print(aggregate(cbind(greedy = gap_greedy, topz = gap_topz, scs = gap_scs) ~ design, A, max), digits = 3)

rowsB <- list()
for (dn in c("base", "unequal", "highcorr", "clustercv")) for (K in c(100L, 300L)) for (rep in 1:10) {
  ds <- designs[[dn]]; lam <- as.integer(round(0.1 * K))
  d <- gen_scs(K = K, nk = 20L, p = 5L, regime = "sparse", pi_active = 0.1, delta = ds$delta,
               rho_x = ds$rho_x, unequal_nk = ds$unequal)
  if (ds$clcov) d <- add_cluster_cov(d)
  R <- reduce(d$X, d$y, d$cluster)
  Sg <- greedy(R, lam); Sz <- topz(R, lam); Ss <- scs_sel(d, lam)
  rowsB[[length(rowsB) + 1]] <- data.frame(design = dn, K, gain_scs_minus_greedy = gain(R, Ss) - gain(R, Sg),
    gain_scs_minus_topz = gain(R, Ss) - gain(R, Sz), same_greedy = setequal(Sg, Ss), same_topz = setequal(Sz, Ss))
}
B <- do.call(rbind, rowsB)
cat("\nPart B (K=100,300; lambda=0.1K): SCS vs greedy / top-|z| (gain diff > 0 means SCS better)\n")
print(aggregate(cbind(same_greedy, same_topz, scs_better_greedy = gain_scs_minus_greedy > 1e-8,
                      greedy_better_scs = gain_scs_minus_greedy < -1e-8) ~ design + K, B, mean), digits = 3)
saveRDS(list(A = A, B = B), "reviews/v2/greedy_sim.rds")
