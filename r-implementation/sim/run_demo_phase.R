# run_demo_phase.R — END-TO-END DEMO of the redesigned sim pipeline + F1 "advantage" phase diagram.
# Proves: DGP -> MIO (turbo) + LMM-Gaussian (lme4) baseline -> tidy results -> ggplot phase tile.
# Coarse grid for speed; the full study scales this up via run_sims.R.
setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
source("scs.R"); source("sim/dgp.R"); source("sim/theme_scs.R")
suppressMessages({ library(lme4); library(dplyr) })

set.seed(20260619)
Ks   <- c(10, 30, 100)
pis  <- c(0.05, 0.20, 0.40, 0.60, 0.80)
R    <- 8L
nk <- 30L; p <- 5L; delta <- 3; rho_x <- 0.5; sigma <- 1

gamma_blup_lmm <- function(d) {
  df <- data.frame(y = d$y, d$X, cluster = d$cluster)
  fm <- suppressWarnings(suppressMessages(
    lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  re <- ranef(fm)$cluster[, 1]                       # BLUP random intercepts, ordered by level
  setNames(re, levels(d$cluster))[as.character(seq_len(d$K))]
}

rows <- list(); ix <- 0
for (K in Ks) for (pa in pis) {
  eM <- eL <- numeric(0)
  for (r in seq_len(R)) {
    d <- gen_scs(K, nk, p, regime = "sparse", pi_active = pa, delta = delta,
                 rho_x = rho_x, sigma = sigma, seed = 1000 * K + 10 * round(100*pa) + r)
    lam <- length(d$active)                          # oracle cardinality: isolates estimation quality
    fit <- scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r)
    g_mio <- fit$gamma[as.character(seq_len(K))]; g_mio[is.na(g_mio)] <- 0
    g_lmm <- gamma_blup_lmm(d)
    eM <- c(eM, sqrt(sum((g_mio        - d$gamma)^2)))
    eL <- c(eL, sqrt(sum((as.numeric(g_lmm) - d$gamma)^2)))
  }
  ix <- ix + 1
  rows[[ix]] <- data.frame(K = K, pi = pa,
                           err_mio = mean(eM), err_lmm = mean(eL),
                           advantage = log2(mean(eL) / mean(eM)))
}
res <- bind_rows(rows)
dir.create("sim/figs", showWarnings = FALSE)
write.csv(res, "sim/figs/demo_phase.csv", row.names = FALSE)
cat("=== demo phase results (advantage = log2(err_LMM/err_MIO); >0 => L0 wins) ===\n")
print(res, row.names = FALSE)

lim <- c(-1, 1) * max(abs(res$advantage))
pF1 <- scs_phase_tile(res, x = "pi", y = "K", value = "advantage", diverging = TRUE, limit = lim,
                      title = "When does L0 cluster-selection beat shrinkage?",
                      subtitle = sprintf("gamma-recovery advantage, sparse regime (delta=%g, n_k=%d, p=%d), R=%d reps", delta, nk, p, R),
                      xlab = "fraction of clusters deviating (pi)",
                      value_name = "L0 advantage\nlog2(err_LMM/err_MIO)")
scs_save(pF1, "sim/figs/F1_phase_demo.png", w = 6.8, h = 4.4)
cat("\nSaved sim/figs/F1_phase_demo.png\n")
