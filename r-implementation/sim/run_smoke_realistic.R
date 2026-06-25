# run_smoke_realistic.R -- realistic-regime smoke tests: "where does L0 cluster-selection beat shrinkage?"
# Lambda tuned by WITHIN-cluster K-fold CV (row-wise holdout; every cluster appears in train & test).
#   Rationale: in the DGP gamma_k _|_ X, so CLUSTER-blocked CV cannot tune
#   lambda (new-cluster MSE is flat in lambda and ~tied L0 vs LMM by construction). We therefore tune
#   within-cluster and report the agreed 3-metric panel below. Cluster-blocked/new-cluster MSE is the
#   tie metric -> deferred to the full study (run_sims.R), not computed here.
# Panel (no single headline): selection F1/FDR, gamma-recovery error, within-cluster predictive MSE.
# Baseline = LMM-Gaussian BLUP (lme4). Set SMOKE_QUICK=TRUE for a 1-cell validation run.
suppressWarnings(suppressMessages({
  setwd("C:/Users/tomch/AIProjects/research/MIO LMM")
  source("scs.R"); source("sim/dgp.R"); library(lme4); library(dplyr)
}))
QUICK <- isTRUE(as.logical(Sys.getenv("SMOKE_QUICK", "FALSE")))

## baseline LMM: BLUP random intercepts + fixed effects (mirrors run_demo_phase.R) -----------------
lmm_fit <- function(X, y, cluster) {
  df <- data.frame(y = y, X, cluster = cluster)
  fm <- suppressWarnings(suppressMessages(
    lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  re <- ranef(fm)$cluster[, 1]
  list(blup = setNames(re, rownames(ranef(fm)$cluster)), fixef = lme4::fixef(fm))
}

## within-cluster K-fold CV: choose lambda minimizing row-holdout MSE; also LMM row-holdout MSE --------
cv_within <- function(d, grid, nfold = 5L, seed = 1L) {
  set.seed(seed); n <- length(d$y); fold <- sample(rep_len(seq_len(nfold), n))
  Xf <- cbind(1, d$X)
  cvm <- numeric(length(grid)); cnt <- numeric(length(grid)); lmm_err <- numeric(0)
  for (f in seq_len(nfold)) {
    tr <- fold != f; te <- !tr; cl_tr <- droplevels(d$cluster[tr])
    for (gi in seq_along(grid)) {
      lam_f <- min(grid[gi], nlevels(cl_tr) - 1L)
      fit <- tryCatch(scs_fit(d$X[tr, , drop = FALSE], d$y[tr], cl_tr, lambda = lam_f,
                              n_restart = 2L, seed = 1L), error = function(e) NULL)
      if (is.null(fit)) next
      g <- fit$gamma[as.character(d$cluster[te])]; g[is.na(g)] <- 0
      pred <- as.numeric(Xf[te, , drop = FALSE] %*% fit$fixef) + g
      cvm[gi] <- cvm[gi] + mean((d$y[te] - pred)^2); cnt[gi] <- cnt[gi] + 1
    }
    lf <- tryCatch(lmm_fit(d$X[tr, , drop = FALSE], d$y[tr], cl_tr), error = function(e) NULL)
    if (!is.null(lf)) {
      gb <- lf$blup[as.character(d$cluster[te])]; gb[is.na(gb)] <- 0
      pm <- as.numeric(Xf[te, , drop = FALSE] %*% lf$fixef) + gb
      lmm_err <- c(lmm_err, mean((d$y[te] - pm)^2))
    }
  }
  cvm <- cvm / pmax(cnt, 1)
  list(lambda = grid[which.min(cvm)], mse_l0 = min(cvm),
       mse_lmm = if (length(lmm_err)) mean(lmm_err) else NA_real_)
}

## one cell = R reps at fixed (K, pi, delta) ---------------------------------------------------------
eval_cell <- function(K, nk, p, pi_a, delta, rho_x, sigma, R, nfold, grid_fn) {
  acc <- vector("list", R)
  for (r in seq_len(R)) {
    d <- gen_scs(K, nk, p, regime = "sparse", pi_active = pi_a, delta = delta,
                 rho_x = rho_x, sigma = sigma,
                 seed = 7000000L + 10000L * K + 100L * round(100 * pi_a) + 10L * delta + r)
    cv  <- cv_within(d, grid_fn(K), nfold, seed = r)
    lam <- cv$lambda
    fitL <- scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r)
    sL  <- sel_scores(match(fitL$selected, levels(d$cluster)), d$active, K)
    gL  <- fitL$gamma[as.character(seq_len(K))]; gL[is.na(gL)] <- 0
    lm0 <- tryCatch(lmm_fit(d$X, d$y, d$cluster), error = function(e) NULL)
    if (is.null(lm0)) next
    blup <- lm0$blup[as.character(seq_len(K))]; blup[is.na(blup)] <- 0
    selM <- order(abs(blup), decreasing = TRUE)[seq_len(max(1L, lam))]  # LMM: flag top-lambda by |BLUP|
    sM   <- sel_scores(selM, d$active, K)
    acc[[r]] <- data.frame(
      K = K, pi = pi_a, delta = delta, lambda = lam, Ka = length(d$active),
      f1_l0 = sL$f1, fdr_l0 = sL$fdr, f1_lmm = sM$f1, fdr_lmm = sM$fdr,
      gerr_l0 = sqrt(sum((gL - d$gamma)^2)), gerr_lmm = sqrt(sum((blup - d$gamma)^2)),
      wmse_l0 = cv$mse_l0, wmse_lmm = cv$mse_lmm)
  }
  bind_rows(acc)
}

## grid of candidate budgets for CV (relative to K, no oracle peek) ----------------------------------
grid_fn <- function(K) unique(pmax(0L, as.integer(round(K * c(0, 0.02, 0.05, 0.10, 0.20, 0.35)))))

if (QUICK) { cells <- expand.grid(K = 100L, pi = 0.10, delta = 3L); R <- 2L; nfold <- 3L
} else      { cells <- expand.grid(K = c(100L, 200L), pi = c(0.05, 0.15), delta = c(2L, 3L)); R <- 10L; nfold <- 5L }
nk <- 30L; p <- 5L; rho_x <- 0.5; sigma <- 1

out <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  cat(sprintf("[cell %d/%d] K=%d pi=%.2f delta=%d  ", i, nrow(cells), cells$K[i], cells$pi[i], cells$delta[i]))
  t0 <- Sys.time()
  out[[i]] <- eval_cell(cells$K[i], nk, p, cells$pi[i], cells$delta[i], rho_x, sigma, R, nfold, grid_fn)
  cat(sprintf("(%.1fs)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res <- bind_rows(out)
summ <- res %>% group_by(K, pi, delta) %>%
  summarise(across(c(lambda, Ka, f1_l0, f1_lmm, fdr_l0, fdr_lmm,
                     gerr_l0, gerr_lmm, wmse_l0, wmse_lmm), mean), .groups = "drop") %>%
  mutate(f1_gain = f1_l0 - f1_lmm, gerr_ratio = gerr_lmm / gerr_l0, wmse_gain = wmse_lmm - wmse_l0)

dir.create("sim", showWarnings = FALSE)
write.csv(res,  if (QUICK) "sim/smoke_quick_raw.csv"     else "sim/smoke_realistic_raw.csv",     row.names = FALSE)
write.csv(summ, if (QUICK) "sim/smoke_quick_summary.csv" else "sim/smoke_realistic_summary.csv", row.names = FALSE)
cat("\n=== realistic smoke summary  (f1_gain>0, gerr_ratio>1, wmse_gain>0  => L0 wins) ===\n")
print(as.data.frame(summ), row.names = FALSE, digits = 3)
