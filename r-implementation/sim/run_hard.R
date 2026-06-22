# run_hard.R -- driver for the harder-regime robustness study.
# Reuses scs.R + lme4 + sel_scores(); same CV/stability machinery as run_sims.R but on the gen_hard()
# regimes. Sharding: env SIM_CELLS="start-end" (cid slice) and SIM_REPS_RANGE="a-b" (rep slice), exactly
# as run_sims.R; outputs sim/results_hard/cell_<id>[<repstart>].csv.
#
# Regimes / cells (all K in {100,300}, pi=0.1 where sparse, R=100):
#   boundary : sparse, delta in {0.25,0.5,0.75,1}            -- weak signals near detection boundary
#   nullsd   : sparse, delta=2, null_sd in {0.10,0.25}       -- "true zero" not clean
#   dense_t  : dense t_3 effects, icc in {0.3,0.5,0.7}        -- heavy-tailed dense (misspecification stress)
#   sizecorr : sparse, delta in {1,2}, size_effect=TRUE       -- deviating clusters smaller
#   misspec  : sparse, delta=2, drop_x=TRUE                   -- omitted cluster-correlated covariate
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp_hard.R"); library(lme4); library(dplyr)
}))
QUICK <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
RNGkind("Mersenne-Twister")

lmm_blup <- function(X, y, cluster) {
  df <- data.frame(y = y, X, cluster = cluster)
  fm <- suppressWarnings(suppressMessages(lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  v <- ranef(fm)$cluster[, 1]; setNames(v, rownames(ranef(fm)$cluster))
}

# one replicate for a hard-regime cell
run_rep_hard <- function(cell, r) {
  d <- gen_hard(cell$K, cell$nk, cell$p, regime = cell$regime, pi_active = cell$pi, delta = cell$delta,
                icc = cell$icc, rho_x = cell$rho_x, err = cell$err, nu = cell$nu,
                null_sd = cell$null_sd, size_effect = cell$size_effect,
                seed = as.integer((1e6 * cell$cid) %% 2e9) + r)
  K <- d$K
  Xfit <- if (isTRUE(cell$drop_x)) d$X[, -d$p, drop = FALSE] else d$X   # misspecified FE: omit last covariate
  grid <- unique(pmax(0L, as.integer(round(K * c(0, 0.02, 0.05, 0.10, 0.20, 0.35)))))
  o <- list(rep = r)

  # CV over budget (within-cluster, simplified): pick lambda minimizing within-cluster holdout MSE
  set.seed(r); n <- length(d$y); fold <- sample(rep_len(seq_len(5L), n)); Xf1 <- cbind(1, Xfit)
  cvm <- numeric(length(grid)); cnt <- numeric(length(grid))
  for (f in 1:5) {
    tr <- fold != f; te <- !tr; cl_tr <- droplevels(d$cluster[tr])
    for (gi in seq_along(grid)) {
      lam <- min(grid[gi], nlevels(cl_tr) - 1L)
      fit <- tryCatch(scs_fit(Xfit[tr, , drop = FALSE], d$y[tr], cl_tr, lambda = lam, n_restart = 2L, seed = 1L),
                      error = function(e) NULL)
      if (is.null(fit)) next
      g <- fit$gamma[as.character(d$cluster[te])]; g[is.na(g)] <- 0
      cvm[gi] <- cvm[gi] + mean((d$y[te] - (as.numeric(Xf1[te, , drop = FALSE] %*% fit$fixef) + g))^2)
      cnt[gi] <- cnt[gi] + 1
    }
  }
  cvm <- cvm / pmax(cnt, 1); cvm[cnt == 0] <- Inf        # (#6) never select an all-failed budget
  if (all(!is.finite(cvm))) stop("all CV fits failed")
  lam <- grid[which.min(cvm)]

  t0 <- Sys.time()                                        # (#7) timing of final fit only (CV excluded)
  fitL <- scs_fit(Xfit, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r)
  o$runtime_fit_l0 <- as.numeric(difftime(Sys.time(), t0, units = "secs")); o$lambda <- lam
  gL <- fitL$gamma[as.character(seq_len(K))]; gL[is.na(gL)] <- 0
  lm_ok <- TRUE                                           # (#8) record LMM failures
  blup <- tryCatch({ v <- lmm_blup(Xfit, d$y, d$cluster)[as.character(seq_len(K))]; v[is.na(v)] <- 0; v },
                   error = function(e) { lm_ok <<- FALSE; rep(0, K) })
  o$fail_lmm <- !lm_ok

  # estimation error vs the STRUCTURAL gamma (always reported)
  o$gerr_l0 <- sqrt(sum((gL - d$gamma)^2)); o$gerr_lmm <- sqrt(sum((blup - d$gamma)^2))
  # (#4) under misspecification the fitted cluster estimand absorbs the omitted covariate's between-cluster
  # part beta_p * mk[,p]; report error vs that contaminated target too, so the metric is interpretable.
  if (isTRUE(cell$drop_x)) {
    gmm <- d$gamma + d$beta[d$p] * d$mk[, d$p]
    o$gerr_l0_eff <- sqrt(sum((gL - gmm)^2)); o$gerr_lmm_eff <- sqrt(sum((blup - gmm)^2))
  }

  if (cell$regime == "sparse") {
    # selection vs TRUE active set (the +/- delta clusters; null_sd clusters are NOT active)
    o <- c(o, with(sel_scores(match(fitL$selected, levels(d$cluster)), d$active, K),
                   list(f1_l0 = f1, fdr_l0 = fdr, tpr_l0 = recall)))
    selM <- if (lam > 0) order(abs(blup), decreasing = TRUE)[seq_len(lam)] else integer(0)
    o <- c(o, with(sel_scores(selM, d$active, K), list(f1_lmm = f1, fdr_lmm = fdr, tpr_lmm = recall)))
    st <- tryCatch(stability(fitL, B = if (QUICK) 20L else 50L, pi_thr = 0.8, seed = r), error = function(e) NULL)
    if (!is.null(st)) o <- c(o, with(sel_scores(match(st$flagged, levels(d$cluster)), d$active, K),
                                     list(f1_stab = f1, fdr_stab = fdr, n_stab = length(st$flagged))),
                             list(mb_bound = st$expected_false_flags_bound))
  } else {
    # dense_t misspecification check: report which method recovers dense effects better (no selection metric)
    o$icc_true <- cell$icc
  }
  as.data.frame(o, stringsAsFactors = FALSE)
}

## ---- hard-regime cell grid ----------------------------------------------------------------------
ref <- list(nk = 30L, p = 5L, rho_x = 0.5, err = "normal", nu = 3, icc = NA_real_,
            null_sd = 0, size_effect = FALSE, drop_x = FALSE, pi = 0.1)
mk <- function(...) { a <- list(...); modifyList(c(list(K = 100L, regime = "sparse", delta = 2), ref), a) }

cells <- do.call(rbind, lapply(list(
  # boundary signals
  lapply(c(100L, 300L), function(K) lapply(c(0.25, 0.5, 0.75, 1.0), function(de) mk(K = K, delta = de))),
  # small-nonzero nulls
  lapply(c(100L, 300L), function(K) lapply(c(0.10, 0.25), function(ns) mk(K = K, delta = 2, null_sd = ns))),
  # heavy-tailed dense
  lapply(c(100L, 300L), function(K) lapply(c(0.3, 0.5, 0.7), function(ic) mk(K = K, regime = "dense_t", icc = ic, pi = NA_real_, delta = NA_real_))),
  # size/effect coupling
  lapply(c(100L, 300L), function(K) lapply(c(1, 2), function(de) mk(K = K, delta = de, size_effect = TRUE))),
  # misspecified fixed effects
  lapply(c(100L, 300L), function(K) lapply(1, function(.) mk(K = K, delta = 2, drop_x = TRUE)))
), function(grp) do.call(rbind, lapply(unlist(grp, recursive = FALSE), as.data.frame))))
cells$cid <- seq_len(nrow(cells)); NTOT <- nrow(cells)
Rfor <- function(K) if (QUICK) 3L else 100L

sel <- Sys.getenv("SIM_CELLS", "")
if (QUICK) cells <- cells[c(1, 9, 13, 19, 23), ] else if (nzchar(sel)) {
  rng <- as.integer(strsplit(sel, "-")[[1]]); cells <- cells[cells$cid >= rng[1] & cells$cid <= rng[2], ] }
repsel <- Sys.getenv("SIM_REPS_RANGE", "")
cat(sprintf("run_hard: %d/%d cells\n", nrow(cells), NTOT))
dir.create("sim/results_hard", showWarnings = FALSE, recursive = TRUE)
writeLines(c(capture.output(sessionInfo()), paste("RNGkind:", paste(RNGkind(), collapse = ","))),
           "sim/results_hard/_sessionInfo.txt")

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]; t0 <- Sys.time(); Rcell <- Rfor(cell$K)
  rseq <- if (nzchar(repsel)) { v <- as.integer(strsplit(repsel, "-")[[1]]); seq(v[1], min(v[2], Rcell)) } else seq_len(Rcell)
  reps <- lapply(rseq, function(r) tryCatch(cbind(cid = cell$cid, run_rep_hard(cell, r)), error = function(e) NULL))
  ok <- Filter(Negate(is.null), reps); res <- bind_rows(ok)
  meta <- cell[rep(1, nrow(res)), c("K","nk","p","rho_x","regime","pi","delta","icc","err","nu","null_sd","size_effect","drop_x")]
  out <- cbind(meta, res); out$reps_ok <- nrow(res); out$reps_target <- Rcell   # (#5) record completion
  fname <- if (nzchar(repsel)) sprintf("sim/results_hard/cell_%04d%03d.csv", cell$cid, min(rseq)) else
           sprintf("sim/results_hard/cell_%04d.csv", cell$cid)
  write.csv(out, fname, row.names = FALSE)
  cat(sprintf("  [%d/%d cid=%d K=%d %s d=%s ns=%s sz=%s ms=%s] %d/%d reps %.1fs\n", i, nrow(cells), cell$cid,
              cell$K, cell$regime, cell$delta, cell$null_sd, cell$size_effect, cell$drop_x,
              nrow(res), length(rseq), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("done -> sim/results_hard/cell_*.csv\n")
