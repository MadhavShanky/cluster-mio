# run_sims.R -- FULL simulation driver.
# Produces: the DGP (sim/dgp.R), the summary table, and the sparsity x K x signal phase map.
# Design notes: right-sized CORE grid + 1-D sweeps (not a full factorial); lambda=0
# allowed (no forced false flag under the null); WITHIN-cluster MSE restored; group-lasso baseline that
# tunes its OWN threshold; regime Sx (prediction-only Z) for the new-cluster gate; robustness sweeps
# (unequal n_k, heavy-tailed errors); failed-fit counts + sessionInfo/RNGkind recorded.
#
# SHARDING: env SIM_CELLS="start-end" runs a cid slice -> sim/results/cell_<id>.csv. SIM_QUICK=TRUE = smoke.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R"); library(lme4); library(dplyr); library(Matrix)
  HAVE_GLMNET <- requireNamespace("glmnet", quietly = TRUE)
  HAVE_RPART  <- requireNamespace("rpart",  quietly = TRUE)
}))
QUICK <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
RNGkind("Mersenne-Twister")

## ---- baselines ----------------------------------------------------------------------------------
lmm_fit <- function(X, y, cluster) {
  df <- data.frame(y = y, X, cluster = cluster)
  fm <- suppressWarnings(suppressMessages(lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  vc <- as.data.frame(VarCorr(fm))
  list(blup = setNames(ranef(fm)$cluster[, 1], rownames(ranef(fm)$cluster)), fixef = lme4::fixef(fm),
       tau2 = vc$vcov[vc$grp == "cluster"], sigma2 = vc$vcov[vc$grp == "Residual"])
}
ols_beta <- function(X, y) { cf <- stats::lm.fit(cbind(1, X), y)$coefficients; cf[is.na(cf)] <- 0; cf[-1] }

# group-lasso baseline: lasso on cluster dummies (X + intercept UNpenalized), CV-tuned -> own sparsity.
glasso_fit <- function(X, y, cluster, K) {
  if (!HAVE_GLMNET) return(NULL)
  D <- Matrix::sparse.model.matrix(~ cluster - 1)
  M <- cbind(Matrix(cbind(1, X), sparse = TRUE), D); p1 <- ncol(X) + 1
  pf <- c(rep(0, p1), rep(1, K))
  cvf <- tryCatch(glmnet::cv.glmnet(M, y, penalty.factor = pf, standardize = FALSE, intercept = FALSE, nfolds = 5),
                  error = function(e) NULL)
  if (is.null(cvf)) return(NULL)
  co <- as.numeric(coef(cvf, s = "lambda.min"))[-1]
  g <- co[(p1 + 1):(p1 + K)]; names(g) <- levels(cluster)
  list(gamma = g, beta = co[2:p1], selected = names(g)[g != 0])
}

# within-cluster K-fold CV: choose lambda (0 ALLOWED); also return L0 & LMM within-cluster test-MSE.
cv_within <- function(d, grid, nfold = 5L, seed = 1L) {
  set.seed(seed); n <- length(d$y); fold <- sample(rep_len(seq_len(nfold), n)); Xf <- cbind(1, d$X)
  cvm <- numeric(length(grid)); cnt <- numeric(length(grid)); lmm_err <- numeric(0)
  for (f in seq_len(nfold)) {
    tr <- fold != f; te <- !tr; cl_tr <- droplevels(d$cluster[tr])
    for (gi in seq_along(grid)) {
      lam <- min(grid[gi], nlevels(cl_tr) - 1L)
      fit <- tryCatch(scs_fit(d$X[tr, , drop = FALSE], d$y[tr], cl_tr, lambda = lam, n_restart = 2L, seed = 1L),
                      error = function(e) NULL)
      if (is.null(fit)) next
      g <- fit$gamma[as.character(d$cluster[te])]; g[is.na(g)] <- 0
      cvm[gi] <- cvm[gi] + mean((d$y[te] - (as.numeric(Xf[te, , drop = FALSE] %*% fit$fixef) + g))^2)
      cnt[gi] <- cnt[gi] + 1
    }
    lf <- tryCatch(lmm_fit(d$X[tr, , drop = FALSE], d$y[tr], cl_tr), error = function(e) NULL)
    if (!is.null(lf)) {
      gb <- lf$blup[as.character(d$cluster[te])]; gb[is.na(gb)] <- 0
      lmm_err <- c(lmm_err, mean((d$y[te] - (as.numeric(Xf[te, , drop = FALSE] %*% lf$fixef) + gb))^2))
    }
  }
  cvm <- cvm / pmax(cnt, 1)
  list(lambda = grid[which.min(cvm)], mse_l0 = min(cvm), mse_lmm = if (length(lmm_err)) mean(lmm_err) else NA)
}

# new-cluster prediction via a Z->gamma_hat tree (regime sparse_x); else X-mean CART (predict.scs default).
ztree_pred <- function(gtrain, Ztr, Zte, te_rows_cluster, teK) {
  if (!HAVE_RPART) return(rep(0, length(te_rows_cluster)))
  colnames(Ztr) <- paste0("z", seq_len(ncol(Ztr))); colnames(Zte) <- colnames(Ztr)
  tr <- rpart::rpart(g ~ ., data = data.frame(g = as.numeric(gtrain), Ztr), method = "anova",
                     control = rpart::rpart.control(maxdepth = 4))
  gK <- setNames(as.numeric(predict(tr, as.data.frame(Zte))), as.character(teK))
  gK[as.character(te_rows_cluster)]
}
cb_mse <- function(d, lam, nfold = 5L, seed = 1L) {
  set.seed(seed); K <- d$K; foldk <- sample(rep_len(seq_len(nfold), K)); useZ <- d$regime == "sparse_x" && !is.null(d$Z)
  el0 <- elmm <- numeric(0)
  for (f in seq_len(nfold)) {
    teK <- which(foldk == f); idx <- d$cluster %in% as.character(teK); trK <- setdiff(seq_len(K), teK)
    if (!any(idx) || all(idx)) next
    Xtr <- d$X[!idx, , drop = FALSE]; ytr <- d$y[!idx]; cltr <- droplevels(d$cluster[!idx])
    Xte <- d$X[idx, , drop = FALSE]; yte <- d$y[idx]
    fit <- tryCatch(scs_fit(Xtr, ytr, cltr, lambda = min(lam, nlevels(cltr) - 1L), n_restart = 2L), error = function(e) NULL)
    if (!is.null(fit)) {
      pr <- if (useZ) as.numeric(cbind(1, Xte) %*% fit$fixef) +
              ztree_pred(fit$gamma[as.character(trK)], d$Z[trK, , drop = FALSE], d$Z[teK, , drop = FALSE], d$cluster[idx], teK)
            else predict(fit, newdata = Xte)
      el0 <- c(el0, mean((yte - pr)^2))
    }
    lm <- tryCatch(lmm_fit(Xtr, ytr, cltr), error = function(e) NULL)
    if (!is.null(lm)) {
      pr <- if (useZ) as.numeric(cbind(1, Xte) %*% lm$fixef) +
              ztree_pred(lm$blup[as.character(trK)], d$Z[trK, , drop = FALSE], d$Z[teK, , drop = FALSE], d$cluster[idx], teK)
            else as.numeric(cbind(1, Xte) %*% lm$fixef)               # new cluster -> gamma = 0
      elmm <- c(elmm, mean((yte - pr)^2))
    }
  }
  c(l0 = if (length(el0)) mean(el0) else NA, lmm = if (length(elmm)) mean(elmm) else NA)
}

## ---- one rep ------------------------------------------------------------------------------------
run_rep <- function(cell, r) {
  d <- gen_scs(cell$K, cell$nk, cell$p, regime = cell$regime, pi_active = cell$pi, delta = cell$delta,
               icc = cell$icc, rho_x = cell$rho_x, err = cell$err, unequal_nk = cell$unequal,
               seed = as.integer((1e6 * cell$cid) %% 2e9) + r)
  K <- d$K; sparse <- cell$regime != "gaussian"
  grid <- unique(pmax(0L, as.integer(round(K * c(0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.5)))))   # includes 0
  o <- list(rep = r)

  t0 <- Sys.time(); cv <- cv_within(d, grid, seed = r); lam <- cv$lambda
  fitL <- scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r)
  o$runtime_l0 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  gL <- fitL$gamma[as.character(seq_len(K))]; gL[is.na(gL)] <- 0
  lm0 <- tryCatch(lmm_fit(d$X, d$y, d$cluster), error = function(e) NULL)
  o$fail_lmm <- is.null(lm0)
  blup <- if (!is.null(lm0)) { v <- lm0$blup[as.character(seq_len(K))]; v[is.na(v)] <- 0; v } else rep(0, K)
  gl <- glasso_fit(d$X, d$y, d$cluster, K); o$fail_glasso <- is.null(gl)
  o$lambda <- lam; o$n_obs <- length(d$y); o$wmse_l0 <- cv$mse_l0; o$wmse_lmm <- cv$mse_lmm

  if (sparse) {
    o <- c(o, with(sel_scores(match(fitL$selected, levels(d$cluster)), d$active, K),
                   list(f1_l0 = f1, fdr_l0 = fdr, tpr_l0 = recall)))
    selM <- if (lam > 0) order(abs(blup), decreasing = TRUE)[seq_len(lam)] else integer(0)  # budget-matched ranking
    o <- c(o, with(sel_scores(selM, d$active, K), list(f1_lmm = f1, fdr_lmm = fdr, tpr_lmm = recall)))
    if (!is.null(gl)) o <- c(o, with(sel_scores(match(gl$selected, levels(d$cluster)), d$active, K),
                                     list(f1_glasso = f1, fdr_glasso = fdr, tpr_glasso = recall)))
    st <- tryCatch(stability(fitL, B = if (QUICK) 20L else 50L, pi_thr = 0.8, seed = r), error = function(e) NULL)
    if (!is.null(st)) o <- c(o, with(sel_scores(match(st$flagged, levels(d$cluster)), d$active, K),
                                     list(f1_stab = f1, fdr_stab = fdr)), list(mb_bound = st$expected_false_flags_bound))
  } else {
    o$icc_true <- cell$icc; o$icc_lmm <- if (!is.null(lm0)) lm0$tau2 / (lm0$tau2 + lm0$sigma2) else NA
  }
  cb <- cb_mse(d, lam, seed = r)
  o$gerr_l0 <- sqrt(sum((gL - d$gamma)^2)); o$gerr_lmm <- sqrt(sum((blup - d$gamma)^2))
  if (!is.null(gl)) o$gerr_glasso <- sqrt(sum((gl$gamma[as.character(seq_len(K))] - d$gamma)^2))
  o$berr_l0 <- sqrt(sum((fitL$fixef[-1] - d$beta)^2)); o$berr_ols <- sqrt(sum((ols_beta(d$X, d$y) - d$beta)^2))
  if (!is.null(lm0)) o$berr_lmm <- sqrt(sum((lm0$fixef[-1] - d$beta)^2))
  o$cbmse_l0 <- cb["l0"]; o$cbmse_lmm <- cb["lmm"]
  as.data.frame(o, stringsAsFactors = FALSE)
}

## ---- CORE grid + targeted 1-D sweeps -------------------------------------------------------------
ref <- list(nk = 30L, p = 5L, rho_x = 0.5)
mkcell <- function(...) { a <- list(...); modifyList(list(K = 100L, nk = ref$nk, p = ref$p, rho_x = ref$rho_x,
  regime = "sparse", pi = 0.1, delta = 2, icc = NA, err = "normal", unequal = FALSE), a) }
core <- do.call(rbind, lapply(c(10, 30, 100, 300, 1000), function(K)
  do.call(rbind, lapply(c(.05, .1, .2, .4, .6, .8), function(pi)
    do.call(rbind, lapply(c(1, 2, 3), function(de) as.data.frame(mkcell(K = K, pi = pi, delta = de))))))))
sweeps <- do.call(rbind, list(                                # off reference (K=100, pi=.1, delta=2)
  as.data.frame(mkcell(nk = 20L)), as.data.frame(mkcell(nk = 50L)), as.data.frame(mkcell(p = 25L)),
  as.data.frame(mkcell(rho_x = 0)), as.data.frame(mkcell(unequal = TRUE)), as.data.frame(mkcell(err = "t"))))
gaus <- do.call(rbind, lapply(c(30, 100, 300), function(K)
  do.call(rbind, lapply(c(.1, .3, .5, .7, .9), function(ic)
    as.data.frame(mkcell(K = K, regime = "gaussian", pi = NA, delta = NA, icc = ic))))))
sx <- do.call(rbind, lapply(c(30, 100, 300), function(K)
  do.call(rbind, lapply(c(.1, .2, .4), function(pi) as.data.frame(mkcell(K = K, regime = "sparse_x", pi = pi, delta = 2))))))
cells <- rbind(core, sweeps, gaus, sx); cells$cid <- seq_len(nrow(cells)); NTOT <- nrow(cells)
Rfor <- function(K) if (QUICK) 3L else if (K >= 1000L) 100L else 200L   # K=1000 reps halved (cell-means stable)

sel <- Sys.getenv("SIM_CELLS", "")
if (QUICK) cells <- cells[c(8, which(cells$regime == "sparse_x")[1], which(cells$regime == "gaussian")[1],
                            which(cells$unequal)[1]), ] else if (nzchar(sel)) {
  rng <- as.integer(strsplit(sel, "-")[[1]]); cells <- cells[cells$cid >= rng[1] & cells$cid <= rng[2], ] }
# REP-SHARDING: SIM_REPS_RANGE="start-end" runs only rep indices start..end of each cid (seeds keyed on r,
# so disjoint shards = disjoint reproducible draws). Output goes to cell_<cid><start>.csv (digit-only ->
# still matched by the aggregator's ^cell_[0-9]+\.csv$; scenario meta cols pool reps across shards exactly).
repsel <- Sys.getenv("SIM_REPS_RANGE", "")
cat(sprintf("run_sims: %d/%d cells (R=200, K1000->100; core %d + sweeps %d + gaus %d + sx %d) glmnet=%s\n",
            nrow(cells), NTOT, nrow(core), nrow(sweeps), nrow(gaus), nrow(sx), HAVE_GLMNET))
dir.create("sim/results", showWarnings = FALSE, recursive = TRUE)
writeLines(c(capture.output(sessionInfo()), paste("RNGkind:", paste(RNGkind(), collapse = ","))),
           "sim/results/_sessionInfo.txt")

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]; t0 <- Sys.time(); Rcell <- Rfor(cell$K)
  rseq <- if (nzchar(repsel)) { v <- as.integer(strsplit(repsel, "-")[[1]]); seq(v[1], min(v[2], Rcell)) } else seq_len(Rcell)
  reps <- lapply(rseq, function(r) tryCatch(cbind(cid = cell$cid, run_rep(cell, r)), error = function(e) NULL))
  ok <- Filter(Negate(is.null), reps); res <- bind_rows(ok)
  meta <- cell[rep(1, nrow(res)), c("K","nk","p","rho_x","regime","pi","delta","icc","err","unequal")]
  out <- cbind(meta, res); out$reps_ok <- nrow(res); out$reps_target <- Rcell
  fname <- if (nzchar(repsel)) sprintf("sim/results/cell_%04d%03d.csv", cell$cid, min(rseq)) else
           sprintf("sim/results/cell_%04d.csv", cell$cid)
  write.csv(out, fname, row.names = FALSE)
  cat(sprintf("  [%d/%d cid=%d K=%d %s reps %d-%d] %d/%d reps, %.1fs\n", i, nrow(cells), cell$cid, cell$K,
              cell$regime, min(rseq), max(rseq), nrow(res), length(rseq), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("done -> sim/results/cell_*.csv (+ _sessionInfo.txt)\n")
