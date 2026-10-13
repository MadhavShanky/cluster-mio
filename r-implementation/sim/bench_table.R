# bench_table.R -- single clean benchmark table requested by the editor (ref2 editor item):
# separate SINGLE FIT / CV PATH / STABILITY LAYER / COMPETITORS under IDENTICAL conditions, across K.
# Each method is timed on the SAME generated instance; SCS optimality is checked against brute force where
# the support is enumerable, and against L0Learn on the common unpenalized objective R(S).
#
# SHARDING (for cluster parallelism): the work is nK x REPS cells. env BENCH_TASK="t" runs ONLY cell t
# (1-based, gi-major: t=1..REPS -> grid row 1, etc.) and writes sim/results_bench/shard_<gi>_<rep>.csv.
# With BENCH_TASK unset the script runs every cell sequentially and then aggregates. Either way the final
# table is built from the shards by sim/bench_table_aggregate.R (sourced here in the full-run case).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R"); library(lme4); library(Matrix)
  HAVE_GLMNET <- requireNamespace("glmnet",  quietly = TRUE)
  HAVE_L0     <- requireNamespace("L0Learn", quietly = TRUE)
}))
RNGkind("Mersenne-Twister")
QUICK <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
B_STAB <- if (QUICK) 10L else 50L
REPS   <- if (QUICK) 2L else 10L
dir.create("sim/results_bench", showWarnings = FALSE, recursive = TRUE)

## ---- competitor wrappers (copied from sim/run_sims.R so the table matches the paper's baselines) --
lmm_fit <- function(X, y, cluster) {
  df <- data.frame(y = y, X, cluster = cluster)
  fm <- suppressWarnings(suppressMessages(lmer(y ~ . - cluster - y + (1 | cluster), data = df, REML = TRUE)))
  invisible(fm)
}
glasso_fit <- function(X, y, cluster, K) {
  if (!HAVE_GLMNET) return(NULL)
  D <- Matrix::sparse.model.matrix(~ cluster - 1)
  M <- cbind(Matrix(cbind(1, X), sparse = TRUE), D); p1 <- ncol(X) + 1
  pf <- c(rep(0, p1), rep(1, K))
  tryCatch(glmnet::cv.glmnet(M, y, penalty.factor = pf, standardize = FALSE, intercept = FALSE, nfolds = 5),
           error = function(e) NULL)
}
# SCS CV path: 5-fold within-cluster CV over the budget grid (SCS-only; this is the cost of CHOOSING lambda)
scs_cv_path <- function(X, y, cluster, grid, nfold = 5L, seed = 1L) {
  set.seed(seed); n <- length(y); fold <- sample(rep_len(seq_len(nfold), n)); Xf <- cbind(1, X)
  cvm <- numeric(length(grid)); cnt <- numeric(length(grid))
  for (f in seq_len(nfold)) {
    tr <- fold != f; te <- !tr; cl_tr <- droplevels(cluster[tr])
    for (gi in seq_along(grid)) {
      lam <- min(grid[gi], nlevels(cl_tr) - 1L)
      fit <- tryCatch(scs_fit(X[tr, , drop = FALSE], y[tr], cl_tr, lambda = lam, n_restart = 2L, seed = 1L),
                      error = function(e) NULL)
      if (is.null(fit)) next
      g <- fit$gamma[as.character(cluster[te])]; g[is.na(g)] <- 0
      cvm[gi] <- cvm[gi] + mean((y[te] - (as.numeric(Xf[te, , drop = FALSE] %*% fit$fixef) + g))^2); cnt[gi] <- cnt[gi] + 1
    }
  }
  grid[which.min(cvm / pmax(cnt, 1))]
}

## ---- common unpenalized objective R(S) and brute force -------------------------------------------
refit_resid <- function(My, MD, S) { if (!length(S)) return(sum(My^2)); sum(stats::lm.fit(MD[, S, drop = FALSE], My)$residuals^2) }
brute_best <- function(My, MD, K, lambda) {
  combs <- utils::combn(K, lambda); best <- Inf
  for (j in seq_len(ncol(combs))) { r <- refit_resid(My, MD, combs[, j]); if (r < best) best <- r }; best
}

bench_one <- function(K, p, nk, lambda, seed) {
  d <- gen_scs(K, nk, p, regime = "sparse", pi_active = 0.05, delta = 2, rho_x = 0.5, seed = seed)
  X <- cbind(1, d$X); D <- model.matrix(~ d$cluster - 1)
  Mfun <- function(V) V - X %*% solve(crossprod(X), crossprod(X, V)); My <- as.numeric(Mfun(d$y)); MD <- Mfun(D)
  grid <- unique(pmax(0L, as.integer(round(K * c(0, 0.02, 0.05, 0.10, 0.20)))))

  tm <- function(expr) as.numeric(system.time(expr)["elapsed"])
  t_single <- tm(fit <- scs_fit(d$X, d$y, d$cluster, lambda = lambda, n_restart = 4L, seed = seed))
  t_cvpath <- tm(scs_cv_path(d$X, d$y, d$cluster, grid, seed = seed))
  t_stab   <- tm(stability(fit, B = B_STAB, pi_thr = 0.8, seed = seed))
  t_lmm    <- tryCatch(tm(lmm_fit(d$X, d$y, d$cluster)), error = function(e) NA_real_)
  t_glasso <- if (HAVE_GLMNET) tryCatch(tm(glasso_fit(d$X, d$y, d$cluster, K)), error = function(e) NA_real_) else NA_real_

  # L0Learn at matched cardinality + objective comparison on R(S)
  S_scs <- match(fit$selected, levels(d$cluster)); obj_scs <- refit_resid(My, MD, S_scs)
  t_l0 <- NA_real_; obj_l0 <- NA_real_
  if (HAVE_L0) {
    t_l0 <- tm(fc <- tryCatch(L0Learn::L0Learn.fit(MD, My, penalty = "L0", maxSuppSize = min(K, lambda + 10L),
                                                   nLambda = 200L, intercept = FALSE), error = function(e) NULL))
    if (!is.null(fc)) { g <- as.matrix(coef(fc)); ss <- colSums(g != 0); cand <- which(ss == lambda)
      if (length(cand)) { best <- Inf; for (cc in cand) { rr <- refit_resid(My, MD, which(g[, cc] != 0)); if (rr < best) best <- rr }; obj_l0 <- best } }
  }
  obj_bf <- if (choose(K, lambda) <= 5e5) brute_best(My, MD, K, lambda) else NA_real_

  data.frame(K = K, p = p, nk = nk, lambda = lambda, seed = seed,
             t_single = t_single, t_cvpath = t_cvpath, t_stab = t_stab,
             t_lmm = t_lmm, t_glasso = t_glasso, t_l0learn = t_l0,
             gap_scs_brute = if (is.finite(obj_bf)) obj_scs - obj_bf else NA_real_,
             gap_l0_brute  = if (is.finite(obj_bf) && is.finite(obj_l0)) obj_l0 - obj_bf else NA_real_,
             gap_scs_l0    = if (is.finite(obj_l0)) obj_scs - obj_l0 else NA_real_)
}

## ---- grid: identical (p, nk), increasing K -------------------------------------------------------
grid <- data.frame(K = c(30L, 100L, 300L, 1000L), p = 5L, nk = 30L)
grid$lambda <- pmax(2L, as.integer(round(0.05 * grid$K)))
if (QUICK) grid <- grid[grid$K <= 100L, ]
grid$gi <- seq_len(nrow(grid))

# the full job list (gi-major: row1 reps 1..REPS, row2 reps 1..REPS, ...) -> BENCH_TASK indexes it
jobs <- do.call(rbind, lapply(grid$gi, function(gi) data.frame(gi = gi, r = seq_len(REPS))))
task <- Sys.getenv("BENCH_TASK", "")
sharded <- nzchar(task)
if (sharded) { t <- as.integer(task); stopifnot(t >= 1L, t <= nrow(jobs)); jobs <- jobs[t, , drop = FALSE] }
cat(sprintf("bench_table: %d cell(s) [%s], B_stab=%d, REPS=%d, glmnet=%s L0Learn=%s\n",
            nrow(jobs), if (sharded) paste0("shard ", task) else "full grid", B_STAB, REPS, HAVE_GLMNET, HAVE_L0))

for (jj in seq_len(nrow(jobs))) {
  gi <- jobs$gi[jj]; r <- jobs$r[jj]; g <- grid[grid$gi == gi, ]
  res <- tryCatch(cbind(gi = gi, rep = r, bench_one(g$K, g$p, g$nk, g$lambda, seed = 1000L * gi + r)),
                  error = function(e) { cat(sprintf("  cell gi=%d r=%d FAILED: %s\n", gi, r, conditionMessage(e))); NULL })
  if (!is.null(res)) {
    write.csv(res, sprintf("sim/results_bench/shard_%02d_%02d.csv", gi, r), row.names = FALSE)
    cat(sprintf("  [gi=%d K=%d r=%d] single=%.3f cvpath=%.3f stab=%.3f lmm=%.3f\n",
                gi, g$K, r, res$t_single, res$t_cvpath, res$t_stab, res$t_lmm))
  }
}
writeLines(c(capture.output(sessionInfo()), paste("glmnet:", HAVE_GLMNET, "L0Learn:", HAVE_L0),
             paste("RNGkind:", paste(RNGkind(), collapse = ","))), "sim/results_bench/_bench_table_sessionInfo.txt")

if (!sharded) { cat("full run done; aggregating ->\n"); source("sim/bench_table_aggregate.R") } else
  cat(sprintf("shard %s done -> run sim/bench_table_aggregate.R after all shards finish\n", task))
