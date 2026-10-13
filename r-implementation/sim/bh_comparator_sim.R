# bh_comparator_sim.R -- per-cluster z-test comparators against SCS and SCS+CPSS in the sparse regime.
# Methods per replicate (same data):
#   SCS        : scs_fit at the within-cluster-CV budget (cv_within of run_sims.R, grid incl. 0)
#   SCS+CPSS   : stability(), B = 100 complementary half-samples, pi_thr = 0.8, at the same budget
#   BH(.05/.10): Benjamini-Hochberg on two-sided p of z_k, the standardized residual total of cluster k
#                under the fixed-effects-only OLS fit (the standard per-provider O-minus-E test)
#   Bonf(.05)  : Bonferroni on the same p-values
#   BLUP-top   : the lambda largest |BLUP| (lme4 REML), budget matched to SCS
# Outputs one row per (replicate, method): TP, FP (= V), selected count, plus the CPSS reference.
# Env: SIM_CELL (1..36), SIM_REPS (default 100), SIM_OUT (default sim/results_bh). Seeds: per-rep seed =
# 7e5 * cell + rep, recorded in the output.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R"); library(lme4)
  HAVE_GLMNET <- FALSE; HAVE_RPART <- FALSE
}))
RNGkind("Mersenne-Twister")
ex <- parse("sim/run_sims.R")                       # reuse the paper's lmm_fit and cv_within verbatim
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
                  as.character(e[[2]]) %in% c("lmm_fit", "cv_within")) eval(e)
stopifnot(exists("lmm_fit"), exists("cv_within"))

# SIM_DESIGN = base (36-cell grid below) | clustercv | highcorr (8-cell grids, unequal sizes):
#   clustercv: X gains a cluster-level covariate with correlation 0.7 with the deviation indicator
#   highcorr : rho_x = 0.95 (covariates vary mainly between clusters)
# These are the designs in which the marginal |z| ranking misses the brute-force optimum (Table S1).
design <- Sys.getenv("SIM_DESIGN", "base")
grid_cells <- if (design == "base")
  expand.grid(K = c(100L, 300L), pi = c(0.05, 0.1, 0.2), delta = c(0.5, 1, 2), unequal = c(FALSE, TRUE)) else
  expand.grid(K = c(100L, 300L), pi = c(0.1, 0.2), delta = c(1, 2), unequal = TRUE)
rho_x <- if (design == "highcorr") 0.95 else 0.5
add_cluster_cov <- function(d, rho = 0.7) {
  z <- as.numeric(d$gamma != 0); w <- rho * scale(z) + sqrt(1 - rho^2) * rnorm(length(z))
  d$X <- cbind(d$X, w[as.integer(d$cluster)]); d }
cid  <- as.integer(Sys.getenv("SIM_CELL", "1")); cell <- grid_cells[cid, ]
R    <- as.integer(Sys.getenv("SIM_REPS", "100"))
Bst  <- as.integer(Sys.getenv("SIM_B", "100"))
outd <- Sys.getenv("SIM_OUT", "sim/results_bh")
cat(sprintf("cell %d: K=%d pi=%.2f delta=%.2f unequal=%s R=%d B=%d\n", cid, cell$K, cell$pi, cell$delta, cell$unequal, R, Bst))

zscores <- function(X, y, cl) {
  Xf <- cbind(1, X); f <- lm.fit(Xf, y); r <- f$residuals
  XtXi <- chol2inv(qr.R(qr(Xf))); G <- rowsum(Xf, cl); nk <- as.numeric(table(cl))
  den <- nk - rowSums((G %*% XtXi) * G); s2 <- sum(r^2) / (length(y) - ncol(Xf))
  as.numeric(rowsum(r, cl)) / sqrt(s2 * den)
}
score <- function(sel, active) c(TP = length(intersect(sel, active)), FP = length(setdiff(sel, active)), n_sel = length(sel))

rows <- vector("list", R)
for (r in seq_len(R)) {
  seed <- as.integer(7e5 * cid + r + switch(design, base = 0, clustercv = 5e7, highcorr = 1e8))
  d <- gen_scs(K = cell$K, nk = 20L, p = 5L, regime = "sparse", pi_active = cell$pi, delta = cell$delta,
               rho_x = rho_x, unequal_nk = cell$unequal, seed = seed)
  if (design == "clustercv") d <- add_cluster_cov(d)
  K <- d$K; act <- d$active
  grid <- unique(pmax(0L, as.integer(round(K * c(0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.5)))))
  lam <- cv_within(d, grid, seed = r)$lambda
  out <- list()
  if (lam > 0) {
    fit <- scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r)
    out$SCS <- score(match(fit$selected, levels(d$cluster)), act)
    st <- stability(fit, B = Bst, pi_thr = 0.8, seed = r, pairing = "complementary")
    out$SCS_CPSS <- c(score(match(st$flagged, levels(d$cluster)), act), ref = st$expected_false_flags_bound)
  } else {
    out$SCS <- score(integer(0), act); out$SCS_CPSS <- c(score(integer(0), act), ref = 0)
  }
  pz <- 2 * pnorm(-abs(zscores(d$X, d$y, d$cluster)))
  out$BH05 <- score(which(p.adjust(pz, "BH") <= 0.05), act)
  out$BH10 <- score(which(p.adjust(pz, "BH") <= 0.10), act)
  out$Bonf05 <- score(which(p.adjust(pz, "bonferroni") <= 0.05), act)
  lm0 <- tryCatch(lmm_fit(d$X, d$y, d$cluster), error = function(e) NULL)
  if (!is.null(lm0) && lam > 0) {
    b <- lm0$blup[as.character(seq_len(K))]; b[is.na(b)] <- 0
    out$BLUP_top <- score(order(abs(b), decreasing = TRUE)[seq_len(lam)], act)
  }
  rows[[r]] <- do.call(rbind, lapply(names(out), function(m) data.frame(
    design = design, cell = cid, K = cell$K, pi = cell$pi, delta = cell$delta, unequal = cell$unequal, rep = r, seed = seed,
    lambda = lam, n_active = length(act), method = m, TP = out[[m]][["TP"]], FP = out[[m]][["FP"]],
    n_sel = out[[m]][["n_sel"]], ref = if ("ref" %in% names(out[[m]])) out[[m]][["ref"]] else NA_real_)))
  if (r %% 10 == 0) cat(sprintf("  rep %d/%d done\n", r, R))
}
res <- do.call(rbind, rows)
dir.create(outd, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(outd, sprintf("cell_%02d.csv", cid)), row.names = FALSE)
cat("saved", file.path(outd, sprintf("cell_%02d.csv", cid)), "\n")
