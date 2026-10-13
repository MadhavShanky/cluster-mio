# subsample_scheme_sim.R -- observation-level vs cluster-stratified subsampling in stability selection.
# Tests the enabler (reviews/strengthen_synthesis.md sec.3): cluster-stratified subsampling keeps every
# cluster eligible in every subsample, so (a) small active clusters are no longer capped below pi_thr by
# 1-(1-frac)^{n_k}, (b) the MB/CPSS denominator K is legitimate, (c) the watch-list is less ridge-fragile.
#
# Per (rep, scheme, eps): run the EXACT stability() (B, frac=0.5, pi_thr) under scheme in {obs,cluster_strat}
# and eps in {paper default, 0}. Record realized E[V] (false flags = flagged inactive), TPR overall and by
# cluster-size tercile, q, MB reference, watch-list size, and Jaccard of the flagged set vs the obs-scheme
# run at the same eps (how much stratification changes the list) and vs eps=1e-4 same scheme (ridge fragility).
# DGP: gen_scs regime="sparse"; unequal_nk exposes the small-cluster cap. K limited to {100,300} (stability
# at K=1000 is the known >12h cost).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R")
}))
RNGkind("Mersenne-Twister")
QUICK  <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
B      <- if (QUICK) 20L else as.integer(Sys.getenv("SIM_B", "100"))
FRAC   <- 0.5; PITHR <- 0.8; PI_ACTIVE <- 0.1; P <- 5L; NK <- 30L
EPS    <- c(1e-4, 0); EPSREF <- 1e-4
SCHEMES <- c("obs", "cluster_strat")

grid <- expand.grid(K = c(100L, 300L), unequal = c(FALSE, TRUE), delta = c(1, 2),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cid <- seq_len(nrow(grid))
Rfor <- function(K) if (QUICK) 3L else if (K >= 300L) 60L else 100L
jacc <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

run_rep <- function(cell, r) {
  d <- gen_scs(cell$K, NK, P, regime = "sparse", pi_active = PI_ACTIVE, delta = cell$delta,
               rho_x = 0.5, unequal_nk = cell$unequal,
               seed = as.integer((1e6 * cell$cid) %% 2e9) + r)
  K <- d$K; lev <- levels(d$cluster); active <- as.integer(d$active); s <- length(active)
  nk <- d$nk_vec
  rk <- rank(nk, ties.method = "first")     # rank-terciles: distinct breaks even when all n_k equal (balanced)
  terc <- cut(rk, quantile(rk, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("small","mid","large"))
  lam <- min(s, K - 1L)
  fit <- tryCatch(scs_fit(d$X, d$y, d$cluster, lambda = lam, eps = EPSREF, n_restart = 4L, seed = r),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  flagstore <- list()                         # key "scheme|eps" -> integer flagged indices
  rows <- list()
  for (sch in SCHEMES) for (eps in EPS) {
    # refit selection ridge at this eps so the stability base selector matches (fit object carries mu)
    f0 <- if (eps == EPSREF) fit else tryCatch(scs_fit(d$X, d$y, d$cluster, lambda = lam, eps = eps,
                                                       n_restart = 4L, seed = r), error = function(e) NULL)
    if (is.null(f0)) next
    st <- tryCatch(stability(f0, B = B, frac = FRAC, pi_thr = PITHR, scheme = sch, n_restart = 2L, seed = r),
                   error = function(e) NULL)
    if (is.null(st)) next
    fl <- match(st$flagged, lev)
    flagstore[[paste(sch, eps, sep = "|")]] <- fl
    tp <- intersect(fl, active); fp <- setdiff(fl, active)
    tpr_t <- sapply(c("small","mid","large"), function(tt) {
      at <- active[terc[active] == tt]; if (length(at) == 0) NA_real_ else mean(at %in% fl) })
    rows[[length(rows) + 1L]] <- data.frame(
      cid = cell$cid, K = K, unequal = cell$unequal, delta = cell$delta, s = s, rep = r,
      scheme = sch, eps = eps, B_used = st$B_used, q = st$q, mb_ref = st$expected_false_flags_bound,
      n_watch = length(fl), V_false = length(fp), tpr = if (s > 0) length(tp) / s else NA_real_,
      tpr_small = tpr_t[["small"]], tpr_mid = tpr_t[["mid"]], tpr_large = tpr_t[["large"]],
      n_active_small = sum(terc[active] == "small"), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$jacc_vs_obs <- mapply(function(sch, eps)
    jacc(flagstore[[paste(sch, eps, sep = "|")]], flagstore[[paste("obs", eps, sep = "|")]]),
    out$scheme, out$eps)
  out$jacc_vs_epsref <- mapply(function(sch, eps)
    jacc(flagstore[[paste(sch, eps, sep = "|")]], flagstore[[paste(sch, EPSREF, sep = "|")]]),
    out$scheme, out$eps)
  out
}

sel <- Sys.getenv("SIM_CELLS", "")
cells <- if (nzchar(sel)) { rng <- as.integer(strsplit(sel, "-")[[1]]); grid[grid$cid >= rng[1] & grid$cid <= rng[2], ] } else grid
repsel <- Sys.getenv("SIM_REPS_RANGE", "")
dir.create("sim/results_scheme", showWarnings = FALSE, recursive = TRUE)
cat(sprintf("subsample_scheme_sim: %d cells, B=%d\n", nrow(cells), B))
for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]; t0 <- Sys.time(); Rc <- Rfor(cell$K)
  rseq <- if (nzchar(repsel)) { v <- as.integer(strsplit(repsel, "-")[[1]]); seq(v[1], min(v[2], Rc)) } else seq_len(Rc)
  res <- do.call(rbind, lapply(rseq, function(r) tryCatch(run_rep(cell, r), error = function(e) NULL)))
  fname <- if (nzchar(repsel)) sprintf("sim/results_scheme/cell_%04d%03d.csv", cell$cid, min(rseq)) else
           sprintf("sim/results_scheme/cell_%04d.csv", cell$cid)
  write.csv(res, fname, row.names = FALSE)
  cat(sprintf("  [cid=%d K=%d unequal=%s delta=%g] reps %d-%d, %d rows, %.1fs\n",
              cell$cid, cell$K, cell$unequal, cell$delta, min(rseq), max(rseq), nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("done -> sim/results_scheme/cell_*.csv\n")
