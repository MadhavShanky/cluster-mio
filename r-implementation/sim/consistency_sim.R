# consistency_sim.R -- validates the L0 support-consistency / sure-screening proposition
# (reviews/strengthen_synthesis.md sec.2A; Codex proof in reviews/proof_codex_out.txt).
# Estimand: signed mean-shift support S* = {k : gamma_k != 0}. For the GLOBAL-ish L0 estimator we measure
#   exact recovery  P(hatS = S*)   at lambda = s = |S*|,
#   sure screening  P(S* subset hatS) at lambda = s and lambda = s + slack,
# as functions of K, signal delta, cluster size, with the ridge eps and restart count swept WITHIN each rep
# (same data) so the ridge/restart Jaccard measures the near-tie fragility (C3b) at fixed signal.
#
# DGP: gen_scs() regime="sparse" (gamma_k = +/- delta on round(pi*K) clusters, gamma _|_ X) -- encodes the
#   exact mean-shift-support assumption the theorem is about.
# Sharding: SIM_CELLS="a-b" selects cell ids; SIM_REPS_RANGE="a-b" a rep slice (seeds keyed on r). Output
#   sim/results_consistency/cell_<cid>[<repstart>].csv .  Aggregate: rbind all cell_*.csv.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R")
}))
RNGkind("Mersenne-Twister")
QUICK <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))

EPS   <- c(0, 1e-6, 1e-4, 1e-2)     # ridge sweep; eps=1e-4 is the paper default (ridge reference)
EPSREF <- 1e-4
NRES  <- c(1L, 4L)                  # local-search restarts; 4 is the paper default (restart reference)
NRREF <- 4L
PI_ACTIVE <- 0.1; P <- 5L; NK <- 30L

## ---- cell grid: (K, imbalance, delta) ; ridge/restart/lambda swept inside each rep ----
grid <- expand.grid(K = c(100L, 300L, 1000L), unequal = c(FALSE, TRUE), delta = c(1, 2, 3),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cid <- seq_len(nrow(grid))
Rfor <- function(K) if (QUICK) 3L else if (K >= 1000L) 100L else 200L

jacc <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

run_rep <- function(cell, r) {
  d <- gen_scs(cell$K, NK, P, regime = "sparse", pi_active = PI_ACTIVE, delta = cell$delta,
               rho_x = 0.5, unequal_nk = cell$unequal,
               seed = as.integer((1e6 * cell$cid) %% 2e9) + r)
  K <- d$K; active <- as.integer(d$active); s <- length(active)
  slack <- max(2L, as.integer(round(0.5 * s)))
  lambdas <- c(exact = s, screen = s + slack)
  mmin <- min(d$nk_vec); snr <- cell$delta^2 * mmin / log(K)   # theory transition coordinate

  # fit under every (eps, n_restart, lambda) on the SAME data; store selected index sets
  rows <- list()
  selstore <- list()                      # key "eps|nr|lamname" -> integer selected indices
  for (eps in EPS) for (nr in NRES) for (li in seq_along(lambdas)) {
    lamname <- names(lambdas)[li]; lam <- min(lambdas[[li]], K - 1L)
    fit <- tryCatch(scs_fit(d$X, d$y, d$cluster, lambda = lam, eps = eps,
                            n_restart = nr, seed = r), error = function(e) NULL)
    sel <- if (is.null(fit)) integer(0) else match(fit$selected, levels(d$cluster))
    selstore[[paste(eps, nr, lamname, sep = "|")]] <- sel
    sc <- sel_scores(sel, active, K)
    rows[[length(rows) + 1L]] <- data.frame(
      cid = cell$cid, K = K, unequal = cell$unequal, delta = cell$delta, s = s, mmin = mmin, snr = snr,
      rep = r, eps = eps, n_restart = nr, lambda_kind = lamname, lambda = lam, n_sel = length(sel),
      exact = as.integer(setequal(sel, active)),
      screen = as.integer(length(active) > 0 && all(active %in% sel)),
      f1 = sc$f1, fdr = sc$fdr, tpr = sc$recall, n_false = sc$fp,
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  # fragility Jaccards (same data): selected set vs the ridge-reference and restart-reference
  out$jacc_vs_epsref <- mapply(function(eps, nr, lk)
    jacc(selstore[[paste(eps, nr, lk, sep = "|")]], selstore[[paste(EPSREF, nr, lk, sep = "|")]]),
    out$eps, out$n_restart, out$lambda_kind)
  out$jacc_vs_nrref <- mapply(function(eps, nr, lk)
    jacc(selstore[[paste(eps, nr, lk, sep = "|")]], selstore[[paste(eps, NRREF, lk, sep = "|")]]),
    out$eps, out$n_restart, out$lambda_kind)
  out
}

## ---- driver ----
sel <- Sys.getenv("SIM_CELLS", "")
cells <- if (nzchar(sel)) { rng <- as.integer(strsplit(sel, "-")[[1]]); grid[grid$cid >= rng[1] & grid$cid <= rng[2], ] } else grid
repsel <- Sys.getenv("SIM_REPS_RANGE", "")
dir.create("sim/results_consistency", showWarnings = FALSE, recursive = TRUE)
cat(sprintf("consistency_sim: %d cells (eps x restart x lambda swept within rep)\n", nrow(cells)))

for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]; t0 <- Sys.time(); Rc <- Rfor(cell$K)
  rseq <- if (nzchar(repsel)) { v <- as.integer(strsplit(repsel, "-")[[1]]); seq(v[1], min(v[2], Rc)) } else seq_len(Rc)
  res <- do.call(rbind, lapply(rseq, function(r) tryCatch(run_rep(cell, r), error = function(e) NULL)))
  fname <- if (nzchar(repsel)) sprintf("sim/results_consistency/cell_%04d%03d.csv", cell$cid, min(rseq)) else
           sprintf("sim/results_consistency/cell_%04d.csv", cell$cid)
  write.csv(res, fname, row.names = FALSE)
  cat(sprintf("  [cid=%d K=%d unequal=%s delta=%g] reps %d-%d, %d rows, %.1fs\n",
              cell$cid, cell$K, cell$unequal, cell$delta, min(rseq), max(rseq), nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("done -> sim/results_consistency/cell_*.csv\n")
