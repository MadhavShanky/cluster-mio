# restart_robustcore_sim.R -- restart-jitter sensitivity and the robust core (C3b fragility, sec.2D/4 of
# reviews/strengthen_synthesis.md). The NY watch-list moved 14<->4 under a tiny ridge change because the
# local-search path is near-tie fragile across stability subsamples. Two questions:
#   (1) does raising n_restart damp that jitter (watch-list stabilizes as restarts grow)?
#   (2) is the "robust core" (selection freq >= 0.9) the stable, defensible report?
# Per (rep, n_restart in {2,5,10,20}): run stability() (obs scheme, paper default) and record watch-list
# size (freq>=pi_thr), robust-core size (freq>=0.9), TPR, false flags V, max null-cluster freq, and Jaccard
# of the watch-list vs the n_restart=20 run on the SAME data (convergence of the list with more restarts).
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R")
}))
RNGkind("Mersenne-Twister")
QUICK <- isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))
B     <- if (QUICK) 20L else as.integer(Sys.getenv("SIM_B", "80"))
FRAC  <- 0.5; PITHR <- 0.8; CORE <- 0.9; PI_ACTIVE <- 0.1; P <- 5L; NK <- 30L
NRES  <- c(2L, 5L, 10L, 20L); NRREF <- 20L

grid <- expand.grid(K = c(100L, 300L), unequal = c(FALSE, TRUE), delta = c(1, 2),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$cid <- seq_len(nrow(grid))
Rfor <- function(K) if (QUICK) 3L else if (K >= 300L) 40L else 80L
jacc <- function(a, b) { u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

run_rep <- function(cell, r) {
  d <- gen_scs(cell$K, NK, P, regime = "sparse", pi_active = PI_ACTIVE, delta = cell$delta,
               rho_x = 0.5, unequal_nk = cell$unequal,
               seed = as.integer((1e6 * cell$cid) %% 2e9) + r)
  K <- d$K; lev <- levels(d$cluster); active <- as.integer(d$active); s <- length(active)
  inactive <- setdiff(seq_len(K), active)
  lam <- min(s, K - 1L)
  fit <- tryCatch(scs_fit(d$X, d$y, d$cluster, lambda = lam, n_restart = 4L, seed = r), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  watchstore <- list(); rows <- list()
  for (nr in NRES) {
    st <- tryCatch(stability(fit, B = B, frac = FRAC, pi_thr = PITHR, scheme = "obs", n_restart = nr, seed = r),
                   error = function(e) NULL)
    if (is.null(st)) next
    freq <- st$freq                                   # named, sorted, over all K
    fi <- match(names(freq), lev)
    watch <- fi[freq >= PITHR]; core <- fi[freq >= CORE]
    watchstore[[as.character(nr)]] <- watch
    null_freq <- freq[fi %in% inactive]
    rows[[length(rows) + 1L]] <- data.frame(
      cid = cell$cid, K = K, unequal = cell$unequal, delta = cell$delta, s = s, rep = r,
      n_restart = nr, B_used = st$B_used, q = st$q, mb_ref = st$expected_false_flags_bound,
      n_watch = length(watch), n_core = length(core),
      tpr = if (s > 0) length(intersect(watch, active)) / s else NA_real_,
      V_false = length(setdiff(watch, active)),
      max_null_freq = if (length(null_freq)) max(null_freq) else 0,
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  out <- do.call(rbind, rows)
  out$jacc_watch_vs_nrref <- sapply(out$n_restart, function(nr)
    jacc(watchstore[[as.character(nr)]], watchstore[[as.character(NRREF)]]))
  out
}

sel <- Sys.getenv("SIM_CELLS", "")
cells <- if (nzchar(sel)) { rng <- as.integer(strsplit(sel, "-")[[1]]); grid[grid$cid >= rng[1] & grid$cid <= rng[2], ] } else grid
repsel <- Sys.getenv("SIM_REPS_RANGE", "")
dir.create("sim/results_restart", showWarnings = FALSE, recursive = TRUE)
cat(sprintf("restart_robustcore_sim: %d cells, B=%d, n_restart in {%s}\n", nrow(cells), B, paste(NRES, collapse=",")))
for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]; t0 <- Sys.time(); Rc <- Rfor(cell$K)
  rseq <- if (nzchar(repsel)) { v <- as.integer(strsplit(repsel, "-")[[1]]); seq(v[1], min(v[2], Rc)) } else seq_len(Rc)
  res <- do.call(rbind, lapply(rseq, function(r) tryCatch(run_rep(cell, r), error = function(e) NULL)))
  fname <- if (nzchar(repsel)) sprintf("sim/results_restart/cell_%04d%03d.csv", cell$cid, min(rseq)) else
           sprintf("sim/results_restart/cell_%04d.csv", cell$cid)
  write.csv(res, fname, row.names = FALSE)
  cat(sprintf("  [cid=%d K=%d unequal=%s delta=%g] reps %d-%d, %d rows, %.1fs\n",
              cell$cid, cell$K, cell$unequal, cell$delta, min(rseq), max(rseq), nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
cat("done -> sim/results_restart/cell_*.csv\n")
