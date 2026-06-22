# bench_l0learn.R -- head-to-head benchmark vs L0Learn on the same
# instances, with auditable provenance. We compare, on identical cluster-selection instances:
#   SCS (our solver)  vs  L0Learn (Hazimeh-Mazumder coordinate descent + local combinatorial search)
#   vs  brute force (exact, where C(K,lambda) is enumerable).
# Fair common scoring: every method returns a support S of size lambda; ALL returned supports are
# RESCORED by the SAME unpenalized least-squares objective -- the residual of the fixed-effect-residualized
# response on the residualized cluster dummies in S, R(S) = || M y - (M D_S) gammahat_S ||^2 with
# M = I - X(X'X)^-1 X'. (SCS selects under a small ridge mu and refits unshrunk; L0Learn uses its own L2;
# rescoring on the common unpenalized objective removes that mismatch.) The winner attains the smaller
# R(S). Reports per (K,p,lambda): objective gap to brute optimum, support Jaccard, wall-clock.
#
# Gurobi NOTE: Gurobi is not re-run here (no license on this node). This script provides
# the fresh, license-free head-to-head against the closest open competitor, L0Learn.
suppressWarnings(suppressMessages({
  setwd(Sys.getenv("MIO_DIR", "C:/Users/tomch/AIProjects/research/MIO LMM"))
  source("scs.R"); source("sim/dgp.R")
  HAVE_L0 <- requireNamespace("L0Learn", quietly = TRUE)
}))
RNGkind("Mersenne-Twister")
if (!HAVE_L0) cat("WARNING: L0Learn not installed; SCS-vs-brute only.\n")

# residualized objective for a support S (vector of cluster indices, 1..K)
refit_resid <- function(My, MD, S) {
  if (length(S) == 0) return(sum(My^2))
  Ds <- MD[, S, drop = FALSE]
  fit <- stats::lm.fit(Ds, My); res <- fit$residuals
  sum(res^2)
}
brute_best <- function(My, MD, K, lambda) {            # exact best support of size lambda
  combs <- utils::combn(K, lambda); best <- Inf; bestS <- NULL
  for (j in seq_len(ncol(combs))) { S <- combs[, j]; r <- refit_resid(My, MD, S); if (r < best) { best <- r; bestS <- S } }
  list(obj = best, S = bestS)
}
jacc <- function(a, b) { a <- as.integer(a); b <- as.integer(b); u <- length(union(a, b)); if (u == 0) 1 else length(intersect(a, b)) / u }

bench_one <- function(K, p, nk, lambda, pi_active, seed) {
  d <- gen_scs(K, nk, p, regime = "sparse", pi_active = pi_active, delta = 2, rho_x = 0.5, seed = seed)
  X <- cbind(1, d$X); y <- d$y                          # intercept + covariates as the always-in block
  D <- model.matrix(~ d$cluster - 1)                    # K cluster dummies
  Mfun <- function(V) V - X %*% solve(crossprod(X), crossprod(X, V))
  My <- as.numeric(Mfun(y)); MD <- Mfun(D)

  # --- SCS ---
  t0 <- Sys.time(); fitL <- scs_fit(d$X, d$y, d$cluster, lambda = lambda, n_restart = 4L, seed = seed)
  t_scs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  S_scs <- match(fitL$selected, levels(d$cluster)); obj_scs <- refit_resid(My, MD, S_scs)

  # --- L0Learn on the residualized cluster-selection problem ---
  obj_l0 <- NA_real_; t_l0 <- NA_real_; S_l0 <- integer(0); jac_l0b <- NA_real_
  if (HAVE_L0) {
    t0 <- Sys.time()
    fitc <- tryCatch(L0Learn::L0Learn.fit(MD, My, penalty = "L0", maxSuppSize = min(K, lambda + 10L),
                                          nLambda = 200L, intercept = FALSE), error = function(e) NULL)
    t_l0 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(fitc)) {
      # (#1) score ONLY path solutions of EXACTLY size lambda (else the gap vs exact-size brute is
      # meaningless); record NA if L0Learn's path never realizes that cardinality.
      gmat <- as.matrix(coef(fitc)); supp_sizes <- apply(gmat != 0, 2, sum)
      cand <- which(supp_sizes == lambda)
      if (length(cand)) {
        best <- Inf
        for (cc in cand) { S <- which(gmat[, cc] != 0); rr <- refit_resid(My, MD, S); if (rr < best) { best <- rr; S_l0 <- S } }
        obj_l0 <- best
      }
    }
  }

  # --- brute force (only if enumerable) ---
  obj_bf <- NA_real_; S_bf <- integer(0)
  if (choose(K, lambda) <= 5e5) { bf <- brute_best(My, MD, K, lambda); obj_bf <- bf$obj; S_bf <- bf$S }

  data.frame(K = K, p = p, nk = nk, lambda = lambda, pi_active = pi_active, seed = seed,
             obj_scs = obj_scs, obj_l0learn = obj_l0, obj_brute = obj_bf,
             gap_scs = if (is.finite(obj_bf)) obj_scs - obj_bf else NA_real_,
             gap_l0learn = if (is.finite(obj_bf) && is.finite(obj_l0)) obj_l0 - obj_bf else NA_real_,
             jacc_scs_brute = if (length(S_bf)) jacc(S_scs, S_bf) else NA_real_,
             jacc_l0_brute = if (length(S_bf) && length(S_l0)) jacc(S_l0, S_bf) else NA_real_,
             jacc_scs_l0 = if (length(S_l0)) jacc(S_scs, S_l0) else NA_real_,
             t_scs = t_scs, t_l0learn = t_l0)
}

## grid: small (enumerable, includes brute) + large (SCS vs L0Learn timing/agreement only)
grid <- rbind(
  expand.grid(K = c(12L, 16L, 20L), p = 5L, nk = 30L, lambda = c(2L, 3L, 4L), pi_active = 0.15),
  expand.grid(K = c(100L, 300L),    p = 5L, nk = 30L, lambda = c(5L, 10L),    pi_active = 0.10),
  expand.grid(K = 1000L,            p = 5L, nk = 30L, lambda = 10L,            pi_active = 0.05))
REPS <- if (isTRUE(as.logical(Sys.getenv("SIM_QUICK", "FALSE")))) 2L else 20L

out <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  for (r in seq_len(REPS)) {
    out[[length(out) + 1]] <- tryCatch(bench_one(g$K, g$p, g$nk, g$lambda, g$pi_active, seed = 1000L * i + r),
                                       error = function(e) NULL)
  }
  cat(sprintf("  [%d/%d] K=%d lambda=%d done\n", i, nrow(grid), g$K, g$lambda))
}
res <- do.call(rbind, Filter(Negate(is.null), out))
dir.create("sim/results_bench", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "sim/results_bench/bench_l0learn.csv", row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), paste("L0Learn:", HAVE_L0)), "sim/results_bench/_sessionInfo.txt")
cat("done -> sim/results_bench/bench_l0learn.csv\n")
