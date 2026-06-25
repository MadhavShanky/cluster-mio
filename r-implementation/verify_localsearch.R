# Validate the SOLVER ENGINE we intend to port to C++ (no Gurobi):
#   warm-start by gradient ranking  ->  best-improvement swap local search.
# Ground truth = brute-force best subset of size lambda over the profiled objective c(S).
# If local search matches brute force ~always, the non-Gurobi port is low-risk; OA branch-and-cut
# only needs to confirm the certificate.

set.seed(7)

prof_obj <- function(Xfull, Y, cols, mu) {
  Xs <- Xfull[, cols, drop = FALSE]
  G  <- crossprod(Xs) + mu * diag(ncol(Xs))
  a  <- Y - Xs %*% solve(G, crossprod(Xs, Y))
  as.numeric(crossprod(Y, a)) / (2 * length(Y))
}

solve_localsearch <- function(Xfull, Y, fixed_cols, clust_cols, lambda, mu, n_restart = 3) {
  K <- length(clust_cols)
  obj_of <- function(sel) prof_obj(Xfull, Y, c(fixed_cols, clust_cols[sel]), mu)
  # warm start: gradient ranking at fixed-only fit
  Xs <- Xfull[, fixed_cols, drop = FALSE]
  a  <- Y - Xs %*% solve(crossprod(Xs) + mu * diag(length(fixed_cols)), crossprod(Xs, Y))
  score <- as.numeric(crossprod(Xfull[, clust_cols], a))^2
  starts <- list(order(score, decreasing = TRUE)[1:lambda])
  if (n_restart > 1) for (r in 2:n_restart) starts[[r]] <- sample(K, lambda)  # random restarts
  best_sel <- NULL; best_obj <- Inf
  for (sel0 in starts) {
    sel <- sel0; cur <- obj_of(sel)
    repeat {
      improved <- FALSE; bestswap <- NULL; bestswapobj <- cur
      out <- setdiff(seq_len(K), sel)
      for (i in sel) for (j in out) {            # swap i (in) <-> j (out)
        cand <- c(setdiff(sel, i), j)
        o <- obj_of(cand)
        if (o < bestswapobj - 1e-12) { bestswapobj <- o; bestswap <- cand; improved <- TRUE }
      }
      if (!improved) break
      sel <- bestswap; cur <- bestswapobj
    }
    if (cur < best_obj) { best_obj <- cur; best_sel <- sel }
  }
  list(sel = sort(best_sel), obj = best_obj)
}

run_one <- function(K, p_fix, nk, lambda, mu, x_corr) {
  n <- K * nk; label <- rep(1:K, each = nk)
  A <- model.matrix(~ factor(label) - 1)
  X <- matrix(rnorm(n * p_fix), n, p_fix)
  if (x_corr) X <- X + 0.8 * matrix(rnorm(K * p_fix), K, p_fix)[label, ]
  X <- cbind(1, X); pf <- ncol(X)
  gamma <- rep(0, K); active <- sample(K, lambda); gamma[active] <- rnorm(lambda, 0, 3)
  Y <- X %*% rnorm(pf) + A %*% gamma + rnorm(n)
  Xfull <- cbind(X, A); fixed_cols <- 1:pf; clust_cols <- (pf + 1):(pf + K)

  combs <- combn(K, lambda); best <- Inf; bestset <- NULL
  for (j in 1:ncol(combs)) {
    o <- prof_obj(Xfull, Y, c(fixed_cols, clust_cols[combs[, j]]), mu)
    if (o < best) { best <- o; bestset <- combs[, j] }
  }
  ls <- solve_localsearch(Xfull, Y, fixed_cols, clust_cols, lambda, mu)
  list(match = setequal(bestset, ls$sel), gap = (ls$obj - best) / abs(best))
}

bench <- function(K, p_fix, nk, lambda, mu, x_corr, reps = 50) {
  r <- replicate(reps, run_one(K, p_fix, nk, lambda, mu, x_corr), simplify = FALSE)
  cat(sprintf("K=%d pf=%d nk=%d lam=%d mu=%g corr=%s -> match %.2f  meanGap %.2e  maxGap %.2e\n",
              K, p_fix, nk, lambda, mu, x_corr,
              mean(sapply(r, `[[`, "match")), mean(sapply(r, `[[`, "gap")), max(sapply(r, `[[`, "gap"))))
}

bench(10, 3, 40, 3, 1.0, FALSE)
bench(10, 3, 40, 3, 1.0, TRUE)
bench(12, 5, 20, 4, 1.0, TRUE)
bench(14, 8, 50, 5, 1.0, TRUE)      # "High" sim setting size
bench(12, 3, 8,  4, 0.01, TRUE)     # small nk, tiny mu (paper's 1e-4-ish), strong coupling
